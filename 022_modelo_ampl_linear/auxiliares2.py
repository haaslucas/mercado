from pyomo.environ import *
import numpy as np
import pandas as pd
from pyomo.util.infeasible import log_infeasible_constraints
from pyomo.core.expr.calculus.derivatives import differentiate
#from pyomo.contrib.latex_printer import latex_printer
import logging
import pandas as pd
import pickle
from pathlib import Path
import re

class _StationarityRule:
    def __init__(self, lagrangian_name, var_component_name):
        self.lagrangian_name = lagrangian_name
        self.var_component_name = var_component_name

    def __call__(self, m, *index):
        lagrangian = getattr(m, self.lagrangian_name)
        var_component = getattr(m, self.var_component_name)
        idx = index if len(index) > 1 else index[0]
        var_data = var_component[idx]
        #identify_variables of lagrangian
        deriv = get_symbolic_derivative(lagrangian, var_data)
        return deriv == 0

class _ComplementarityRule:
    def __init__(self, con_name, dual_name, modeloEmBlocos=False):
        self.con_name = con_name
        self.dual_name = dual_name
        self.modeloEmBlocos = modeloEmBlocos

    def __call__(self, m, *index):
        if self.modeloEmBlocos:
            constraint = getattr(m, self.con_name.split('.')[-1])
        else:
            constraint = getattr(m, self.con_name)
        dual_var = getattr(m, self.dual_name)

        idx = index if len(index) > 1 else index[0]
        con_data = constraint[idx]
        expr, is_equality = canonical_form(con_data)
        if is_equality:
            return Constraint.Skip

        dual_data = dual_var[idx]
        return expr * dual_data == 0

#solver_name = "conopt"  # "highs", "cbc",  "couenne", "bonmin", "ipopt", "scip", or "gcg".
#solver = SolverFactory(solver_name+"nl", executable=modules.find(solver_name), solve_io="nl")


#############################################################################################################################
#############################################################################################################################

# função auxiliar para exportar resultados do Pyomo para Excel e pickle
def _sheet_ok(name: str, used: set):
    
    '''
    Garante que o nome da planilha não ultrapasse 31 caracteres e não se repita.
    Se necessário, adiciona um sufixo numérico para evitar duplicatas.
    '''
    
    base = str(name)[:31]
    if base.lower() not in used:
        used.add(base.lower())
        return base
    i = 1
    base_trunc = base[:29]
    while f"{base_trunc}_{i}".lower() in used:
        i += 1
    new = f"{base_trunc}_{i}"
    used.add(new.lower())
    return new


def _domain_name(component):
    """
    Retorna um nome de domínio sem levantar AttributeError,
    compatível com Var, IndexedVar, Param ou IndexedParam.
    """
    # para Var escalar ou Param escalar
    try:
        return component.domain.getname(False)
    except AttributeError:
        # IndexedVar / IndexedParam
        if hasattr(component, "domain_type"):
            return component.domain_type.__name__
        return ""


def _export_pyomo_to_parquet(book: dict, path: str | Path):
    """
    Transforms and exports Pyomo model results to a Parquet file.

    The format is a wide DataFrame with hours as columns and component indices as rows.
    """
    all_dfs = []
    for comp_name, df in book.items():
        if df.empty or 'Key' not in df.columns:
            continue
        
        value_col = 'Value' if 'Value' in df.columns else 'Body'
        if value_col not in df.columns:
            continue
        
        if comp_name in ["DUALS", "OF"] or df['Key'].iloc[0] == '(scalar)':
            continue

        tipo = df['tipo'].iloc[0]
        dftemp = df.set_index('Key')[[value_col]]
        first_key = dftemp.index[0]

        df_wide = None
        if isinstance(first_key, tuple):
            if len(first_key) >= 2:
                try:
                    dftemp.index = pd.MultiIndex.from_tuples(dftemp.index)
                    if not pd.api.types.is_numeric_dtype(dftemp.index.get_level_values(-1)):
                        continue
                    
                    df_wide = dftemp[value_col].unstack(level=-1)
                    
                    if isinstance(df_wide.index, pd.MultiIndex):
                        df_wide.index = df_wide.index.map(lambda x: f"{comp_name}{x}")
                    else:
                        df_wide.index = df_wide.index.map(lambda x: f"{comp_name}[{x}]")
                except (ValueError, TypeError):
                    continue
        else:
            if pd.api.types.is_numeric_dtype(dftemp.index):
                df_wide = dftemp[value_col].to_frame().T
                df_wide.index = [comp_name]

        if df_wide is not None and not df_wide.empty:
            df_wide['tipo'] = tipo
            all_dfs.append(df_wide)
            
    if not all_dfs:
        print("Nenhuma variável indexada por tempo encontrada para exportar para Parquet.")
        return

    final_df = pd.concat(all_dfs)
    final_df.index.name = 'componente'
    
    # Move 'tipo' column to the front
    if 'tipo' in final_df.columns:
        cols = final_df.columns.tolist()
        cols.insert(0, cols.pop(cols.index('tipo')))
        final_df = final_df.loc[:, cols]

    final_df.to_parquet(path)


def export_pyomo(
    model,
    nome_arquivo_resultados: str = "resultados",
    to_xlsx: bool = True,
    to_pkl: bool = True,
    to_parquet: bool = False,
    xlsx_engine: str = "xlsxwriter",
    extra_dfs: dict = None,
):
    """
    Exporta Objetivos, Variáveis, Parâmetros, Restrições e Duais.

    Flags `to_xlsx`, `to_pkl` e `to_parquet` controlam quais arquivos serão criados.
    Sempre devolve um dicionário {sheet_name: DataFrame}.

    .. seealso:: :func:`~_sheet_ok`, :func:`~_export_pyomo_to_parquet`
    """
    if not any([to_xlsx, to_pkl, to_parquet]):
        print("Nenhuma opção de exportação de arquivo selecionada. A função 'export_pyomo' não será executada.")
        return {}
    
    xlsx_path: str | Path = "./resultados/" + nome_arquivo_resultados + ".xlsx"
    pkl_path: str | Path = "./resultados/" + nome_arquivo_resultados + ".pkl"
    parquet_path: str | Path = "./resultados/" + nome_arquivo_resultados + ".parquet"
    
    used, book = set(), {}
    dual_map = getattr(model, "dual", None)

    # ---------- Objetivos ----------
    for obj in model.component_objects(Objective, active=True):
        book[_sheet_ok(f"{obj.name}", used)] = pd.DataFrame([{
            "Nome": obj.name,
            "Sentido": "min" if obj.is_minimizing() else "max",
            "Valor": value(obj),
        }])

    # ---------- Variáveis ----------
    for var in model.component_objects(Var, active=True):
        dom_name = _domain_name(var)
        rows = []
        if var.is_indexed():
            for key, v in var.items():
                rows.append((key, v.lb, v.value, v.ub))
        else:
            v = var
            rows.append(("(scalar)", v.lb, v.value, v.ub))
        df = pd.DataFrame(rows, columns=["Key", "Lower", "Value", "Upper"])
        df["Domain"] = dom_name
        if var.name.startswith(('mu_', 'lambda_')):
            df["tipo"] = "Dual"
        else:
            df["tipo"] = "Var"
        book[_sheet_ok(var.name, used)] = df

    # ---------- Parâmetros ----------
    for par in model.component_objects(Param, active=True):
        dom_name = _domain_name(par)
        rows = []
        if par.is_indexed():
            for key in par:
                rows.append((key, value(par[key])))
        else:
            rows.append(("(scalar)", value(par)))
        df = pd.DataFrame(rows, columns=["Key", "Value"])
        df["Mutable"] = par.mutable
        df["Domain"] = dom_name
        df["tipo"] = "Param"
        book[_sheet_ok(f"{par.name}", used)] = df

    # ---------- Restrições + Duais ----------
    duals_rows = []
    for cons in model.component_objects(Constraint, active=True):
        con_name = cons.name
        if cons.is_indexed():
            items = list(cons.values())
            if not items:
                continue

            keys = [c.index() for c in items]
            lowers = [value(c.lower) for c in items]
            bodies = [value(c.body) for c in items]
            uppers = [value(c.upper) for c in items]

            if dual_map is not None:
                duals = [dual_map.get(c) for c in items]
                duals_rows.extend([(con_name, k, d) for k, d in zip(keys, duals) if d is not None])
            else:
                duals = [None] * len(items)

            df = pd.DataFrame({
                "Key": keys,
                "Lower": lowers,
                "Body": bodies,
                "Upper": uppers,
                "Dual": duals,
            })
        else:
            # Handling for scalar constraint
            c = cons
            dual = dual_map.get(c) if dual_map is not None else None
            rows = [("(scalar)", value(c.lower), value(c.body), value(c.upper), dual)]
            if dual is not None:
                duals_rows.append((con_name, "(scalar)", dual))
            df = pd.DataFrame(rows, columns=["Key", "Lower", "Body", "Upper", "Dual"])

        df["tipo"] = "Constraint"
        book[_sheet_ok(f"{con_name}", used)] = df

    if duals_rows:
        book[_sheet_ok("DUALS", used)] = pd.DataFrame(
            duals_rows, columns=["Constraint", "Key", "Dual"]
        )

    # Add extra dataframes if provided
    if extra_dfs:
        for name, df in extra_dfs.items():
            book[_sheet_ok(name, used)] = df

    import re 
    
    MAX_LEN = 31
    INVALID_CHARS_RE = re.compile(r'[:\\/?*\[\]]')

    def make_valid_sheet_name(name, used):
        s = str(name)

        # Remove/aplica regras do Excel
        s = INVALID_CHARS_RE.sub('_', s)
        s = s.strip().strip("'")  # sem espaços nas pontas e sem apóstrofo no início/fim

        if not s:
            s = "Sheet"

        # Corta para 31
        s = s[:MAX_LEN]

        # Garante unicidade adicionando sufixo _1, _2, ...
        base = s
        i = 1
        while s in used:
            suffix = f"_{i}"
            s = base[:MAX_LEN - len(suffix)] + suffix
            i += 1

        used.add(s)
        return s

    # Uso
    used_names = set()
    name_map = {}  # opcional: mapeamento original -> final

    if to_xlsx:
        with pd.ExcelWriter(xlsx_path, engine=xlsx_engine) as writer:
            for sheet, df in book.items():
                valid = make_valid_sheet_name(sheet, used_names)
                name_map[sheet] = valid
                df.to_excel(writer, sheet_name=valid, index=False)

    if to_pkl:
        with open(pkl_path, "wb") as f:
            pickle.dump(book, f, protocol=pickle.HIGHEST_PROTOCOL)

    if to_parquet:
        _export_pyomo_to_parquet(book, parquet_path)

    destinos = []
    if to_xlsx:
        destinos.append(str(xlsx_path))
    if to_pkl:
        destinos.append(str(pkl_path))
    if to_parquet:
        destinos.append(str(parquet_path))
    print("Arquivos gerados: " + ", ".join(destinos) if destinos else
          "Nenhum arquivo gerado.")
    return book

#############################################################################################################################
#############################################################################################################################

def listar_solvers():
    """
    Lista os solvers disponíveis no Pyomo e retorna um DataFrame com suas descrições.
    """
    from pyomo.opt import SolverFactory
    import numpy as np
    import pandas as pd

    # Obtém a lista de solvers disponíveis
    solver_dict = SolverFactory._doc
    solvers = pd.DataFrame.from_dict(solver_dict, orient='index', columns=['Description'])
    solvers['error'] = np.zeros(len(solvers))

    return solvers


#############################################################################################################################
#############################################################################################################################

from pyomo.common.collections import ComponentSet
from pyomo.environ import Var, Constraint, value, Piecewise

def degrees_of_freedom(m):
    # variáveis que NÃO estão fixadas
    n_var = sum(1 for v in m.component_data_objects(Var) if not v.fixed)
    # só contamos igualdades ativas (Ipopt transforma ≤ / ≥ em desigualdades, então não entram aqui)
    n_eq = sum(1 for c in m.component_data_objects(Constraint, active=True)
               if c.equality)
    return n_var - n_eq, n_var, n_eq

#dof, nvar, neq = degrees_of_freedom(m)
#print(f"DOF = {dof}  (vars={nvar},  equalities={neq})")


#############################################################################################################################
#############################################################################################################################

def canonical_form(con):
    """
    Retorna a forma canônica de uma restrição como uma tupla (expr, is_equality).
    A forma canônica é expr <= 0 ou expr == 0.
    """

    if con.lower is not None and con.upper is not None: # É igualdade?
        # igualdade: body == lower (e upper) -> body - lower == 0
        return con.body - con.lower, True
    elif con.lower is None and con.upper is not None: # É do tipo corpo <= upper?
        # corpo <= upper -> body - upper <= 0
        return con.body - con.upper, False
    elif con.lower is not None and con.upper is None: # É do tipo corpo >= lower?
        # corpo >= lower -> lower - body <= 0
        return con.lower - con.body, False
    else:
        raise ValueError(f"Restrição '{con.name}' sem limites definidos (impossível)")
    
    
    
def canonical_form_list(model):
    """
    Retorna uma lista de tuplas (expr, is_equality) para cada restrição ativa do modelo.

    .. seealso:: :func:`~canonical_form`
    """
    return [canonical_form(c) for c in model.component_data_objects(Constraint, active=True)]

def canonical_form_df(model):
    """
    Retorna um DataFrame com as restrições ativas do modelo, suas expressões canônicas e se são igualdades.

    .. seealso:: :func:`~canonical_form`
    """
    data = []
    for c in model.component_data_objects(Constraint, active=True):
        expr, is_eq = canonical_form(c)
        data.append((c.name, expr, is_eq))
    return pd.DataFrame(data, columns=["Constraint", "Expression", "IsEquality"])

#############################################################################################################################
#############################################################################################################################

def create_dual_variables(model, constraint_list=None):
    """Cria variáveis duais para um conjunto específico de restrições.

    Se `constraint_list` for fornecido, cria duais apenas para essas restrições.
    Caso contrário, o comportamento padrão é criar para todas as restrições ativas no modelo.

    Para cada restrição, cria uma variável dual correspondente.
    - Para restrições de igualdade (==), o dual é irrestrito (lambda_nome_da_restricao).
    - Para restrições de desigualdade (<=, >=), o dual é não-negativo (mu_nome_da_restricao).

    Se a restrição for indexada, a variável dual também será indexada pelo mesmo conjunto.

    Componentes criados no modelo:
    --------------------------------

    Para cada restrição `con` na lista de restrições, uma nova variável dual é adicionada:

    - **Variáveis (Var):**
        - ``lambda_{con.name}``: Criada para restrições de igualdade, com domínio `Reals`.
        - ``mu_{con.name}``: Criada para restrições de desigualdade, com domínio `NonNegativeReals`.

    Args:
        model: O modelo Pyomo.
        constraint_list (list, optional): Uma lista de componentes de restrição para
                                          os quais as variáveis duais devem ser criadas.
                                          Se None, todas as restrições ativas são usadas.

    Returns:
        Um dicionário mapeando nomes de restrições para seus nomes de variáveis duais correspondentes.
    """
    constraint_dual_map = {}
    if constraint_list is None:
        constraint_list = model.component_objects(Constraint, active=True)

    for con in constraint_list:
        # Precisamos verificar um objeto de dados para ver se é uma igualdade
        a_con_data = con if not con.is_indexed() else next(con.values())

        if a_con_data.equality:
            prefix = 'lambda_'
            domain = Reals
        else:
            prefix = 'mu_'
            domain = NonNegativeReals

        dual_var_name = prefix + con.name

        if con.is_indexed():
            index_set = con.index_set()
            dual_var = Var(index_set, domain=domain)
        else:
            dual_var = Var(domain=domain)

        setattr(model, dual_var_name, dual_var)
        constraint_dual_map[con.name] = dual_var_name

    return constraint_dual_map

#############################################################################################################################
#############################################################################################################################



#############################################################################################################################
#############################################################################################################################

from pyomo.core.expr.visitor import clone_expression   # garante cópia independente
def duplicate_scalar(model, old_con, new_name):
    """
    Cria dentro de `model` uma cópia escalar de `old_con`
    e devolve o objeto recém-criado.
    """
    # A expressão relacional completa (==, <=, >=) está em old_con.expr
    new_expr = clone_expression(old_con.expr)
    new_con  = Constraint(expr=new_expr)
    setattr(model, new_name, new_con)   # registra no modelo
    return new_con


def build_lagrangian_expression(model, constraint_dual_map, modeloEmBlocos=False):
    r"""
    Constrói a expressão da Lagrangiana de forma genérica.

    Lagr(x,λ,μ) = f(x) + Σ λ_i*h_i(x) + Σ μ_j*g_j(x)

    Args:
        model: O modelo Pyomo.
        constraint_dual_map: Um dicionário mapeando nomes de restrições
                             a nomes de variáveis duais.

    Returns:
        Uma expressão Pyomo para a Lagrangiana.

    .. seealso:: :func:`~canonical_form`
    """
    # Encontra a função objetivo ativa
    obj = next(model.component_objects(Objective, active=True), None)
    if obj is None:
        raise ValueError("O modelo não possui uma função objetivo ativa.")

    # Começa com a expressão da função objetivo
    term = obj.expr

    # Adiciona os termos das restrições
    for con_name, dual_name in constraint_dual_map.items():
        if modeloEmBlocos:
            constraint = getattr(model, con_name.split('.')[-1])
        else:
            constraint = getattr(model, con_name)   
        dual_var = getattr(model, dual_name)

        if not constraint.is_indexed():
            # Restrição escalar
            expr, _ = canonical_form(constraint)
            term += dual_var * expr
        else:
            # Restrição indexada
            for idx in constraint:
                con_data = constraint[idx]
                expr, _ = canonical_form(con_data)
                term += dual_var[idx] * expr

    return term


def add_lagrangian_derivatives_from_constraints(model, lagrangian, primal_vars):
    """Adiciona ao modelo as derivadas da Lagrangiana em relação às variáveis primais
    (condições de estacionariedade do KKT).

    Para cada variável em `primal_vars`, calcula a derivada da Lagrangiana
    e a adiciona como uma restrição de igualdade a zero.

    Componentes criados no modelo:
    --------------------------------

    Para cada variável `var` em `primal_vars`, a seguinte restrição é adicionada:

    - **Restrições (Constraint):**
        - ``dLagr_{var.name}``: A condição de estacionariedade :math:`\\frac{\\partial \\mathcal{L}}{\\partial \\text{var}} = 0`.
          O nome da restrição é gerado a partir do nome da variável, com caracteres especiais removidos.

    Args:
        model: O modelo Pyomo.
        lagrangian: A expressão da Lagrangiana.
        primal_vars: Uma lista de componentes de variáveis primais (Var) para as quais
                     as derivadas serão calculadas.
    """
    #lagrangian = build_lagrangian_expression(model, constraint_dual_map)

    # Itera sobre os componentes Var, não sobre os dados de variáveis
    lista = []
    for var_component in primal_vars:
        if var_component.is_indexed():
            for idx in var_component.index_set():
                lista.append(var_component[idx])
        else:
            lista.append(var_component)

    dLagr = differentiate(lagrangian, wrt_list=lista, mode='reverse_symbolic')

# Criar restrições definidas como igualadas a zero
    for var_index in range(len(lista)):
        deriv = dLagr[var_index]
        if deriv is not None:
            con_name = f'dLagr_{lista[var_index].name}'.replace('[', '_').replace(']', '').replace(',', '_').replace(' ', '').replace('(', '').replace(')', '')
            new_con = Constraint(expr=deriv == 0)
            setattr(model, con_name, new_con)

def get_symbolic_derivative(expr, wrt):
    """
    Calcula a derivada simbólica de uma expressão Pyomo em relação a uma variável.
    Usa o modo 'sympy' para diferenciação simbólica.

    Args:
        expr: A expressão Pyomo a ser diferenciada.
        wrt: A variável Pyomo em relação à qual derivar.

    Returns:
        A derivada da expressão como uma nova expressão Pyomo.
    """
    return differentiate(expr, wrt = wrt, mode = 'reverse_symbolic') # sympy, reverse_symbolic, reverse_numeric


def add_complementarity_slackness_conditions(model, constraint_dual_map, modeloEmBlocos=False):
    r"""Adiciona as condições de folga complementar ao modelo.

    Para cada restrição de desigualdade g(x) <= 0, adiciona a restrição
    g(x) * mu = 0, onde mu é a variável dual associada.

    Componentes criados no modelo:
    --------------------------------

    Para cada restrição de desigualdade `con` no `constraint_dual_map`, a seguinte restrição é adicionada:

    - **Restrições (Constraint):**
        - ``comp_slack_{con.name}``: A condição de folga complementar :math:`g(x) \cdot \mu = 0`.

    Args:
        model: O modelo Pyomo.
        constraint_dual_map: Um dicionário mapeando nomes de restrições
                             a nomes de variáveis duais.
        modeloEmBlocos (bool): Indica se o modelo está estruturado em blocos.

    .. seealso:: :func:`~canonical_form`
    """
    for con_name, dual_name in constraint_dual_map.items():
        if modeloEmBlocos:
            constraint = getattr(model, con_name.split('.')[-1])
        else:
            constraint = getattr(model, con_name)
        dual_var = getattr(model, dual_name)

        if not constraint.is_indexed():
            # Restrição escalar
            expr, is_equality = canonical_form(constraint)
            if not is_equality:
                comp_con_name = f'comp_slack_{con_name}'
                new_con = Constraint(expr=expr * dual_var == 0)
                setattr(model, comp_con_name, new_con)
        else:
            # Restrição indexada
            comp_con_name = f'comp_slack_{con_name}'
            index_set = constraint.index_set()

            rule = _ComplementarityRule(con_name, dual_name, modeloEmBlocos)
            new_con = Constraint(index_set, rule=rule)
            setattr(model, comp_con_name, new_con)

    
    
def add_dual_non_negativity_constraints(model, constraint_dual_map):
    """Adiciona restrições explícitas de não-negatividade para variáveis duais 'mu'.

    Para cada restrição de desigualdade, o Pyomo já define o domínio da variável
    dual associada (:math:`\mu`) como `NonNegativeReals`. Esta função adiciona
    uma restrição explícita :math:`\mu \geq 0`, o que geralmente é redundante
    mas pode ser útil para certas ferramentas de análise ou solvers.

    Componentes criados no modelo:
    --------------------------------

    Para cada variável dual `mu` correspondente a uma restrição de desigualdade:

    - **Restrições (Constraint):**
        - ``nonneg_{mu.name}``: A restrição de não-negatividade :math:`\mu \geq 0`.

    :param model: O modelo Pyomo.
    :param constraint_dual_map: Um dicionário mapeando nomes de restrições
                                a nomes de variáveis duais.
    """
    for con_name, dual_name in constraint_dual_map.items():
        if dual_name.startswith('mu_'):
            dual_var = getattr(model, dual_name)
            con_new_name = f'nonneg_{dual_name}'

            if dual_var.is_indexed():
                index_set = dual_var.index_set()
                # Usamos uma factory para capturar a variável de loop corretamente
                def rule_factory(var):
                    def rule(m, *idx):
                        return var[idx] >= 0
                    return rule
                new_con = Constraint(index_set, rule=rule_factory(dual_var))
            else:
                new_con = Constraint(expr=dual_var >= 0)
            
            setattr(model, con_new_name, new_con)




def add_complementarity_slackness_condition_linear(model, constraint_dual_map, modeloEmBlocos=False):
    r"""Adiciona as condições de folga complementar linearizadas ao modelo usando Big-M.

    Substitui a restrição de complementaridade não linear :math:`g(x) \cdot \mu = 0` por
    um conjunto de restrições lineares, utilizando uma variável binária auxiliar (:math:`y`) e
    um parâmetro Big-M (:math:`M`). As restrições adicionadas são:

    .. math::
        -g(x) \leq M \cdot (1 - y)
        \mu \leq M \cdot y

    Componentes criados no modelo:
    --------------------------------

    Para cada restrição de desigualdade `con`, os seguintes componentes são adicionados
    (com nomes baseados no nome da restrição e índice):

    - **Variáveis (Var):**
        - ``y_{con.name}_{index}``: Variável binária auxiliar para a linearização.
    - **Parâmetros (Param):**
        - ``M_{con.name}_{index}``: Parâmetro Big-M para a linearização.
    - **Restrições (Constraint):**
        - ``linear_comp1_{con.name}_{index}``: A restrição :math:`-g(x) \leq M \cdot (1 - y)`.
        - ``linear_comp2_{con.name}_{index}``: A restrição :math:`\mu \leq M \cdot y`.

    Args:
        model: O modelo Pyomo.
        constraint_dual_map: Um dicionário mapeando nomes de restrições
                             a nomes de variáveis duais.
        modeloEmBlocos (bool): Indica se o modelo está estruturado em blocos.

    .. seealso:: :func:`~canonical_form`, :func:`~add_linear_complementarity`
    """
    for con_name, dual_name in constraint_dual_map.items():
        if modeloEmBlocos:
            constraint = getattr(model, con_name.split('.')[-1])
        else:
            constraint = getattr(model, con_name)
        dual_var = getattr(model, dual_name)

        # Determinar se a restrição é escalar ou indexada
        if not constraint.is_indexed():
            # Restrição escalar
            expr, is_equality = canonical_form(constraint)
            if not is_equality:
                add_linear_complementarity(model, con_name, dual_name, expr, dual_var)
        else:
            # Restrição indexada
            index_set = constraint.index_set()
            for idx in index_set:
                con_data = constraint[idx]
                expr, is_equality = canonical_form(con_data)
                if not is_equality:
                     add_linear_complementarity(model, con_name, dual_name, expr, dual_var[idx], idx)

def add_linear_complementarity(model, con_name, dual_name, expr, dual_var, index=None):
    """
    Adiciona as restrições lineares para uma única condição de complementaridade.

    Args:
        model: O modelo Pyomo.
        con_name: Nome da restrição original.
        dual_name: Nome da variável dual.
        expr: Expressão da restrição original.
        dual_var: Variável dual.
        index: Índice, se a restrição for indexada (opcional).
    """
    identifier = f"{con_name}_{index}" if index is not None else con_name
    #remove special characters from identifier
    identifier = identifier.replace('[', '').replace(']', '').replace(', ', '').replace(',', '').replace(' ', '').replace('(', '').replace(')', '')
    # Criação da variável binária
    y_name = f'y_{identifier}'
    setattr(model, y_name, Var(within=Binary))
    y = getattr(model, y_name)

    # Criação do parâmetro M
    M_name = f'M_{identifier}'
    setattr(model, M_name, Param(initialize=1000))  # Defina um valor grande para M
    M = getattr(model, M_name)

    # Adicionar as restrições lineares
    # Restrição 1: expr <= M * (1 - y)
    constraint1_name = f'linear_comp1_{identifier}'
    constraint1 = Constraint(expr=(expr <= M * (1 - y)))
    setattr(model, constraint1_name, constraint1)

    # Restrição 2: dual_var <= M * y
    constraint2_name = f'linear_comp2_{identifier}'
    constraint2 = Constraint(expr=(dual_var <= M * y))
    setattr(model, constraint2_name, constraint2)
    
    a=1
    



def add_complementarity_slackness_condition_linear_optimized(model, constraint_dual_map, modeloEmBlocos=False):
    r"""
    Adiciona as condições de folga complementar linearizadas ao modelo usando o método Big-M.

    A condição de folga complementar, parte das condições KKT, estabelece que para cada
    restrição de desigualdade na forma canônica :math:`g(x) \leq 0`, o produto da
    expressão da restrição e sua variável dual correspondente (:math:`\mu`) deve ser zero:
    :math:`g(x) \cdot \mu = 0`. Isso implica que:
    - Se a restrição não for ativa (:math:`g(x) < 0`), sua variável dual deve ser zero (:math:`\mu = 0`).
    - Se a variável dual for positiva (:math:`\mu > 0`), a restrição deve ser ativa (:math:`g(x) = 0`).

    Esta função substitui a restrição de complementaridade não linear :math:`g(x) \cdot \mu = 0` por
    um conjunto de restrições lineares, utilizando uma variável binária auxiliar (:math:`y`) e
    um parâmetro Big-M (:math:`M`). As restrições adicionadas são:
    
    .. math::
        -g(x) \leq M \cdot (1 - y)
        \mu \leq M \cdot y

    Esta abordagem transforma um problema não linear em um Problema Linear Misto-Inteiro (MILP),
    que pode ser resolvido por solvers como Gurobi ou CPLEX.

    Esta versão é otimizada para desempenho, pré-calculando as expressões canônicas e
    criando conjuntos de índices específicos para as restrições de desigualdade,
    o que evita iteração desnecessária e a avaliação repetida de expressões.

    Componentes criados no modelo:
    --------------------------------

    Para cada restrição de desigualdade `con`, os seguintes componentes são adicionados:
    
    - **Conjuntos (Set):**
        - ``ineq_indices_{con.name}``: Um conjunto contendo apenas os índices das restrições de desigualdade,
          para otimizar a criação de componentes indexados.
    - **Variáveis (Var):**
        - ``y_{con.name}``: Variáveis binárias auxiliares indexadas pelo novo conjunto de índices.
    - **Parâmetros (Param):**
        - ``M_{con.name}``: Parâmetros Big-M, também indexados pelo novo conjunto.
    - **Restrições (Constraint):**
        - ``linear_comp1_{con.name}``: A restrição :math:`-g(x) \leq M \cdot (1 - y)`.
        - ``linear_comp2_{con.name}``: A restrição :math:`\mu \leq M \cdot y`.

    Args:
        model: O modelo Pyomo.
        constraint_dual_map: Um dicionário mapeando nomes de restrições
                             a nomes de variáveis duais.
        modeloEmBlocos (bool): Indica se o modelo está estruturado em blocos,
                             afetando como as restrições são acessadas.

    .. seealso:: :func:`~canonical_form`
    """
    for con_name, dual_name in constraint_dual_map.items():
        if modeloEmBlocos:
            constraint = getattr(model, con_name.split('.')[-1])
        else:
            constraint = getattr(model, con_name)
        
        dual_var = getattr(model, dual_name)
        base_name = con_name.replace('.', '_')
        MMM = 10000  # Defina um valor grande para M
        if not constraint.is_indexed():
            expr, is_equality = canonical_form(constraint)
            if not is_equality:
                y_name = f'y_{base_name}'
                M_param_name = f'M_{base_name}'
                
                model.add_component(y_name, Var(within=Binary))
                y = getattr(model, y_name)
                
                # SEM mutable=True, como na original
                model.add_component(M_param_name, Param(initialize=MMM))
                M = getattr(model, M_param_name)
                
                model.add_component(f'comp1_{base_name}', 
                                  Constraint(expr=-expr <= M * (1 - y)))
                model.add_component(f'comp2_{base_name}', 
                                  Constraint(expr=dual_var <= M * y))
        else:
            # MUDANÇA CRÍTICA: Usar canonical_form para filtrar, não .equality
            inequality_data = []
            for idx in constraint.index_set():
                con_data = constraint[idx]
                expr, is_equality = canonical_form(con_data)
                if not is_equality:
                    # Armazenar tanto o índice quanto a expressão já calculada
                    inequality_data.append((idx, expr))
            
            if not inequality_data:
                continue

            # Extrair apenas os índices para criar os componentes indexados
            inequality_indices_list = [idx for idx, _ in inequality_data]
            
            # Criar um Set explícito para os índices de desigualdade
            index_set_name = f'ineq_indices_{base_name}'
            # É crucial definir a dimensão correta para o conjunto
            dim = constraint.index_set().dimen
            inequality_set = Set(initialize=inequality_indices_list, dimen=dim)
            model.add_component(index_set_name, inequality_set)

            # Criar um dicionário para armazenar as expressões pré-calculadas
            expr_dict = {idx: expr for idx, expr in inequality_data}

            y_name = f'y_{base_name}'
            model.add_component(y_name, Var(inequality_set, within=Binary))
            y = getattr(model, y_name)

            M_name = f'M_{base_name}'
            # SEM mutable=True, como na original
            model.add_component(M_name, Param(inequality_set, initialize=MMM))
            M = getattr(model, M_name)

            # Usar as expressões pré-calculadas nas regras
            def linear_comp1_rule(m, *idx):
                # Usar a expressão já calculada, não recalcular
                key = idx if len(idx) > 1 else idx[0]
                expr = expr_dict[key]
                return -expr <= M[key] * (1 - y[key])

            def linear_comp2_rule(m, *idx):
                key = idx if len(idx) > 1 else idx[0]
                return dual_var[key] <= M[key] * y[key]
            
            model.add_component(f'comp1_{base_name}', 
                              Constraint(inequality_set, rule=linear_comp1_rule))
            model.add_component(f'comp2_{base_name}', 
                              Constraint(inequality_set, rule=linear_comp2_rule))




def _expr_to_latex(expr, symbol_map):
    """Converte uma expressão Pyomo para uma string LaTeX."""
    s = str(expr)
    s = re.sub(r'\w+\.', '', s)  # Remove prefixos de modelo como 'm.'
    s = s.replace('\'', '')

    # Ordena por comprimento para evitar substituições parciais
    sorted_items = sorted(symbol_map.items(), key=lambda item: len(item[0]), reverse=True)
    for pyomo_name, latex_symbol in sorted_items:
        s = re.sub(r'\b' + re.escape(pyomo_name) + r'\b', latex_symbol, s)

    # Formata índices e operadores
    s = re.sub(r'([a-zA-Z0-9_]+)\s*\[([^\]]+)\]', r'\1_{\2}', s)
    s = s.replace('*', ' \\cdot ')
    s = s.replace('==', ' = ')
    s = s.replace('<=', ' \\leq ')
    s = s.replace('>=', ' \\geq ')
    return s

def _generate_latex_table(model, symbol_map):
    def escape_latex(s):
        if s is None:
            return ''
        return str(s).replace('_', '\\_').replace('%', '\\%').replace('&', '\\&')

    params = list(model.component_objects(Param, active=True))
    variables = list(model.component_objects(Var, active=True))

    if not params and not variables:
        return ""

    table_rows = []

    # Header
    table_rows.append(r"\begin{table}[H]")
    table_rows.append(r"    \centering")
    table_rows.append(r"    \small")
    table_rows.append(fr"    \caption{{Variáveis e Parâmetros do Modelo {escape_latex(model.name)}}}")
    table_rows.append(r"    \renewcommand{\arraystretch}{1.4}")
    table_rows.append(fr"    \label{{tb:{escape_latex(model.name).replace(' ', '')}VariaveisParametros}}")
    table_rows.append(r"    \begin{tabular}{m{3cm}|m{2.5cm}m{9cm}}")
    table_rows.append(r"\hline")
    table_rows.append(r"\textbf{Tipo} & \textbf{Símbolo} & \textbf{Descrição} \\\\")
    table_rows.append(r"\hline")

    # Parameters
    if params:
        first_param = params[0]
        pyomo_name = first_param.name
        symbol = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
        doc = escape_latex(first_param.doc)
        if first_param.is_indexed():
            index_set_name = symbol_map.get(first_param.index_set().name, escape_latex(first_param.index_set().name))
            doc += f" (indexado por ${index_set_name}$)"
        
        table_rows.append(fr"\multirow{{{len(params)}}}{{*}}{{Parâmetros}}")
        table_rows.append(f"  & ${symbol}$ & {doc} \\\\")

        for p in params[1:]:
            pyomo_name = p.name
            symbol = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
            doc = escape_latex(p.doc)
            if p.is_indexed():
                index_set_name = symbol_map.get(p.index_set().name, escape_latex(p.index_set().name))
                doc += f" (indexado por ${index_set_name}$)"
            table_rows.append(f"  & ${symbol}$ & {doc} \\\\")
        table_rows.append(r"\hline")

    # Variables
    if variables:
        first_var = variables[0]
        pyomo_name = first_var.name
        symbol = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
        doc = escape_latex(first_var.doc)
        if first_var.is_indexed():
            index_set_name = symbol_map.get(first_var.index_set().name, escape_latex(first_var.index_set().name))
            doc += fr", ${symbol}_{{i}},\; i\in {index_set_name}$"
        else:
            domain = _domain_name(first_var)
            if "Reals" in domain:
                doc += fr" (${symbol} \in \mathbb{{R}}$)"

        table_rows.append(fr"\multirow{{{len(variables)}}}{{*}}{{Variáveis de decisão}}")
        table_rows.append(f"  & ${symbol}$ & {doc} \\\\")

        for v in variables[1:]:
            pyomo_name = v.name
            symbol = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
            doc = escape_latex(v.doc)
            if v.is_indexed():
                index_set_name = symbol_map.get(v.index_set().name, escape_latex(v.index_set().name))
                doc += fr", ${symbol}_{{i}},\; i\in {index_set_name}$"
            else:
                domain = _domain_name(v)
                if "Reals" in domain:
                    doc += fr" (${symbol} \in \mathbb{{R}}$)"
            table_rows.append(f"  & ${symbol}$ & {doc} \\\\")
        table_rows.append(r"\hline")

    # Footer
    table_rows.append(r"\end{tabular}")
    table_rows.append(r"\end{table}")

    return "\n".join(table_rows)


def modelo_to_latex(model, symbol_map, filename="model.tex"):
    """Gera um arquivo .tex a partir de um modelo Pyomo completo.

    A função gera seções para conjuntos, parâmetros, variáveis,
    a função objetivo e as restrições, incluindo metadados como
    documentação (doc) e unidades (units), e usando um mapa de símbolos
    para a nomenclatura em LaTeX.

    :param model: O objeto do modelo Pyomo.
    :param symbol_map: Dicionário mapeando nomes de componentes para símbolos LaTeX.
    :param filename: O nome do arquivo .tex a ser gerado (dentro da pasta 'resultados').
    :type filename: str
    
    .. seealso:: :func:`~_generate_latex_table`, :func:`~_expr_to_latex`
    """
    
    def escape_latex(s):
        if s is None:
            return ''
        return str(s).replace('_', '\\_').replace('%', '\\%').replace('&', '\\&')

    latex_parts = []

    latex_parts.append(_generate_latex_table(model, symbol_map))

    # 1. Conjuntos
    latex_parts.append("\\section*{Conjuntos}")
    for s in model.component_objects(Set, active=True):
        pyomo_name = s.name
        name = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
        #latex_parts.append(f"\\subsection*{{{name}}}")
        if s.doc:
            latex_parts.append(f"{escape_latex(s.doc)}\\\\")
        values = ", ".join(map(str, list(s)[:5]))
        if len(list(s)) > 5:
            values += ", \\dots"
        latex_parts.append(f"\\( {name} = \\{{{values}\\}} \\)")

    # 2. Parâmetros
    latex_parts.append("\\section*{Parâmetros}")
    for p in model.component_objects(Param, active=True):
        pyomo_name = p.name
        name = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
        #latex_parts.append(f"\\subsection*{{{name}}}")
        if p.doc:
            latex_parts.append(f"{escape_latex(p.doc)}\\\\")
        units = getattr(p, 'units', None)
        if units:
            latex_parts.append(f"Unidades: {escape_latex(units)}\\\\")
        if p.is_indexed():
            index_set_pyomo_name = p.index_set().name
            index_set_name = symbol_map.get(index_set_pyomo_name, escape_latex(index_set_pyomo_name))
            latex_parts.append(f"Parâmetro indexado por \\({index_set_name}\\).")
        else:
            latex_parts.append(f"\\( {name} = {value(p)} \\)")

    # 3. Variáveis
    latex_parts.append("\\section*{Variáveis}")
    for v in model.component_objects(Var, active=True):
        pyomo_name = v.name
        name = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
        #latex_parts.append(f"\\subsection*{{{name}}}")
        if v.doc:
            latex_parts.append(f"{escape_latex(v.doc)}\\\\")
        units = getattr(v, 'units', None)
        if units:
            latex_parts.append(f"Unidades: {escape_latex(units)}\\\\")
        domain = escape_latex(_domain_name(v))
        if v.is_indexed():
            index_set_pyomo_name = v.index_set().name
            index_set_name = symbol_map.get(index_set_pyomo_name, escape_latex(index_set_pyomo_name))
            latex_parts.append(f"\\( {name}_{{i}} \\) para \\( i \\in {index_set_name} \\) (Domínio: {domain})")
        else:
            latex_parts.append(f"\\( {name} \\) (Domínio: {domain})")

    # 4. Função Objetivo
    obj = next(model.component_objects(Objective, active=True), None)
    if obj is not None:
        pyomo_name = obj.name
        name = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
        latex_parts.append("\\section*{Função Objetivo}")
        #latex_parts.append(f"\\subsection*{{{name}}}")
        if obj.doc:
            latex_parts.append(f"{escape_latex(obj.doc)}\\\\")
        sense = "Maximizar" if obj.sense == maximize else "Minimizar"
        expr_str = _expr_to_latex(obj.expr, symbol_map)
        latex_parts.append(f"\\begin{{equation*}}")
        latex_parts.append(f"\\text{{{sense}}}: \\quad {expr_str}")
        latex_parts.append(f"\\end{{equation*}}")

    # 5. Restrições
    latex_parts.append("\\section*{Restrições}")
    for c in model.component_objects(Constraint, active=True):
        pyomo_name = c.name
        name = symbol_map.get(pyomo_name, escape_latex(pyomo_name))
        #latex_parts.append(f"\\subsection*{{{name}}}")
        if c.doc:
            latex_parts.append(f"{escape_latex(c.doc)}\\\\")
        if c.is_indexed():
            try:
                c_data = next(c.values())
                expr_str = _expr_to_latex(c_data.expr, symbol_map)
                idx_tuple = c_data.index()
                idx_repr = str(idx_tuple) if isinstance(idx_tuple, tuple) else str(idx_tuple)
                idx_repr = idx_repr.replace("'", "")
                latex_parts.append(f"\\begin{{equation*}}")
                latex_parts.append(f"\\text{{Exemplo para índice {idx_repr}:}} \\\\")
                latex_parts.append(f"{expr_str}")
                latex_parts.append(f"\\end{{equation*}}")
            except StopIteration:
                latex_parts.append(f"Restrição indexada mas vazia.")
        else:
            expr_str = _expr_to_latex(c.expr, symbol_map)
            latex_parts.append(f"\\begin{{equation*}}")
            latex_parts.append(expr_str)
            latex_parts.append(f"\\end{{equation*}}")

    full_latex_string = "\n".join(latex_parts)

    latex_template = fr"""\documentclass{{article}}
\usepackage{{amsmath}}
\usepackage[utf8]{{inputenc}}
\usepackage{{geometry}}
\usepackage{{longtable}}
\usepackage{{multirow}}
\usepackage{{float}}
\geometry{{a4paper, margin=1in}}
\begin{{document}}
\title{{{escape_latex(model.name)}}}
\author{{Gerado Automaticamente}}
\maketitle
{full_latex_string}
\end{{document}}
"""
    
    # Garante que o diretório de resultados exista
    Path("./resultados").mkdir(exist_ok=True)
    filepath = Path("./resultados") / filename
    
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(latex_template)
    
    print(f"Arquivo LaTeX gerado: {filepath}")



def analyze_model_gurobi(model):
    """
    Analyzes a Pyomo model using gurobi_modelanalyzer and saves a complete report.

    Args:
        model: The Pyomo model to analyze.
    """
    try:
        import gurobipy as gp
        import gurobi_modelanalyzer as gma
        import os
    except ImportError:
        print("Could not import Gurobi modules. Please install 'gurobipy' and 'gurobi-modelanalyzer'.")
        return

    filename = "temp_model_for_analyzer.mps"
    output_dir = Path("./resultados/model-analyser")

    # Create the output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n--- Iniciando análise do modelo com Gurobi Model Analyzer ---")
    try:
        # Write the model to a file, keeping variable/constraint names
        print(f"Escrevendo modelo para o arquivo temporário '{filename}'...")
        model.write(filename, io_options={"symbolic_solver_labels": True})

        # Read the model with Gurobi
        print("Lendo o modelo com GurobiPy...")
        m_gurobi = gp.read(filename)

        # Analyze the model and generate a report
        print(f"Analisando o modelo e salvando o relatório em '{output_dir}'...")
        #gma.analyze(m_gurobi, expltype="ROWS")
        
        #gma.matrix_bitmap(m_gurobi)
        gma.angle_explain(m_gurobi)
        gma.kappa_explain(m_gurobi)

        print(f"--- Análise concluíd. Relatório salvo em '{output_dir}' ---")

    except Exception as e:
        print(f"Ocorreu um erro durante a análise do modelo: {e}")
    
    #finally:
        # Clean up the temporary file
        #if os.path.exists(filename):
            #os.remove(filename)
            #print(f"Arquivo temporário '{filename}' removido.")


def add_piecewise_linear_approximation(model, P_to_linearize, P_max_func, Y, base_name, output_var_name):
    """Adiciona componentes ao modelo para representar uma aproximação linear por partes de :math:`P^2`.

    Esta função implementa a linearização da função quadrática :math:`P^2` usando uma
    aproximação linear por partes, baseada no método descrito na dissertação de
    Possagnolo (páginas 149-150). A função é projetada para funcionar com
    variáveis Pyomo indexadas.

    Para uma variável :math:`P`, a aproximação é dada por:

    .. math::
        P^2 \\approx \\sum_{y=1}^{Y} m_y \\cdot \\Delta_y

    Onde:
    - :math:`Y` é o número de segmentos lineares.
    - :math:`\\Delta_y` são as variáveis de segmento de potência.
    - :math:`m_y` é a inclinação de cada segmento.

    Componentes criados no modelo:
    --------------------------------

    Para um `base_name` fornecido, os seguintes componentes são adicionados ao `model`:

    - **Parâmetros (Param):**
        - ``delta_bar_{base_name}`` (:math:`\\bar{\\Delta}`): A largura de cada segmento da aproximação linear.
        - ``m_{base_name}`` (:math:`m_y`): As inclinações de cada segmento linear.
    - **Variáveis (Var):**
        - ``delta_{base_name}`` (:math:`\\Delta_y`): Variáveis de segmento de potência para a aproximação.
        - ``{output_var_name}``: A variável de saída que contém o valor da aproximação de :math:`P^2`.
    - **Restrições (Constraint):**
        - ``p_decomposition_{base_name}``: Decompõe :math:`P` na soma das variáveis de segmento :math:`\\Delta_y`.
        - ``p_squared_approx_constr_{base_name}``: Define a variável de saída como a soma ponderada dos segmentos.

    :param model: O modelo Pyomo.
    :type model: pyomo.environ.ConcreteModel
    :param P_to_linearize: A variável Pyomo indexada a ser linearizada (ex: `model.Pg`).
    :type P_to_linearize: pyomo.environ.Var
    :param P_max_func: Uma função que recebe um índice do `index_set` da variável e
                       retorna o valor máximo (:math:`P_{max}`) para essa variável.
    :type P_max_func: function
    :param Y: O número de segmentos lineares a serem usados na aproximação.
    :type Y: int
    :param base_name: Um prefixo a ser usado nos nomes dos componentes criados
                      para evitar colisões de nomes.
    :type base_name: str
    :param output_var_name: O nome da nova variável que conterá a aproximação linear por partes.
    :type output_var_name: str
    """
    var_indices = P_to_linearize.index_set()

    # --- Linearization Parameters ---
    delta_bar_name = f'delta_bar_{base_name}'
    model.add_component(delta_bar_name, Param(var_indices, initialize=lambda model, *idx: P_max_func(*idx) / Y))
    delta_bar = getattr(model, delta_bar_name)

    Y_set_name = f'Y_set_{base_name}'
    if not hasattr(model, Y_set_name):
        model.add_component(Y_set_name, RangeSet(1, Y))
    Y_set = getattr(model, Y_set_name)

    # Slopes for each segment (m_y = (2y - 1) * delta_bar)
    m_name = f'm_{base_name}'
    model.add_component(m_name, Param(
        var_indices, Y_set,
        initialize=lambda model, *idx: (2 * idx[-1] - 1) * delta_bar[idx[:-1]]
    ))
    m_param = getattr(model, m_name)

    # --- Variables ---
    # Power segment variables (delta)
    delta_name = f'delta_{base_name}'
    model.add_component(delta_name, Var(
        var_indices, Y_set,
        within=Reals,
        doc="Power segment variables"
    ))
    delta_var = getattr(model, delta_name)

    # New variable for the approximation of P^2
    model.add_component(output_var_name, Var(var_indices, within=NonNegativeReals))
    p_squared_approx_var = getattr(model, output_var_name)

    # --- Constraints ---
    # Decompose P into its linear segments
    def p_decomposition_rule(model, *idx):
        return P_to_linearize[idx] == sum(delta_var[idx + (y,)] for y in Y_set)
    p_decomposition_name = f'p_decomposition_{base_name}'
    model.add_component(p_decomposition_name, Constraint(var_indices, rule=p_decomposition_rule))

    # Define the new variable as the approximation expression
    def p_squared_approx_rule(model, *idx):
        return p_squared_approx_var[idx] == sum(m_param[idx + (y,)] * delta_var[idx + (y,)] for y in Y_set)
    p_squared_approx_constr_name = f'p_squared_approx_constr_{base_name}'
    model.add_component(p_squared_approx_constr_name, Constraint(var_indices, rule=p_squared_approx_rule))

    # Explicit bound constraints for delta variables for correct KKT formulation
    def delta_lower_bound_rule(model, *idx):
        return delta_var[idx] >= 0
    delta_lower_bound_name = f'delta_lower_bound_{base_name}'
    model.add_component(delta_lower_bound_name, Constraint(var_indices, Y_set, rule=delta_lower_bound_rule))

    def delta_upper_bound_rule(model, *idx):
        return delta_var[idx] <= delta_bar[idx[:-1]]
    delta_upper_bound_name = f'delta_upper_bound_{base_name}'
    model.add_component(delta_upper_bound_name, Constraint(var_indices, Y_set, rule=delta_upper_bound_rule))

# --- New function to handle piecewise linear approximation for signed variables ---
def add_piecewise_linear_approximation_signed(model, P_to_linearize, P_max_abs_func, Y, base_name, output_var_name):
    """Adiciona componentes ao modelo para uma aproximação linear por partes de :math:`P^2`,
    lidando com casos em que :math:`P` pode ser negativo.

    Esta função lineariza a expressão quadrática :math:`P^2` para uma variável :math:`P` que pode
    assumir valores positivos e negativos. A abordagem se baseia na decomposição de :math:`P`
    em suas partes positiva e negativa e, em seguida, na aplicação de uma aproximação
    linear por partes sobre o valor absoluto de :math:`P`.

    A metodologia é a seguinte:
    1. A variável :math:`P` é decomposta em duas novas variáveis não negativas, :math:`P^+` e :math:`P^-`:
       
       .. math::
          P = P^+ - P^-
          
       onde :math:`P^+ \\ge 0` e :math:`P^- \\ge 0`.
       
    2. O valor absoluto de :math:`P` é então expresso como a soma dessas duas variáveis:
       
       .. math::
          |P| = P^+ + P^-
          
    3. A função quadrática :math:`P^2` é igual a :math:`|P|^2`. A aproximação linear por partes
       é aplicada a :math:`|P|`. A variável :math:`|P|` é dividida em :math:`Y` segmentos
       lineares, cada um com uma largura :math:`\\bar{\\Delta}`.
       
       .. math::
          |P| = \\sum_{y=1}^{Y} \\Delta_y
          
       onde :math:`0 \\le \\Delta_y \\le \\bar{\\Delta}`.
       
    4. A aproximação de :math:`P^2` é então calculada como uma soma ponderada dos segmentos :math:`\\Delta_y`:
    
       .. math::
          P^2 \\approx |P|^2 \\approx \\sum_{y=1}^{Y} m_y \\cdot \\Delta_y
    
       As inclinações :math:`m_y` de cada segmento são calculadas para aproximar a curva :math:`x^2`, onde :math:`x = |P|`:
       
       .. math::
          m_y = (2y - 1) \\cdot \\bar{\\Delta}
          \\quad \\text{com} \\quad
          \\bar{\\Delta} = \\frac{P_{max}^{abs}}{Y}

    Componentes criados no modelo:
    --------------------------------

    Para um `base_name` fornecido, os seguintes componentes são adicionados ao `model`:

    - **Variáveis (Var):**
        - ``P_plus_{base_name}`` (:math:`P^+`): Variável auxiliar para a parte positiva de :math:`P`.
        - ``P_minus_{base_name}`` (:math:`P^-`): Variável auxiliar para a parte negativa de :math:`P`.
        - ``abs_P_{base_name}``: Variável auxiliar representando o valor absoluto de :math:`P`, i.e., :math:`|P|`.
        - ``delta_{base_name}`` (:math:`\\Delta_y`): Variáveis de segmento de potência para a aproximação linear.
        - ``{output_var_name}``: A variável de saída que contém o valor da aproximação de :math:`P^2`.

    - **Parâmetros (Param):**
        - ``delta_bar_{base_name}`` (:math:`\\bar{\\Delta}`): A largura de cada segmento da aproximação linear.
        - ``m_{base_name}`` (:math:`m_y`): As inclinações de cada segmento linear.

    - **Restrições (Constraint):**
        - ``p_decomposition_constraint_{base_name}``: Garante a decomposição :math:`P = P^+ - P^-`.
        - ``abs_p_constraint_{base_name}``: Define :math:`|P| = P^+ + P^-`.
        - ``abs_p_decomposition_{base_name}``: Decompõe :math:`|P|` na soma das variáveis de segmento :math:`\\Delta_y`.
        - ``p_squared_approx_constr_{base_name}``: Define a variável de saída como a soma ponderada dos segmentos.

    :param model: O modelo Pyomo.
    :type model: pyomo.environ.ConcreteModel
    :param P_to_linearize: A variável Pyomo a ser linearizada (ex: `model.Pij`).
    :type P_to_linearize: pyomo.environ.Var
    :param P_max_abs_func: Uma função que retorna o valor absoluto máximo de `P_to_linearize`.
    :type P_max_abs_func: function
    :param Y: O número de segmentos lineares a serem usados na aproximação.
    :type Y: int
    :param base_name: Um prefixo para os nomes dos componentes criados para evitar colisões de nomes.
    :type base_name: str
    :param output_var_name: O nome da nova variável que conterá a aproximação de :math:`P^2`.
    :type output_var_name: str
    """
    var_indices = P_to_linearize.index_set()

    # --- Auxiliary variables for P = P+ - P- ---
    P_plus_name = f'P_plus_{base_name}'
    P_minus_name = f'P_minus_{base_name}'
    model.add_component(P_plus_name, Var(var_indices, within=NonNegativeReals))
    model.add_component(P_minus_name, Var(var_indices, within=NonNegativeReals))
    P_plus = getattr(model, P_plus_name)
    P_minus = getattr(model, P_minus_name)

    # Constraint: P = P+ - P-
    def p_decomposition_constraint_rule(model, *idx):
        return P_to_linearize[idx] == P_plus[idx] - P_minus[idx]
    p_decomposition_constraint_name = f'p_decomposition_constraint_{base_name}'
    model.add_component(p_decomposition_constraint_name, Constraint(var_indices, rule=p_decomposition_constraint_rule))

    # --- Piecewise Linear Approximation for |P| ---
    # We approximate |P| = P+ + P-
    abs_P_name = f'abs_P_{base_name}'
    model.add_component(abs_P_name, Var(var_indices, within=NonNegativeReals))
    abs_P = getattr(model, abs_P_name)

    # Constraint: |P| = P+ + P-
    def abs_p_rule(model, *idx):
        return abs_P[idx] == P_plus[idx] + P_minus[idx]
    abs_p_name = f'abs_p_constraint_{base_name}'
    model.add_component(abs_p_name, Constraint(var_indices, rule=abs_p_rule))

    # --- Parameters for Piecewise Approximation of |P|^2 ---
    delta_bar_name = f'delta_bar_{base_name}'
    # P_max_abs_func should return the max absolute value of P_to_linearize
    model.add_component(delta_bar_name, Param(var_indices, initialize=lambda model, *idx: P_max_abs_func(*idx) / Y))
    delta_bar = getattr(model, delta_bar_name)

    Y_set_name = f'Y_set_{base_name}'
    if not hasattr(model, Y_set_name):
        model.add_component(Y_set_name, RangeSet(1, Y))
    Y_set = getattr(model, Y_set_name)

    # Slopes for each segment of x^2, where x = |P|
    # m_y = (2y - 1) * delta_bar
    m_name = f'm_{base_name}'
    model.add_component(m_name, Param(
        var_indices, Y_set,
        initialize=lambda model, *idx: (2 * idx[-1] - 1) * delta_bar[idx[:-1]]
    ))
    m_param = getattr(model, m_name)

    # --- Variables for Piecewise Approximation ---
    # Segment variables for |P| (delta)
    delta_name = f'delta_{base_name}'
    model.add_component(delta_name, Var(
        var_indices, Y_set,
        within=Reals,
        doc="Piecewise segment variables for |P|"
    ))
    delta_var = getattr(model, delta_name)

    # New variable for the approximation of P^2 = |P|^2
    model.add_component(output_var_name, Var(var_indices, within=NonNegativeReals))
    p_squared_approx_var = getattr(model, output_var_name)

    # --- Constraints for Piecewise Approximation ---
    # Decompose |P| into its linear segments
    def abs_p_decomposition_rule(model, *idx):
        return abs_P[idx] == sum(delta_var[idx + (y,)] for y in Y_set)
    abs_p_decomposition_name = f'abs_p_decomposition_{base_name}'
    model.add_component(abs_p_decomposition_name, Constraint(var_indices, rule=abs_p_decomposition_rule))

    # Define the new variable as the approximation expression for |P|^2
    def p_squared_approx_rule(model, *idx):
        return p_squared_approx_var[idx] == sum(m_param[idx + (y,)] * delta_var[idx + (y,)] for y in Y_set)
    p_squared_approx_constr_name = f'p_squared_approx_constr_{base_name}'
    model.add_component(p_squared_approx_constr_name, Constraint(var_indices, rule=p_squared_approx_rule))

    # Explicit bound constraints for delta variables for correct KKT formulation
    def delta_lower_bound_rule(model, *idx):
        return delta_var[idx] >= 0
    delta_lower_bound_name = f'delta_lower_bound_{base_name}'
    model.add_component(delta_lower_bound_name, Constraint(var_indices, Y_set, rule=delta_lower_bound_rule))

    def delta_upper_bound_rule(model, *idx):
        return delta_var[idx] <= delta_bar[idx[:-1]]
    delta_upper_bound_name = f'delta_upper_bound_{base_name}'
    model.add_component(delta_upper_bound_name, Constraint(var_indices, Y_set, rule=delta_upper_bound_rule))

from pyomo.environ import *

def fix_binaries(model: ConcreteModel) -> ConcreteModel:
    """
    Itera sobre todas as variáveis de um modelo Pyomo resolvido, identifica
    as que são binárias e as fixa em seu valor da solução atual.

    Esta função é o primeiro passo da técnica "fix-and-solve" para extrair
    valores duais de um problema MILP.

    Args:
        model (pyomo.ConcreteModel): Um modelo Pyomo que já foi resolvido
                                     com sucesso e contém os valores da solução.

    Returns:
        pyomo.ConcreteModel: O mesmo objeto de modelo, mas com as variáveis
                             binárias agora fixadas, pronto para uma segunda resolução.
    """
    fixed_vars_count = 0
    print(f"--- Iniciando a rotina para fixar variáveis binárias no modelo '{model.name}' ---")

    # Itera sobre todos os componentes do tipo Var (variável) no modelo
    for var_object in model.component_objects(Var, active=True):
        # Itera sobre cada índice da variável (funciona para vars únicas ou indexadas)
        for index in var_object:
            variable = var_object[index]

            # A verificação chave: o domínio da variável é binário?
            if variable.is_binary():
                try:
                    # Pega o valor da solução atual (e.g., 0.0 ou 1.0)
                    solution_value = value(variable)
                    
                    # Fixa a variável neste valor.
                    # A partir daqui, o solver a tratará como uma constante.
                    variable.fix(solution_value)
                    
                    print(f"  -> Variável '{variable.name}' fixada em: {solution_value}")
                    fixed_vars_count += 1
                
                except (ValueError, TypeError):
                    # Medida de segurança caso a variável não tenha um valor na solução
                    print(f"  -> AVISO: Não foi possível fixar '{variable.name}'. Valor da solução não encontrado.")
    
    print(f"--- Processo concluído. Total de {fixed_vars_count} variáveis binárias fixadas. ---")
    
    return model


def add_piecewise_tangent_approximation(model, P_to_linearize, P_max_func, Y, base_name, output_var_name):
    """
    Adiciona componentes para uma aproximação PWL de P^2 usando o método do envelope de tangentes. (v2 - Corrigido)
    """
    var_indices = P_to_linearize.index_set()

    # --- Parâmetros da Aproximação ---
    
    Y_set_name = f'Y_set_{base_name}'
    if not hasattr(model, Y_set_name):
        model.add_component(Y_set_name, RangeSet(0, Y))
    Y_set = getattr(model, Y_set_name)

    tangent_points_name = f'tangent_points_{base_name}'
    def tangent_points_rule(model, *args): # <<< CORREÇÃO AQUI
        # O último argumento é do Y_set, o resto é de var_indices
        k = args[-1]
        var_idx = args[:-1]
        p_max = P_max_func(*var_idx)
        return p_max * k / Y
    model.add_component(tangent_points_name, Param(var_indices, Y_set, initialize=tangent_points_rule))
    tangent_points = getattr(model, tangent_points_name)

    slope_name = f'slope_{base_name}'
    def slope_rule(model, *args): # <<< CORREÇÃO AQUI
        return 2 * tangent_points[args]
    model.add_component(slope_name, Param(var_indices, Y_set, initialize=slope_rule))
    m_param = getattr(model, slope_name)

    intercept_name = f'intercept_{base_name}'
    def intercept_rule(model, *args): # <<< CORREÇÃO AQUI
        pk = tangent_points[args]
        return -pk**2
    model.add_component(intercept_name, Param(var_indices, Y_set, initialize=intercept_rule))
    c_param = getattr(model, intercept_name)

    # --- Variável de Saída ---
    
    if not hasattr(model, output_var_name):
        model.add_component(output_var_name, Var(var_indices, within=NonNegativeReals))
    p_squared_approx_var = getattr(model, output_var_name)

    # --- Restrição do Envelope ---
    
    def pw_tangent_rule(model, *args): # <<< CORREÇÃO AQUI
        # O último argumento é do Y_set (k), o resto é de var_indices (idx)
        k = args[-1]
        var_idx = args[:-1]
        
        P = P_to_linearize[var_idx]
        P_sq_approx = p_squared_approx_var[var_idx]
        m = m_param[args] # Acessa com o tuple completo (idx..., k)
        c = c_param[args]
        return P_sq_approx >= m * P + c
        
    pw_tangent_constr_name = f'pw_tangent_constr_{base_name}'
    model.add_component(pw_tangent_constr_name, Constraint(var_indices, Y_set, rule=pw_tangent_rule))

# Assumimos que a função 'add_piecewise_tangent_approximation' está disponível no mesmo escopo.

def add_piecewise_tangent_approximation_signed(model, P_to_linearize, P_max_abs_func, Y, base_name, output_var_name):
    """
    Cria uma aproximação PWL para P^2 para variáveis com sinal, usando um método KKT-robusto.

    A função decompõe P em P+ e P-, define |P| = P+ + P-, e então aplica
    a aproximação por envelope de tangentes em |P| para calcular |P|^2, que é igual a P^2.
    """
    var_indices = P_to_linearize.index_set()

    # --- Decomposição da variável com sinal ---

    # 1. Variáveis auxiliares para P = P+ - P-
    P_plus_name = f'P_plus_{base_name}'
    P_minus_name = f'P_minus_{base_name}'
    model.add_component(P_plus_name, Var(var_indices, within=NonNegativeReals))
    model.add_component(P_minus_name, Var(var_indices, within=NonNegativeReals))
    P_plus = getattr(model, P_plus_name)
    P_minus = getattr(model, P_minus_name)

    # 2. Restrição de decomposição: P = P+ - P-
    def p_decomposition_rule(model, *idx):
        return P_to_linearize[idx] == P_plus[idx] - P_minus[idx]
    p_decomposition_name = f'p_decomposition_constraint_{base_name}'
    model.add_component(p_decomposition_name, Constraint(var_indices, rule=p_decomposition_rule))

    # --- Aproximação de |P|^2 ---

    # 3. Variável para o valor absoluto: |P|
    # Esta variável será o alvo da nossa linearização.
    abs_P_name = f'abs_P_{base_name}'
    model.add_component(abs_P_name, Var(var_indices, within=NonNegativeReals))
    abs_P = getattr(model, abs_P_name)

    # 4. Restrição de definição do valor absoluto: |P| = P+ + P-
    def abs_p_rule(model, *idx):
        return abs_P[idx] == P_plus[idx] + P_minus[idx]
    abs_p_name = f'abs_p_constraint_{base_name}'
    model.add_component(abs_p_name, Constraint(var_indices, rule=abs_p_rule))

    # 5. Aplicar a aproximação de tangentes em |P| para obter |P|^2
    # A função que já corrigimos faz todo o trabalho aqui.
    # Note que a variável de entrada é 'abs_P' e o P_max é o mesmo 'P_max_abs'.
    add_piecewise_tangent_approximation(
        model=model,
        P_to_linearize=abs_P,  # Linearizamos |P|
        P_max_func=P_max_abs_func, # O limite de |P|
        Y=Y,
        # Usamos um 'internal_base_name' para evitar colisão de nomes de componentes
        base_name=f'{base_name}_abs',
        output_var_name=output_var_name # O resultado final é P^2
    )
