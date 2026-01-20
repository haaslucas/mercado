from pyomo.environ import *
import pandas as pd
import numpy as np
import sys
sys.path.append('.')
from src.auxiliares_kkt import *


def agregator_bilevel(
        prelog=False,
        poslog=False,
        prelog_simples=False,
        analyze_model=False,
        PW=5, #[5,20]
        pw_type='SOS2', # ['INC', 'SOS2', 'BIGM_BIN', 'BIGM_SOS1', 'DCC', 'DLOG', 'MC', ]
        pw_constr_type_I='EQ', #['EQ','LB']
        pw_constr_type_V='EQ',
        mccormick=False,
        partes=1, # [1,2,5,10]
        name = 'agregator3',
        horas=24,
        Y = 10, # anos 
        C_REF=30000/365, #$/MW
        ALPHA=0.85,
        TAXA_MINIMA_ATRATIVIDADE = 0.12, # Taxa para decisão se vale a pena investir no BESS
        MAX_CAP = 1, # capacidade máxima de cada BESS
        W_MAX = 1, #quantidade máxima de BESS
        SBASE = 10, # MVA
        agregator = 'agregator',
        price = 'pi',
        big_m_value=100,
        load_factor = 10,
        fator_sol = 1
    ):

    #horas = 8
    m = ConcreteModel(name='Agregator BESS Bilevel',
            doc='Bilevel model for BESS aggregator profit maximization')

    m.dual = Suffix(direction=Suffix.IMPORT)

    SBASE = SBASE  # Base power in MVA, 
    VBASE = 24.9  # Base voltage in kV
    IBASE = SBASE * 1e3 / VBASE  # Base current in A
    fator = 100
    C_REF = C_REF/SBASE # Custo de referência em $/pu
    ############ PROCESSAMENTO DE DADOS ############
    PASTA = './modelos/data/'
    DLINES  = pd.read_csv(PASTA + 'distribution_line_data.csv') #obs: os valores em PU utilizam SBASE = 1 MVA e VBASE = 24.9 kV
    DLINES = DLINES.set_index(['From', 'To'])

    DLINES['R (pu)'] = DLINES['R (pu)'] / fator * SBASE # resistência
    DLINES['X (pu)'] = DLINES['X (pu)'] / fator * SBASE # reatância

    load = [0.035, 0.030, 0.028, 0.027, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090, 0.095, 0.100, 0.090, 0.080, 0.060, 0.040]
    load = pd.Series(load, index=range(1,25))*load_factor # em MW
    load = load / SBASE  # em pu.
    load_q = load * 0.3287 # PF=0.95
    
    sol = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0, 0.9, 0.8, 0.6, 0.4, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0]
    sol = pd.Series(sol, index=range(1,25)) * fator_sol
    Pnom_pv = {13:1, 23: 0.8, 33: 1.4} 
    wsp = [2, 3, 4, 4, 10, 10, 4, 5, 6, 7, 7, 8, 8, 9, 8, 7, 8, 9, 12, 13, 14, 12, 10,  9]
    wholesale_prices = pd.Series(wsp, index=range(1,25)) # preços de atacado fictícios para 24 horas
    # Create multi-year data
    load_multi_year = {}
    load_q_multi_year = {}
    wsp_multi_year = {}
    load_growth_rate = 1.02 # 2% ao ano
    wsp_growth_rate = 1.05 # 1% ao ano
    for year in range(1, Y + 1):
        for hour in range(1, 25):
            if hour <= horas:
                load_multi_year[(hour, year)] = load.loc[hour] * (load_growth_rate ** (year - 1))
                load_q_multi_year[(hour, year)] = load_q.loc[hour] * (load_growth_rate ** (year - 1))
                wsp_multi_year[(hour, year)] = wholesale_prices.loc[hour] * (wsp_growth_rate ** (year - 1))
    ####################################################
    
    #################### Conjuntos #####################


    m.i = RangeSet(1, 34, name='i', doc='Nós da rede de distribuição')
    m.t = RangeSet(1, horas, name='t', doc='Horas do dia')
    m.y = RangeSet(1, Y, name='y', doc='Anos do horizonte de planejamento')
    m.s = RangeSet(1, 1, name='s', doc='Subestacoes')
    
    #m.j = Set(initialize=[30], name='j', doc='Conjunto de nós candidatos para recursos distribuídos')
    m.pv = Set(initialize=Pnom_pv.keys(), name='pv', doc='Conjunto de nós com geração fotovoltaica')
    m.j = Set(initialize=m.i, name='j', doc='Conjunto de nós candidatos para recursos distribuídos')
    m.l = Set(initialize=[(int(row.name[0]), int(row.name[1])) for _, row in DLINES.iterrows()], 
              name='l', doc='Linhas de distribuição')

    #################### Parametros ####################
    m.Vmin = Param(initialize=0.90, doc='Tensão mínima nos nós da rede de distribuição [pu]')
    m.Vmax = Param(initialize=1.05, doc='Tensão máxima nos nós da rede de distribuição [pu]')
    m.eta_c = Param(initialize=0.95, doc='Eficiência de carga do BESS')
    m.eta_d = Param(initialize=0.95, doc='Eficiência de descarga do BESS')
    m.slack = Param(initialize=1, name='slack', doc='Nó de referência da tensão (slack bus)')
    m.W_max = Param(initialize=W_MAX, name='W_max', doc='Número máximo de DGs a serem instalados')
    m.R     = Param(m.l, initialize=lambda model, i, j: DLINES.loc[(i, j), 'R (pu)'],
                    name='R', doc='Resistência das linhas de distribuição')
    m.X     = Param(m.l, initialize=lambda model, i, j: DLINES.loc[(i, j), 'X (pu)'],
                    name='X', doc='Reatância das linhas de distribuição')
    m.rho   = Param(m.t, m.y, initialize=lambda model, t, y: wsp_multi_year[(t, y)],
                    name='rho', doc='Preço de mercado de energia no atacado')  
    m.Load  = Param(m.i, m.t, m.y, initialize=lambda model, i, t, y: load_multi_year[(t, y)],
                    name='Load', doc='Carga ativa nos nós da rede de distribuição')
    m.Qd    = Param(m.i, m.t, m.y, initialize=lambda model, i, t, y: load_q_multi_year[(t, y)],
                    name='Qd', doc='Carga reativa nos nós da rede de distribuição')

    m.Ppv   = Var(m.pv, m.t, m.y, initialize=lambda m, i, t, y: sol.loc[t] * Pnom_pv[i] / SBASE,
                    name='Ppv', doc='Potência de geração fotovoltaica')
    
    m.Y = RangeSet(1, Y, name='V', doc='Vida útil do BESS em anos')

    m.c_ref  = Param(initialize=C_REF, name='C_REF', doc='Custo de referência para investimento em BESS [$]')
    m.alpha  = Param(initialize=ALPHA, name='ALPHA', doc='Fator de escala para custo de investimento em BESS')
    m.tma    = Param(initialize=TAXA_MINIMA_ATRATIVIDADE, name='TAXA_MINIMA_ATRATIVIDADE',
                    doc='Taxa mínima de atratividade para investimento no BESS')    
    
    opcoes_capacidade = [0.2*MAX_CAP, 0.4*MAX_CAP, 0.6*MAX_CAP, 0.8*MAX_CAP, 1.0*MAX_CAP]

    dict_tamanhos = {}
    dict_custos = {}
    for k, cap_mw in enumerate(opcoes_capacidade, 1):
        tam_pu = cap_mw / SBASE
        custo = C_REF * ((tam_pu / 1.0) ** m.alpha)
        dict_tamanhos[k] = tam_pu
        dict_custos[k] = custo

    m.K_inv = Set(initialize=dict_tamanhos.keys(), doc='Índice das opções de capacidade do BESS')
    m.tab_cap = Param(m.K_inv, initialize=dict_tamanhos, doc='Capacidade de cada opção [pu*h]')
    m.tab_custo = Param(m.K_inv, initialize=dict_custos, doc='Custo de investimento de cada opção [$]')

    m.SOC   = Var(m.j, m.t, m.y, name='SOC', doc='Estado de carga do BESS')
    m.P_ch  = Var(m.j, m.t, m.y, name='P_ch', doc='Potência de carga do BESS')
    m.P_dch = Var(m.j, m.t, m.y, name='P_dch', doc='Potência de descarga do BESS')
    m.Qbess = Var(m.j, m.t, m.y, name='Qbess', doc='Potência reativa do BESS')
    m.Pse   = Var(m.s, m.t, m.y, name='Pse',
                doc='Potência ativa de energia comprada do mercado atacadista')
    m.Qse   = Var(m.s, m.t, m.y, name='Qse',
                doc='Potência reativa de energia comprada do mercado atacadista')
    m.e     = Var(m.i, m.t, m.y, name='e', 
                doc='Tensão ao quadrado nos nós da rede de distribuição')
    
    #m.u_exc = Var(m.j, m.t, m.y, within=Binary, name='u_exc', 
    #            doc='Variável binária para garantir exclusividade de carga/descarga')
    m.Pij   = Var(m.l, m.t, m.y, name='Pij', doc='Fluxo de potência ativa na linha')
    m.Qij   = Var(m.l, m.t, m.y, name='Qij', doc='Fluxo de potência reativa na linha')
    m.w     = Var(m.j, m.y, m.K_inv, within=Binary, name='w', doc='Variável binária de investimento em BESS') # análogo a eq. 1b

    #def z_rule(m, j, y):
    #    current_year_int = y
    #    return sum(m.w[j, y_prime, k] for y_prime in m.y if y_prime <= current_year_int for k in m.K_inv)
    #m.z = Expression(m.j, m.y, rule=z_rule, name='z', doc='Variável binária indicando se BESS está operacional')
    
    def rule_definir_cap(m, j, y):
        return sum(m.w[j, y_prime, k] * m.tab_cap[k] for y_prime in m.y if y_prime <= y for k in m.K_inv)
    m.Cap = Expression(m.j, m.y, rule=rule_definir_cap,
        doc='Definição da capacidade do BESS no nó j', name='Definição de Capacidade')

    #m.PchExc    = Constraint(m.j, m.t, m.y, name='P_ch_EXCLUSIVE',
    #    doc='Garante que P_ch é zero se não estiver carregando',
    #    rule=lambda m, j, t, y: m.P_ch[j, t, y] <= (MAX_CAP/SBASE) * m.u_exc[j, t, y] )

    #m.PdchExc   = Constraint(m.j, m.t, m.y, name='P_dch_EXCLUSIVE',
    #    doc='Garante que P_dch é zero se não estiver descarregando',
    #    rule=lambda m, j, t, y: m.P_dch[j, t, y] <= (MAX_CAP/SBASE)  * (1 - m.u_exc[j, t, y]) )

    m.e_LWBND = Constraint(m.i, m.t, m.y, name='e_LWBND', # eq. 1g
        rule=lambda m, i, t, y: m.e[i, t, y] >= m.Vmin**2)
    m.e_UPBND = Constraint(m.i, m.t, m.y, name='e_UPBND', # eq. 1g
        rule=lambda m, i, t, y: m.e[i, t, y] <= m.Vmax**2)
    
    m.SOCLwBnd = Constraint(m.j, m.t, m.y, name='SOCLwBnd', doc='Lower bound for SOC',
        rule=lambda m, j, t, y: m.SOC[j, t, y] >= 0.1 * m.Cap[j, y] )
    
    m.SOCUpBnd = Constraint(m.j, m.t, m.y, name='SOCUpBnd', doc='Upper bound for SOC',
        rule=lambda m, j, t, y: m.SOC[j, t, y] <= m.Cap[j, y])
    
    m.SOCfinal = Constraint(m.j, m.y, name='SOC0= SOCfinal',
        doc='Garante que o SOC final é igual ao mínimo',
        rule=lambda m, j, y: m.SOC[j, horas, y] == 0.1 * m.Cap[j, y])
    
    m.PchUpBnd = Constraint(m.j, m.t, m.y, name='PchUpBound', # análogo a eq. 1c
        doc='Upper bound for P_ch',
        rule=lambda m, j, t, y: m.P_ch[j, t, y] <= m.Cap[j, y])

    m.PdchUpBnd = Constraint(m.j, m.t, m.y, name='PdchUpBnd', # análogo a eq. 1c
        doc='Upper bound for P_dch',
        rule=lambda m, j, t, y: m.P_dch[j, t, y] <= m.Cap[j, y])

    m.PchLwBnd = Constraint(m.j, m.t, m.y, name='PchLwBnd',
        doc='Lower bound for P_ch',
        rule=lambda m, j, t, y: m.P_ch[j, t, y] >= 0)
    
    m.PdchLwBnd = Constraint(m.j, m.t, m.y, name='PdchLwBnd',
        doc='Lower bound for P_dch',
        rule=lambda m, j, t, y: m.P_dch[j, t, y] >= 0)
    
    #m.PchExc    = Constraint(m.j, m.t, m.y, name='P_ch_EXCLUSIVE',
    #    doc='Garante que P_ch é zero se não estiver carregando',
    #    rule=lambda m, j, t, y: m.P_ch[j, t, y] <= (MAX_CAP/SBASE)) # * m.u_exc[j, t, y]

    #m.PdchExc   = Constraint(m.j, m.t, m.y, name='P_dch_EXCLUSIVE',
    #    doc='Garante que P_dch é zero se não estiver descarregando',
    #    rule=lambda m, j, t, y: m.P_dch[j, t, y] <= (MAX_CAP/SBASE)) #  * (1 - m.u_exc[j, t, y])

    m.eSlack = Constraint(m.t, m.y, name='e_SLACK',
        doc='Slack bus voltage constraint',
        rule=lambda model, t, y: model.e[m.slack, t, y] == 1)

    m.PseUpBnd = Constraint(m.s, m.t, m.y, name='PSE_BOUNDS_UPPER',
        doc='Upper bound for energy purchased from the wholesale market',
        rule=lambda model, s, t, y: model.Pse[s, t, y] <= 50/SBASE)
    
    m.PseLwBnd = Constraint(m.s, m.t, m.y, name='PSE_BOUNDS_LOWER',
        doc='Lower bound for energy purchased from the wholesale market',
        rule=lambda model, s, t, y: model.Pse[s, t, y] >= 0)

    m.QseUpBnd = Constraint(m.s, m.t, m.y, name='QSE_BOUNDS_UPPER',
        doc='Upper bound for reactive energy purchased from the wholesale market',
        rule=lambda model, s, t, y: model.Qse[s, t, y] <= 50/SBASE)
    
    m.QseLwBnd = Constraint(m.s, m.t, m.y, name='QSE_BOUNDS_LOWER',
        doc='Lower bound for reactive energy purchased from the wholesale market',
        rule=lambda model, s, t, y: model.Qse[s, t, y] >= 0)
        
    m.QbessUpBnd = Constraint(m.j, m.t, m.y, name='QbessUpBound', # análogo a eq. 1d
        doc='Upper bound for Qbess',
        rule=lambda m, j, t, y: m.Qbess[j, t, y] <= m.Cap[j, y])

    m.QbessLwBnd = Constraint(m.j, m.t, m.y, name='QbessLwBnd',
        doc='Lower bound for Qbess',
        rule=lambda m, j, t, y: m.Qbess[j, t, y] >= 0)

    def balanco_pot_ativa(m, i, t, y):
        pg = 0
        if i in m.s: pg += m.Pse[i, t, y]
        if i in m.j: pg += m.P_dch[i, t, y] - m.P_ch[i, t, y]
        if i in m.pv: pg += m.Ppv[i, t, y]
        
        outflow = sum(m.Pij[k, j, t, y] for k,j in m.l if k == i)
        inflow  = sum(m.Pij[k, j, t, y] for k,j in m.l if j == i)
        
        return -(pg - m.Load[i, t, y]) == -(outflow - inflow)
    m.BALANCO_POT_ATIVA = Constraint(m.i, m.t, m.y, rule=balanco_pot_ativa, name='balanco_pot_ativa')

    def balanco_pot_reativa(m, i, t, y):
        qg = 0
        if i in m.s: qg += m.Qse[i, t, y]
        if i in m.j: qg += m.Qbess[i, t, y]
        
        outflow = sum(m.Qij[k, j, t, y] for k,j in m.l if k == i)
        inflow  = sum(m.Qij[k, j, t, y] for k,j in m.l if j == i)
        
        return qg - m.Qd[i, t, y] == outflow - inflow
    m.BALANCO_POT_REATIVA = Constraint(m.i, m.t, m.y, rule=balanco_pot_reativa, name='balanco_pot_reativa')

    def voltage_drop_rule(m, i, j, t, y):
        return m.e[j, t, y] == m.e[i, t, y] - 2 * (m.R[i,j] * m.Pij[i,j,t,y] + m.X[i,j] * m.Qij[i,j,t,y])
    m.VOLTAGE_DROP = Constraint(m.l, m.t, m.y, rule=voltage_drop_rule, name='voltage_drop')

    def soc_balance_rule(m, j, t, y):
        if t == 1:
            return m.SOC[j, t, y] == m.Cap[j, y]*0.1 + (m.eta_c * m.P_ch[j,t,y] - (1/m.eta_d) * m.P_dch[j,t,y])
        return m.SOC[j, t, y] == m.SOC[j, t-1, y] + (m.eta_c * m.P_ch[j,t,y] - (1/m.eta_d) * m.P_dch[j,t,y])
    m.SOC_BALANCE = Constraint(m.j, m.t, m.y, rule=soc_balance_rule, name='SOC_BALANCE', doc='Balanço de energia do BESS')
    
    LISTA_VAR_DUAIS = ['Cap','w','u_exc','uIijComp']
    LISTA_PRIMAIS = [v for v in m.component_objects(Var, active=True) if v.name not in  LISTA_VAR_DUAIS]
    #LISTA_PRIMAIS = [v for v in LISTA_PRIMAIS if 'pw' not in v.name]
    LISTA_CONSTRS_DUAIS = ['PchExc','PdchExc','IijPlusComp','IijMinusComp','CapInvestLink']
    LISTA_CONSTRS = [c for c in m.component_objects(Constraint, active=True) if c.name not in LISTA_CONSTRS_DUAIS]
    #LISTA_CONSTRS = [c for c in LISTA_CONSTRS if 'pw' not in c.name]

    DIC_DUAIS = create_dual_variables(m, constraint_list=LISTA_CONSTRS)
    
    m.Pbess = Expression(m.j, m.t, m.y,
        rule=lambda model, j, t, y: model.P_dch[j, t, y] - model.P_ch[j, t, y],
        
        doc='Potência de descarga do BESS para o produto não-linear com o LMP')
    
    if price == 'rho':
        m.pi = Param(m.i, m.t, m.y, initialize=lambda m, i, t, y: m.rho[t, y], name='pi',
            doc='Preços nodais de energia')

    if price == 'pi_rho':
        m.pi = Reference(m.lambda_BALANCO_POT_ATIVA) 
        m.pi_rule = Constraint(m.i, m.t, m.y, name='pi_rule',
            doc='Definição dos preços nodais de energia',
            rule=lambda m, i, t, y: m.pi[i, t, y] == m.rho[t, y])

    if price == 'pi':
        m.pi = Reference(m.lambda_BALANCO_POT_ATIVA)
        #m.pi.store_values({(i, t, y): m.rho[i, t, y] for i in m.I for t in m.T for y in m.Y})
        #create a 
        #m.LmpUpBnd = Constraint(m.i, m.t, m.y, name = 'Up Bound do LMP para linearizacao',
        #    rule=lambda m, i, t, y: m.pi[i, t, y] <= m.rho[t, y] * 1.5) #1.1
        #m.LmpLwBnd = Constraint(m.i, m.t, m.y, name = 'Lw Bound do LMP para linearizacao', 
        #    rule=lambda m, i, t, y: m.pi[i, t, y] >= m.rho[t, y] * 0.2) # 0.9

    m.custo_agregador = Expression(expr=sum(
        (sum( m.pi[j, t, y]* m.Pbess[j, t, y] for t in m.t ) / ( (1 + m.tma) ** (y - 1) ) for j in m.j for y in m.y)),
        doc='Custo operacional do agregador de BESSs no sistema de distribuição'
    )

    m.custos_pse = Expression(expr=sum(
        (sum( m.rho[t, y] * m.Pse[s, t, y] for t in m.t ) / ( (1 + m.tma) ** (y - 1) ) for s in m.s for y in m.y)),
        doc='Custo operacional dos outros agentes no sistema de distribuição relacionado a Pse')
    
    m.custos_ppv = Expression(expr=sum(
        (sum( m.rho[t, y] * m.Ppv[j, t, y] for t in m.t ) / ( (1 + m.tma) ** (y - 1) ) for j in m.pv for y in m.y)),
        doc='Custo operacional dos outros agentes no sistema de distribuição relacionado a Ppv')

    m.custo_outros = Expression(expr=
        m.custos_pse + m.custos_ppv,
        doc='Custo operacional dos outros agentes no sistema de distribuição'
    )


    m.objinferior_expr= Expression(expr=
             m.custo_outros + m.custo_agregador,
            doc='OBJ do problema inferior (agregador + outros agentes)')
    
    
    m.dual_obj_expr = build_dual_objective_expression(m, DIC_DUAIS)
    m.obj_dual = Expression(expr= m.dual_obj_expr,   
            doc='OBJ do problema inferior em termos dos duais')
        
    # ISSUE: o sinal em -m.obj_dual é pq o create_dual_variables é do KKT e não do Dual.
    m.custo_agregador_linear = Expression(expr=-m.obj_dual - m.custo_outros,
                                doc='Custo do agregador em valor presente (expressão linearizada)')

    #m.strong_duality = Constraint(doc='Strong duality constraint', # eq. 3e
    #    rule=lambda model: m.objinferior_expr == m.dual_obj_expr)
    
    def objinferior(model):
        return m.objinferior_expr #m.custo_agregador_linear + m.custo_outros 
    m.ObjInferior = Objective(rule=objinferior, sense=minimize, name='Custo_ISO') # eq. 1a (lower)
    
    m.Lagr = Expression(expr=build_lagrangian_expression(m, DIC_DUAIS))
    #m.Lagr = Expression(expr=build_lagrangian_expression_from_expr(m, DIC_DUAIS, objective_expr=m.PrimalObjExpr))
    add_lagrangian_derivatives_from_constraints(m, m.Lagr, LISTA_PRIMAIS) #deixar dL/dx sem complementares, dá unbounded.
    add_complementarity_slackness_condition_linear_optimized(m, DIC_DUAIS, M_value=big_m_value) #deixar só compl, sem dL/dx, dá unbounded.
    LISTA_CONSTRS = [c for c in m.component_objects(Constraint, active=True)]
    a = pd.DataFrame([v.name for v in LISTA_CONSTRS])
    m.ObjInferior.deactivate()

    m.BESSQuantMax = Constraint(doc='Max quantidade de BESS', name='BESSCapMax',
        rule=lambda model: sum(model.w[j, y, k] for j in model.j for y in model.y for k in model.K_inv) <= model.W_max)


    m.investimento = Expression(expr=sum(
        ( (m.tab_custo[k] * m.w[j, y, k]) / ((1 + m.tma) ** (y - 1)) for j in m.j for y in m.y for k in m.K_inv)),
        doc='Custo de investimento do agregador de BESSs no sistema de distribuição'
    )

    m.ObjSuperior = Objective(expr= m.custo_agregador_linear - m.investimento, sense=maximize, # eq. 1a (upper)
                            doc='Maximização do VPL do agregador')

    m.InvestOnce = Constraint(m.j, name='InvestOnce',
        doc='Ensures that investment happens at most once per location',
        rule=lambda m, j: sum(m.w[j, y, k] for y in m.y for k in m.K_inv) <= 1)

    return m
