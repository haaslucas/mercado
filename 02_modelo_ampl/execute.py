from amplpy import AMPL, Environment
import pathlib
import amplpy

# (1) Onde estão Novo.mod, input.dat, execute.run, *.inc…?
basedir = pathlib.Path(".")   # ajuste o caminho

# (2) Cria objeto AMPL; baixa o binário community se não existir
ampl = AMPL()             # o cwd do AMPL vira basedir
ampl.setOption('solver', 'ipopt')            # ou ipopt, cplex, gurobi…

# (3) EXECUTA exatamente o que você faria no prompt AMPL ↓↓↓
ampl.eval('include execute.run;')             # chama o .run

# -----------------------------------------------
# Se o .run já contém 'solve', acabou!
# Podemos coletar qualquer coisa que ele calculou:
print("Tempo total de solve (s):",
      ampl.getValue('_total_solve_elapsed_time'))

for entity in amplpy.EntityMap(ampl):
    # Exibe o nome e o valor de cada entidade
    print(f"{entity.name}: {entity.value}")

# --- Como pegar a lista de todos os parâmetros ---
param_names = []
try:
    # Para versões mais recentes do amplpy (>= 0.15), getParameters() retorna um iterador de objetos Parameter
    param_objects = ampl.getParameters()
    param_names = [p.name for p in param_objects]
except AttributeError:
    # Para versões mais antigas, getParameters() pode retornar um dicionário
    # onde as chaves são os nomes dos parâmetros.
    try:
        param_map = ampl.getParameters()
        param_names = list(param_map.keys())
    except Exception as e_dict:
        print(f"Erro ao tentar obter parâmetros como dicionário: {e_dict}")
except Exception as e_iter:
    print(f"Erro inesperado ao obter parâmetros: {e_iter}")


print("\nLista de todos os parâmetros do modelo:")
if param_names:
    for name in sorted(param_names): # Ordenado para melhor visualização
        try:
            # Para obter mais informações, como a dimensionalidade:
            p = ampl.getParameter(name)
            ndim = p.numIndices()
            print(f" - {name} (dimensões: {ndim})")
        except Exception as e_detail:
            # Se o parâmetro foi removido ou renomeado dinamicamente e o nome ainda estiver na lista inicial
            print(f" - {name} (não foi possível obter detalhes: {e_detail})")
else:
    print("Nenhum parâmetro encontrado ou erro ao listá-los.")

# Alternativamente, você pode usar o comando AMPL diretamente:
# print("\nLista de parâmetros via comando AMPL 'print _PARMS;':")
# try:
#     ampl.eval("print _PARMS;")
#     print(ampl.getOutput())
# except Exception as e_ampl_cmd:
#     print(f"Erro ao executar 'print _PARMS;': {e_ampl_cmd}")


# parâmetro/var/constraint → DataFrame
dso_costs = ampl.getData('DSO_Costs')
print(dso_costs.head())