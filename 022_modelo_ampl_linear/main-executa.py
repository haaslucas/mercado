'''
Este main resolve o problema do slide da Lecture 7, slide 74, que é um problema de 
oferta estratégica de um grupo gerador em um sistema de energia elétrica.
É um problema binível, onde o problema de nível inferior é o fechamento do mercado
e o problema de nível superior é a oferta estratégica do grupo gerador.
'''

from pyomo.environ import *
#from pyomo.util.infeasible import log_infeasible_constraints
import pandas as pd
import sys
sys.path.append('.')
import ieee34busV10 as ieee34bus

from auxiliares2 import *
import time

BEGIN = time.time()
mpec = ieee34bus.ieee34bus(prelog=True, poslog=True, analyze_model=False, Y = 40) #
#mpec.obj.deactivate()
prelog=True 
poslog=True

solver = SolverFactory('gurobi')
#solver.options['OutputFlag'] = 1 
#solver.options['NonConvex'] = 2
#solver.options['iisfind'] = 1 
#solver.options['Heuristics'] = 0.9

#solver.options['BarHomogeneous'] = 1    
    
#solver.options['Preqlinearize'] = 1
#ver opção pra conica

# Dicionário com opções específicas para o Gurobi
#gurobi_options = {'TimeLimit': 120 } # Limite de 60 segundos

# --- CONFIGURANDO AS OPÇÕES ---


# A OPÇÃO MAIS IMPORTANTE: ativar o detetive de inviabilidade
#solver.options['iisfind'] = 1 
# Usar estratégias padrão para as outras opções
#solver.options['mipfocus'] = 0
#solver.options.nu  mericfocus = 0 # Outra sintaxe válida, mas a de dicionário é mais comum

# Passe o dicionário para a função solve
#solver = SolverFactory('cplex', executable='C:\\AMPL\\cplex.exe') # ipot, conopt, 


END = time.time()
TOTAL_TIME = END - BEGIN
print(f'Tempo para montar o modelo: {TOTAL_TIME:.2f} segundos')    

if prelog:
    logging.basicConfig(level=logging.INFO)
    log_txt = logging.getLogger('infeas') # cria um logger para o diagnóstico
    log_txt.setLevel(logging.INFO)             # garante que msgs INFO sejam emitidas
    log_txt.propagate = False  # evita que as mensagens sejam propagadas para o logger raiz
    file_hdl = logging.FileHandler('./resultados/log_pre_solve.txt', mode='w') # cria um handler para o arquivo de log
    file_hdl.setFormatter(logging.Formatter('%(message)s'))  # só a mensagem
    log_txt.addHandler(file_hdl) # adiciona o handler ao logger
    log_infeasible_constraints(mpec, logger=log_txt, log_variables=True, log_expression =True) # dispara o diagnóstico
    log_txt.removeHandler(file_hdl) #(opcional) remove handler ou fecha se o script continuar rodando
    file_hdl.close() # fecha o handler

results = solver.solve(mpec, tee=True ) # ,options=gurobi_options) , logfile='teste.log'

mpec = fix_binaries(mpec) # fixa as variaveis binarias para evitar que o gurobi altere elas
results_fixed = solver.solve(mpec, tee=True ) 


END = time.time()
TOTAL_TIME = END - BEGIN
print(f'Tempo para resolver: {TOTAL_TIME:.2f} segundos')    

if poslog:
    logging.basicConfig(level=logging.INFO)
    log_txt = logging.getLogger('infeas')
    log_txt.setLevel(logging.INFO)             # garante que msgs INFO sejam emitidas
    file_hdl = logging.FileHandler('./resultados/log_pos_solve.txt', mode='w') # cria um handler para o arquivo de log
    file_hdl.setFormatter(logging.Formatter('%(message)s'))  # só a mensagem
    log_txt.addHandler(file_hdl)
    log_infeasible_constraints(mpec, logger=log_txt, log_variables=True, log_expression =True)
    log_txt.removeHandler(file_hdl)
    file_hdl.close()
    #log_model_infeasibility(m, './resultados/log_pos_solve.txt')
    
# 6. Exportar e analisar os resultados
if results.solver.status == SolverStatus.ok and results.solver.termination_condition == TerminationCondition.optimal:
    print("Solução ótima encontrada.")

    # --- Comparação de Duais ---
    # Extrai os valores de lambda_BPA (variáveis primais do MPEC)
    lambda_bpa_vals = []
    for idx in mpec.lambda_BPA:
        lambda_bpa_vals.append({'Key': idx, 'lambda_BPA_MPEC': value(mpec.lambda_BPA[idx])})
    df_lambda_bpa = pd.DataFrame(lambda_bpa_vals).set_index('Key')

    # Extrai os duais verdadeiros (LMPs) do solver após fixar os binários
    lmp_vals = []
    if hasattr(mpec, 'dual'):
        for idx in mpec.BALANCO_POT_ATIVA:
            # Check if dual exists for the constraint before accessing
            if mpec.BALANCO_POT_ATIVA[idx] in mpec.dual:
                lmp_vals.append({'Key': idx, 'LMP_Solver': mpec.dual[mpec.BALANCO_POT_ATIVA[idx]]})
    df_lmp = pd.DataFrame(lmp_vals).set_index('Key')

    # Juntar para comparação
    df_compare = df_lambda_bpa.join(df_lmp)

    export_pyomo(
        mpec,
        nome_arquivo_resultados="ieee34bus",
        to_pkl=True,
        to_xlsx=True,
        to_parquet=True,
        extra_dfs={'LMP_Comparison': df_compare.reset_index()}
    )
    #print mpec.display() in a txt
    '''with open('./resultados/ieee34bus_display.txt', 'w') as f:
        mpec.display(ostream=f)

    with open('./resultados/ieee34bus_pprint.txt', 'w') as f:
        mpec.pprint(ostream=f)
    #modelo_to_latex(mpec, symbol_map = mpec.SYMBOL_MAP, filename="linearizacao.tex")'''
    print("Resultados do modelo:")
    mpec.ObjInferior.display()


else:
    with open('./resultados/log_infeasible_constraints.txt', 'w') as f:
        print("O solver não encontrou uma solução ótima.")
        log_infeasible_constraints(mpec, log_variables=True)


a=1
