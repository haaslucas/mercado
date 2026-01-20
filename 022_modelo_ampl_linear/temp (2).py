from pyomo.environ import *
import pandas as pd
import sys
import time
import cProfile
import pstats
from datetime import datetime
from pathlib import Path
import inspect
import shutil

sys.path.append('.')
from modelos import agregator_inv_kkt7 as agregator #agregator_inv_dual  agregator_inv_dualV4 agregator_inv_dualV4cones
from src.auxiliares_kkt import *
from src.auxiliares_main import save_warmstart, load_warmstart, configure_gurobi_solver
from src.auxiliares_dual import export_pyomo
from src.auxiliares_results import analyze_complementarity

nome = "ag" #4 é con custo 2000. 5 é com custo 20000

# --- Diretório de Resultados ---
run_dir = Path("./resultados") / f"{nome}_{datetime.datetime.now().strftime('%Y%m%d_%H%M')}"
run_dir.mkdir(parents=True, exist_ok=True)
print(f"--- Resultados serão salvos em: {run_dir} ---")

# Salva uma cópia do script do modelo no início da execução
try:
    source_file = inspect.getsourcefile(agregator)
    if source_file and Path(source_file).exists():
        shutil.copy(source_file, run_dir)
        print(f"--- Cópia do script do modelo '{Path(source_file).name}' salva em '{run_dir}' ---")
except Exception as e:
    print(f"--- AVISO: Não foi possível copiar o arquivo fonte do modelo. Erro: {e} ---")

solve_and_save_relaxed_model = False
use_warmstart = False
warm_start_file = run_dir / f'{nome}_warmstart.pkl'
# --- Início do Profiling e Timers ---
profiler = cProfile.Profile()
profiler.enable()

print("--- Iniciando a construção do modelo ----")
BEGIN_BUILD = time.time()

mpec = agregator.agregator_bilevel()
prelog=False 
poslog=False
prelog_simples=False

solver, log_file_path = configure_gurobi_solver(run_dir)

END_BUILD = time.time()
BUILD_TIME = END_BUILD - BEGIN_BUILD
print(f'--- Tempo para construir o modelo: {BUILD_TIME:.2f} segundos ---')

# --- Opcional: Resolver e salvar a relaxação do modelo ---
if solve_and_save_relaxed_model:
    print("\n" + "="*80)
    print(" Resolvendo e salvando a relaxação do modelo (sem variáveis inteiras)")
    print("="*80)
    
    mpec_relaxed = mpec.clone()
    TransformationFactory('core.relax_integrality').apply_to(mpec_relaxed)

    # Configura um solver para o modelo relaxado, com um log file separado
    relaxed_log_file = run_dir / 'gurobi_relaxed.log'
    solver_relaxed, _ = configure_gurobi_solver(run_dir, options_dict={'LogFile': str(relaxed_log_file)})

    print("\n--- Iniciando a resolução do modelo relaxado ---")
    results_relaxed = solver_relaxed.solve(mpec_relaxed, tee=True)
    
    if results_relaxed and (results_relaxed.solver.status == SolverStatus.ok or results_relaxed.solver.status == SolverStatus.warning):
        print("--- Exportando resultados do modelo relaxado ---")
        relaxed_results_dir_name = f"{run_dir.name}_relaxed"
        export_pyomo(mpec_relaxed, nome_arquivo_resultados=relaxed_results_dir_name, 
                     to_pkl=True, to_json=False, to_xlsx=False)
        print(f"--- Resultados do modelo relaxado salvos em 'resultados/{relaxed_results_dir_name}' ---")
    else:
        print("--- AVISO: A resolução do modelo relaxado falhou ou não encontrou solução viável. ---")
    print("="*80)

# --- Warm Start from Pickle File ---
if use_warmstart:
    load_warmstart(mpec, warm_start_file)

if prelog_simples:
    log_model_equations_simple(mpec, filename='log_equacoes_simples.md', output_dir=run_dir)

if prelog:
    logging.basicConfig(level=logging.INFO)
    log_txt = logging.getLogger('infeas') # create a logger for diagnostics
    log_txt.setLevel(logging.INFO)             # ensure INFO messages are emitted
    log_txt.propagate = False  # prevent messages from propagating to the root logger
    file_hdl = logging.FileHandler(run_dir / 'logPre.txt', mode='w') # create a handler for the log file
    file_hdl.setFormatter(logging.Formatter('%(message)s'))  # just the message
    log_txt.addHandler(file_hdl) # add the handler to the logger
    log_infeasible_constraints(mpec, logger=log_txt, log_variables=True, log_expression =True) # trigger diagnostics
    log_txt.removeHandler(file_hdl) #(optional) remove handler or close if the script continues running
    file_hdl.close() # close the handler

print("\n--- Iniciando a resolução do modelo ---")
BEGIN_SOLVE = time.time()
results = None
try:
    results = solver.solve(mpec, tee=True, warmstart=False) # ,options=gurobi_options) , logfile='test.log'

    if use_warmstart:
        save_warmstart(mpec, warm_start_file)
except (Exception, KeyboardInterrupt) as e:
    print(f"\nWARNING: Solve interrupted or failed: {e}.")
    if use_warmstart:
        print("Attempting to save warmstart from available solution.")
        save_warmstart(mpec, warm_start_file)
    pass

#mpecfixed = mpec.clone() # create a copy of the original model

#mpecfixed = fix_binaries(mpecfixed) # fix binary variables to prevent gurobi from changing them
#results_fixed = solver.solve(mpecfixed, tee=True ) 


END_SOLVE = time.time()
SOLVE_TIME = END_SOLVE - BEGIN_SOLVE
print(f'--- Tempo para resolver o modelo: {SOLVE_TIME:.2f} segundos ---')

if poslog:
    logging.basicConfig(level=logging.INFO)
    log_txt = logging.getLogger('infeas')
    log_txt.setLevel(logging.INFO)             # ensure INFO messages are emitted
    file_hdl = logging.FileHandler(run_dir / 'logPos.txt', mode='w') # create a handler for the log file
    file_hdl.setFormatter(logging.Formatter('%(message)s'))  # just the message
    log_txt.addHandler(file_hdl)
    log_infeasible_constraints(mpec, logger=log_txt, log_variables=True, log_expression =True)
    log_txt.removeHandler(file_hdl)
    file_hdl.close()
if poslog:
    log_model_equations_complete(mpec, filename='log_equacoes_completas.md', output_dir=run_dir)

if results and (results.solver.status == SolverStatus.ok or results.solver.status == SolverStatus.warning):

    print(f"\nSolver finalizou com status '{results.solver.status}' e condição de término '{results.solver.termination_condition}'.")

    # --- Análise das Condições de Complementaridade ---
    print("\n--- Analisando condições de complementaridade (Big-M) ---")
    try:
        comp_analysis_df = analyze_complementarity(mpec)
        if not comp_analysis_df.empty:
            issues_df = comp_analysis_df[comp_analysis_df['issue'] != '']
            if not issues_df.empty:
                #print("AVISO: Foram encontrados problemas potenciais na formulação Big-M:")
                #with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000):
                #    print(issues_df)
                
                # Salvar os problemas encontrados em um arquivo de log
                log_path = run_dir / 'big_m_analysis_issues.log'
                with open(log_path, 'w') as f:
                    f.write("Problemas potenciais encontrados na formulação Big-M:\n\n")
                    f.write(issues_df.to_string())
                print(f"Log de problemas da análise Big-M salvo em: {log_path}")
            else:
                print("Nenhum problema aparente encontrado na formulação Big-M.")
            
            # Salvar a análise completa em um arquivo CSV
            analysis_csv_path = run_dir / 'complementarity_analysis.csv'
            comp_analysis_df.to_csv(analysis_csv_path, index=False)
            print(f"Relatório de análise de complementaridade salvo em: {analysis_csv_path}")
        else:
            print("Nenhuma restrição de complementaridade para analisar foi encontrada.")
    except Exception as e:
        print(f"--- AVISO: Falha ao executar a análise de complementaridade. Erro: {e} ---")


    if results.solver.termination_condition == TerminationCondition.optimal:
        print("Solução ótima encontrada.")
        to_pkl, to_xlsx = True, False
    else:
        # Isso cobre maxTimeLimit, userInterrupt, etc.
        print("Solução ótima não encontrada. Caso tenha sido encontrada uma solução viável, ela será salva.")
        to_pkl, to_xlsx = True, False

    metadata = {
        'model_name': mpec.name,
        'model_source_file': inspect.getsourcefile(agregator),
        'run_datetime': datetime.datetime.now().isoformat(),
        'build_time_seconds': BUILD_TIME,
        'solve_time_seconds': SOLVE_TIME,
        'solver_status': str(results.solver.status),
        'termination_condition': str(results.solver.termination_condition),
        'solver_options': {k: str(v) for k, v in solver.options.items()},
        'gurobi_log_file': log_file_path.name,
    }
    export_pyomo(mpec, nome_arquivo_resultados=run_dir.name, metadata=metadata,
                 to_parquet=False, to_xlsx=to_xlsx, to_pkl=to_pkl, to_json=False, to_html=False)

    # --- Salvar solução para warm start ---
    if use_warmstart:
        save_warmstart(mpec, warm_start_file)

    #mpec.ObjSuperior.display()

else:
    try:
        # Tenta salvar o que for possível mesmo sem solução
        metadata = {
            'model_name': mpec.name,
            'model_source_file': inspect.getsourcefile(agregator),
            'run_datetime': datetime.datetime.now().isoformat(),
            'build_time_seconds': BUILD_TIME,
            'solve_time_seconds': SOLVE_TIME,
            'solver_status': str(results.solver.status) if results else 'failed',
            'termination_condition': str(results.solver.termination_condition) if results else 'failed',
            'solver_options': {k: str(v) for k, v in solver.options.items()},
            'gurobi_log_file': log_file_path.name if 'log_file_path' in locals() else 'not_created',
        }
        export_pyomo(mpec, nome_arquivo_resultados=run_dir.name, metadata=metadata,
                     to_parquet=False, to_xlsx=False, to_pkl=True, to_json=False, to_html=False)
    except Exception as e:
        print(f"Erro ao exportar resultados: {e}")
    with open(run_dir / 'log_infeasible_constraints.txt', 'w') as f:
        print("The solver did not find an optimal solution.")
        log_infeasible_constraints(mpec, log_variables=True)
    if results:
        print(f"Status do solver: {results.solver.status}, Condição de término: {results.solver.termination_condition}")


# --- Fim do Profiling e Análise ---
profiler.disable()
#print("\n--- Análise de Performance da Execução Completa (Construção + Resolução) ---")
stats = pstats.Stats(profiler).sort_stats('cumtime')
#stats.print_stats(20)  # Imprime as 20 funções mais lentas

# Salvar resultados do profiler em .prof e .csv
profiler_prof_path = run_dir / 'model_build.prof'
profiler_csv_path = run_dir / 'model_build_stats.csv'
stats.dump_stats(profiler_prof_path)

# Converter para CSV (na verdade, é um arquivo de texto formatado)
with open(profiler_csv_path, 'w', newline='') as f:
    original_stream = stats.stream
    stats.stream = f
    stats.print_stats()
    stats.stream = original_stream
#print("Iniciando a gerar o html...")
export_pyomo(mpec, nome_arquivo_resultados=run_dir.name, launch_web_viewer=False, to_json=False)

# Salva uma cópia na pasta geral de resultados do modelo, para servir como "latest"
print(f"--- Salvando cópia do webviewer em ./resultados/{nome} ---")
export_pyomo(mpec, nome_arquivo_resultados=nome, launch_web_viewer=False, to_json=False)

print(f"Resultados do profiler salvos em '{profiler_prof_path}' e '{profiler_csv_path}'")
print(f"Para uma análise interativa, use: python -m snakeviz {profiler_prof_path}\n")
a=1

