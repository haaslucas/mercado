from amplpy import AMPL, Environment
import pathlib
# (1) Onde estão Novo.mod, input.dat, execute.run, *.inc…?
basedir = pathlib.Path(".")   # ajuste o caminho

# (2) Cria objeto AMPL; baixa o binário community se não existir
ampl = AMPL()             # o cwd do AMPL vira basedir
ampl.setOption('solver', 'ipopt')            # ou ipopt, cplex, gurobi…


# --- parte inicial do seu run ---------------------------------
ampl.eval("""
reset;
model Novo.mod;
data  input.dat;
""")

# muda parâmetros direto no Python
ampl.getParameter('U').set(1.0)
ampl.getParameter('Remove_Energy_Mkt').set(1)
ampl.getParameter('Remove_Carbon_Mkt').set(1)
ampl.getParameter('Remove_Ramp').set(0)

# executa o restante do script
ampl.eval('include adjust_data.run;')

# … aqui você pode injetar outros blocos AMPL ou
# chamar funções Python para ler arquivos, gerar gráficos etc.

# finalmente resolvemos
ampl.setOption('solver', 'conopt')
ampl.solve()

# coleta resultados
print("DSO Total Costs  :", ampl.getValue('sum {t in T, s in S} (probability[s,t] *Bid[t,s] * P_DSO[t,s]) \
                                         + sum {g in G_D, t in T, s in S}(probability[s,t] *GenCosts_dist[g,t,s]) \
                                         + sum {t in T, s in S} (probability[s,t] *Carbon_Price[t,s] * Carbon_SE[t,s])'))

# ou simplesmente pegar o parâmetro calculado no run file:
print("DSO_Costs table →")
print(ampl.getData('DSO_Costs').head())