import pyomo.environ as pyo
import pyomo.opt as opt
import pandas as pd



pasta = './01_modelo_pyomo/data/'
I_t = pd.read_csv(pasta + 'matriz_incidencia_transmissao.csv')
I_t = I_t.values.tolist()

dlines  = pd.read_csv(pasta + 'distribution_line_data.csv')

dlines
dlines = dlines.set_index(['From', 'To'])

tnodes = pd.read_csv(pasta + 'transmission_node_data.csv')
dnodes = pd.read_csv(pasta + 'distribution_node_data.csv')
ess = pd.read_csv(pasta + 'ess.csv',index_col=0)

dload = pd.read_csv(pasta + 'load_shape_data.csv')
dload_scn = pd.read_csv(pasta + 'scenario_distribution_connected_loads.csv')
dload2 = pd.read_csv(pasta + 'dist_load.csv')

a=1
SBASE = 100               # MVA
VBASE = 24.9              # kV
IBASE = 100 * 1e3 / 24.90 # A
SE_Capacity = 50          # MWh
U=1                       # Rigidez da regulação de emissões de carbono   




g_t_nd = pd.read_csv(pasta + 'transmission_connected_non_dispatchable_generators.csv')
g_t_d  = pd.read_csv(pasta + 'transmission_connected_dispatchable_generators.csv')
g_d_nd = pd.read_csv(pasta + 'dso_connected_non_dispatchable_generators.csv')
g_d_d  = pd.read_csv(pasta + 'dso_connected_dispatchable_generators.csv')
ls = pd.read_csv(pasta + 'ls.csv')



#nó de conexão entre a rede de distribuição e a rede de transmissão
N_INF_d = 1
N_INF_t = 5

um, LLmodel = pyo.ConcreteModel(), pyo.ConcreteModel()

#################### Índices ####################   
b = range(0,5)  # Índice para bloco de oferta  
g = range(0,8)  # Índice para geradores  

i, j = range(0,14), range(0,14)  # Índice para nós  
ij, ki = range(0,40), range(0,40)  # Índice para linhas de distribuição  
l = range(0,20)  # Índice para linhas de transmissão  
#s = range(0,5)  # Índice para cenários  

#################### Conjuntos ####################
#create a set from dlines double index
um.L = pyo.Set(initialize=dlines.index.tolist()) # Set of distribution lines
um.Lfrom = pyo.Set(initialize=dlines.index.get_level_values(0).tolist()) # Set of distribution lines
um.Lto = pyo.Set(initialize=dlines.index.get_level_values(1).tolist()) # Set of distribution lines
um.N = pyo.RangeSet(1,34) # Set of transmission nodes
um.B_g = pyo.Set()          # Set of bids provided by generator g
#um.G_D = pyo.Set(list(g_d_d['Node'].values))      # Set of distribution-connected generators
um.G_D = pyo.Set(initialize=g_d_d['Node'].tolist())      # Set of distribution-connected generators
um.G_T = pyo.RangeSet(1, 5)      # Set of transmission-connected generators
um.G_i = pyo.Set()          # Set of generators connected to node i
um.L_D = pyo.Set()          # Set of load shifting nodes
um.D = pyo.Set()            # Set of distribution lines
um.ess = pyo.RangeSet(1, 2)      # Set of energy storage systems
um.d_nodes = pyo.RangeSet(1, 34) # Set of distribution nodes


um.N_T = pyo.Set()          # Set of transmission nodes
um.N_INF = pyo.Set()          # Set of interface nodes
um.S = pyo.RangeSet(1, 5)      # Set of scenarios
um.T = pyo.RangeSet(1, 24)      # Set of time periods


#################### Parametros ####################

# Parâmetros relacionados a geradores
um.A_i_l = pyo.Param(i,l,initialize= lambda m,col,lin: I_t[lin][col])# Incidence matrix of transmission lines

# Parâmetros relacionados a lances de geradores

um.P_D_g_b = pyo.Param(um.G_D, b, initialize=lambda m, g, b: 1)           # Size of the generator's electricity block bid
um.P_T_g_b = pyo.Param(um.G_D, b, initialize=lambda m, g, b: 1)           # Size of the generator's electricity block bid

# Parâmetros relacionados a geradores
um.S_D_g = pyo.Param(um.G_D, initialize=lambda m, g: 1)                   # Maximum power injection of generators
um.S_T_g = pyo.Param(um.G_D, initialize=lambda m, g: 1)                   # Maximum power injection of generators
um.gamma_D_g = pyo.Param(um.G_D, initialize=lambda m, g: 1)               # Generator's carbon intensity
um.gamma_T_g = pyo.Param(um.G_D, initialize=lambda m, g: 1)               # Generator's carbon intensity
um.nu_D_g_t_s = pyo.Param(um.G_D, um.T, um.S, initialize=lambda m, g, t, s: 1)  # Generator's carbon emission
um.nu_T_g_t_s = pyo.Param(um.G_D, um.T, um.S, initialize=lambda m, g, t, s: 1)  # Generator's carbon emission

# Parâmetros relacionados a demanda e disponibilidade
um.L_D_i_t_s = pyo.Param(i, um.T, um.S, initialize=lambda m, i, t, s: 1)   # Power demand of a distribution node
um.L_T_i_t_s = pyo.Param(i, um.T, um.S, initialize=lambda m, i, t, s: 1)   # Power demand of a transmission node
um.r_g_t_s = pyo.Param(g, um.T, um.S, initialize=lambda m, g, t, s: 1)     # Renewable-based power availability

# Parâmetros relacionados a linhas e capacidades
um.I_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Capacity of distribution lines
um.PF_l = pyo.Param(l, initialize=lambda m, l: 1)                    # Power flow capacity of transmission lines
um.R_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Resistance of distribution lines
um.X_l = pyo.Param(l, initialize=lambda m, l: 1)                     # Reactance of transmission lines
um.X_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Reactance of distribution lines


# Parâmetros de mercado e regulação
um.p_s_t = 1            # Weight of scenario s at period t
um.pi_E_i_t_s = pyo.Param(i, um.T, um.S, initialize=lambda m, i, t, s: 1)  # Bids on the electricity market
um.pi_C_t_s = pyo.Param(um.T, um.S, initialize=lambda m, t, s: 1)          # Bids on the carbon market
um.U = pyo.Param(initialize=lambda m: 1)                             # Rigidity of the carbon emission regulation
um.Lambda = pyo.Param(initialize=lambda m: 1)                        # Average carbon intensity of the entire power system

# Parâmetros das linhas de distribuição dlines = distribution_line_data.csv
um.R = pyo.Param(i, initialize=lambda m,i: dlines['R (pu)'])                 # branch resistance (ohm)
um.X = pyo.Param(i, initialize=lambda m,i: dlines['X (pu)'])                 # branch reactance (ohm)
#um.Z = pyo.Param(i, initialize=lambda m,i: dlines['X (pu)'])                 # branch reactance (ohm)
um.Imax = pyo.Param(i, initialize=lambda m,i: dlines['Line Ampacity (A)'])   # branch maximum current (kA)

'''# Parâmetros relacionados a BESS do sistema de distribuição ESS
um.ESS_P_max_charge = pyo.Param(i, initialize=lambda m,i: dr['Max ESS Charge (MW)']) # Max charging power of BESS
um.ESS_P_max_discha = pyo.Param(i, initialize=lambda m,i: dr['Max ESS Discharge (MW)']) # Max discharging power of BESS
um.ESS_SOC_min = pyo.Param(i, initialize=lambda m,i: dr['Min ESS SOC (MWh)']) # Min SOC of BESS
um.ESS_SOC_max = pyo.Param(i, initialize=lambda m,i: dr['Max ESS SOC (MWh)']) # Max SOC of BESS
um.ESS_SOC_0   = pyo.Param(i, initialize=lambda m,i: dr['Ini ESS SOC (MWh)']) # Initial SOC of BESS'''

um.P_LS_i = pyo.Param(i, initialize=lambda m, i: 1)                  # Maximum nodal load increase/reduction


#################### Variáveis de decisão ####################

# Variáveis relacionadas à distribuição

um.I = pyo.Var(um.L,  um.T, um.S) # Distibution current flow of lines
um.C_D_g_b = pyo.Var(um.T, um.S)           # Distribution-connected generators' bids
um.C_T_g_b = pyo.Var(um.T, um.S)           # Transmission-connected generators' bids

um.rho_D_g_b_t_s = pyo.Var(um.T, um.S) # Accepted bid of block b
um.P_D_g_t_s = pyo.Var(um.G_D, um.T, um.S) # Active power dispatch of distribution-connected generators
um.rho2_SE_t_s = pyo.Var(um.T, um.S) # Carbon allowances absorbed/provided by generator g
um.lambda_t_s = pyo.Var(um.T, um.S) # LMP of energy / Carbon allowance trades in the wholesale market
um.P_SE_i_t_s = pyo.Var(um.T, um.S)   # Active power trades in the wholesale market
um.upsilon_t_s = pyo.Var(um.T, um.S) # LMP of carbon / Carbon allowance trades in the wholesale market

um.Q_D_g_t_s = pyo.Var()   # Generator reactive power dispatch

um.P_ki_t_s = pyo.Var(um.L, um.T, um.S) # Distibution active power flow of comming lines
um.P_ij_t_s = pyo.Var(um.L, um.T, um.S) # Distibution active power flow of going lines
um.Q_ij_t_s = pyo.Var()    # Reactive power flow on distribution lines
um.P_ESS_i_t_s = pyo.Var(um.ess, um.T, um.S) # Power charge/discharge of storage systems
um.P_LS_i_t_s = pyo.Var(um.N, um.T, um.S)  # Nodal load shift
um.V_SQ_i_t_s = pyo.Var()  # Nodal voltage magnitude squared
um.SOC_i_t_s = pyo.Var()   # State of charge of energy storage systems

# Variáveis relacionadas à transmissão
um.P_T_g_t_s = pyo.Var()    # Generator's total power dispatch
um.PF_l_t_s = pyo.Var()     # Power flow on transmission lines
um.Q_SE_i_t_s = pyo.Var()   # Reactive power trades in the wholesale market



um.delta_i_t_s = pyo.Var()  # Voltage angle of transmission nodes
um.nu_T_g_t_s = pyo.Var()   # Generator's carbon emission
um.rho_T_g_b_t_s = pyo.Var() # Accepted bid of block b
um.upsilon_T_g_t_s = pyo.Var() # Carbon allowances absorbed/provided by generator g


# Função Objetivo para o DSO (Upper-level problem)
um.objective = pyo.Objective(
    expr=sum(
        um.p_s_t * (
            # Custos de geração dos geradores conectados à rede de distribuição
            um.C_D_g_b[t, s] * um.rho_D_g_b_t_s[t, s] +
            # Receita de troca de allowances de carbono com o ISO
            um.upsilon_t_s[t, s] * um.rho2_SE_t_s[t, s] +
            # Receita de intercâmbio de energia com o ISO nos nós de interface
            sum( um.P_SE_i_t_s[i, t, s] * um.lambda_t_s[t, s] for i in um.N_INF )
            ) for s in um.S  for t in um.T),
    sense = pyo.minimize
)

"""
Eq. (2): Garante o balanço de potência ativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência, armazenamento e intercâmbio com o sistema de transmissão.
"""

#a = dlines[dlines.index[x][0] == 4 for x in dlines.index]

def active_power_balance_distribution_rule(um, n, t, s, e):


    pot_Geradores   = sum(um.P_D_g_t_s[g, t, s] for g in g_d_d['Node'] if g == n)

    pot_Entrando    = sum(um.P_ij_t_s[n, t, s] if n == dlines.xs(l).index else 0 for l in um.L)

    pot_Saindo      = sum(um.P_ij_t_s[l, t, s]   + dlines['R (pu)'].loc[l] * um.I[l,t,s] for l in um.L )

    pot_ESS         = sum(um.P_ESS_i_t_s[e , t, s] for bess in ess['Node'] if bess == n)

    pot_Intercambio = sum(um.P_SE_i_t_s[i, t, s] for i in um.N_INF if i == n)

    pot_LoadShift   = sum(um.P_LS_i_t_s[i, t, s] for i in ls['Node'] if i == n)

    pot_LoadDemand  = dload2['Load'][(dload2['Node'] == n) & (dload2['Hour'] == t) & (dload2['Scenario'] == s)].iloc[0]

    return pot_Geradores + pot_Entrando - pot_Saindo + pot_ESS + pot_Intercambio == pot_LoadDemand + pot_LoadShift 

um.active_power_balance_distribution = pyo.Constraint(um.N, um.T, um.S, um.ess,  rule=active_power_balance_distribution_rule)




solver = opt.SolverFactory('ipopt')
results = solver.solve(um, tee=True)

# Análise de resultados
for g in um.G_D:
    print(f'A potência gerada pelo gerador {g} é: {um.P_D_g_t_s[g, :, :].value}')

for t in um.T:
    print(f'O preço de mercado na hora {t} é: {um.lambda_t_s[t, :].value}')

#print result of um.objective
print(f'O valor da função objetivo é: {um.objective()}')
a=1