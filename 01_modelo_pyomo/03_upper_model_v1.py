'''
1a tradução de todas equações do upper model
'''


import pyomo.environ as pyo
import pyomo.opt as opt
import pandas as pd

pasta = './01_modelo_pyomo/data/'
I_t = pd.read_csv(pasta + 'matriz_incidencia_transmissao.csv')
I_t = I_t.values.tolist()

dlines  = pd.read_csv(pasta + 'distribution_line_data.csv')
tnodes = pd.read_csv(pasta + 'transmission_node_data.csv')
dnodes = pd.read_csv(pasta + 'distribution_node_data.csv')
dload = pd.read_csv(pasta + 'load_shape_data.csv')
g_t_nd = pd.read_csv(pasta + 'transmission_connected_non_dispatchable_generators.csv')
g_t_d  = pd.read_csv(pasta + 'transmission_connected_dispatchable_generators.csv')
g_d_nd = pd.read_csv(pasta + 'dso_connected_non_dispatchable_generators.csv')
g_d_d  = pd.read_csv(pasta + 'dso_connected_dispatchable_generators.csv',index_col=0)

SBASE = 100               # MVA
VBASE = 24.9              # kV
IBASE = 100 * 1e3 / 24.90 # A
SE_Capacity = 50          # MWh
U=1                       # Rigidez da regulação de emissões de carbono   
#nó de conexão entre a rede de distribuição e a rede de transmissão
N_INF_d = 1
N_INF_t = 5

model, LLmodel = pyo.ConcreteModel(), pyo.ConcreteModel()

#################### Índices ####################   
b = range(0,5)  # Índice para bloco de oferta  
g = range(0,8)  # Índice para geradores  
i, j = range(0,14), range(0,14)  # Índice para nós  
ij, ki = range(0,40), range(0,40)  # Índice para linhas de distribuição  
l = range(0,20)  # Índice para linhas de transmissão  
s = range(0,5)  # Índice para cenários  
t = range(0,24) # Índice para períodos de tempo  

#################### Conjuntos ####################
model.B_g = pyo.Set()          # Set of bids provided by generator g
model.G_D = pyo.RangeSet(1, 5)      # Set of distribution-connected generators
model.G_i = pyo.Set()          # Set of generators connected to node i
model.G_T = pyo.Set()          # Set of transmission-connected generators
model.L_D = pyo.Set()          # Set of load shifting nodes
model.D = pyo.Set()            # Set of distribution lines
model.T = pyo.RangeSet(1,24)            # Set of transmission lines
model.N_D = pyo.Set()          # Set of distribution nodes
model.N_T = pyo.Set()          # Set of transmission nodes
model.N_INF = pyo.Set()        # Set of interface nodes
model.S = pyo.RangeSet(1,5)   
# Set of scenarios
model.T = pyo.Set()            # Set of time periods

model.N = pyo.RangeSet(1, 34) # Set of distribution nodes
model.L = pyo.Set(initialize=[(dlines['From'][i], dlines['To'][i]) for i in range(len(dlines))]) # Set of distribution lines



#################### Parametros ####################

# Parâmetros relacionados a geradores
model.A_i_l = pyo.Param(i,l,initialize= lambda model,col,lin: I_t[lin][col])# Incidence matrix of transmission lines

# model.C_D_g_b = pyo.Param(g, b, initialize=lambda model, g, b: 1)           # Distribution-connected generators' bids
model.C_D_g_b = pyo.Var(model.T, model.S) # Distribution-connected generators' bids
#model.C_T_g_b = pyo.Param(g, b, initialize=lambda model, g, b: 1)           # Transmission-connected generators' bids
model.P_D_g_b = pyo.Param(g, b, initialize=lambda model, g, b: 1)           # Size of the generator's electricity block bid
#model.P_T_g_b = pyo.Param(g, b, initialize=lambda model, g, b: 1)           # Size of the generator's electricity block bid

model.S_D_g = pyo.Param(g, initialize=lambda model, g: 1)                   # Maximum power injection of generators
#model.S_T_g = pyo.Param(g, initialize=lambda model, g: 1)                   # Maximum power injection of generators
model.gamma_D_g = pyo.Param(g, initialize=lambda model, g: 1)               # Generator's carbon intensity
#model.gamma_T_g = pyo.Param(g, initialize=lambda model, g: 1)               # Generator's carbon intensity
model.nu_D_g_t_s = pyo.Param(g, t, s, initialize=lambda model, g, t, s: 1)  # Generator's carbon emission
#model.nu_T_g_t_s = pyo.Param(g, t, s, initialize=lambda model, g, t, s: 1)  # Generator's carbon emission

# Parâmetros relacionados a demanda e disponibilidade
model.L_D_i_t_s = pyo.Param(i, t, s, initialize=lambda model, i, t, s: 1)   # Power demand of a distribution node
#model.L_T_i_t_s = pyo.Param(i, t, s, initialize=lambda model, i, t, s: 1)   # Power demand of a transmission node
model.r_g_t_s = pyo.Param(g, t, s, initialize=lambda model, g, t, s: 1)     # Renewable-based power availability

# Parâmetros relacionados a linhas e capacidades
model.I_ij = pyo.Param(i, initialize=lambda model,i: dlines['Line Ampacity (A)']) # Capacity of distribution lines
model.R_ij = pyo.Param(i, initialize=lambda model,i: dlines['R (pu)'])            # Resistance of distribution lines
model.X_ij = pyo.Param(i, initialize=lambda model,i: dlines['X (pu)'])            # Reactance of distribution lines
#model.PF_l = pyo.Param(l, initialize=lambda model, l: 1)                         # Power flow capacity of transmission lines
#model.X_l = pyo.Param(l, initialize=lambda model, l: 1)                          # Reactance of transmission lines

# Parâmetros relacionados a armazenamento e flexibilidade
model.P_ESS_i = pyo.Param(i, initialize=lambda model, i: 1)                 # Maximum charging and discharging power of storage systems
model.SOC_ESS_i = pyo.Param(i, initialize=lambda model, i: 1)               # Maximum and minimum state of charge of the storage systems
model.SOC_0_i = pyo.Param(i, initialize=lambda model, i: 1)                 # Initial state of charge of the storage systems
model.P_LS_i = pyo.Param(i, initialize=lambda model, i: 1)                  # Maximum nodal load increase/reduction

# Parâmetros de mercado e regulação
model.p_s_t = pyo.Param(s, t, initialize=lambda model, s, t: 1)             # Weight of scenario s at period t
model.pi_E_i_t_s = pyo.Param(i, t, s, initialize=lambda model, i, t, s: 1)  # Bids on the electricity market
model.pi_C_t_s = pyo.Param(t, s, initialize=lambda model, t, s: 1)          # Bids on the carbon market
model.U = pyo.Param(initialize=lambda model: 1)                             # Rigidity of the carbon emission regulation
model.Lambda = pyo.Param(initialize=lambda model: 1)                        # Average carbon intensity of the entire power system



#################### Variáveis de decisão ####################

# Variáveis relacionadas à distribuição
model.P_D_g_t_s = pyo.Var(model.G_D, model.T, model.S)   # Generator active power dispatch
model.Q_D_g_t_s = pyo.Var(model.G_D, model.T, model.S)   # Generator reactive power dispatch
model.P_ki_t_s = pyo.Var(model.d_lines, model.d_lines, model.T, model.S) # Distribution active power flow of coming lines
model.P_ij_t_s = pyo.Var(model.d_lines, model.d_lines, model.T, model.S) # Distribution active power flow of going lines
model.Q_ij_t_s = pyo.Var()    # Reactive power flow on distribution lines
model.P_ESS_i_t_s = pyo.Var(model.N_D, model.T, model.S) # Active power of energy storage systems
model.P_LS_i_t_s = pyo.Var(model.N_D, model.T, model.S) # Active power of load shifting
model.V_SQ_i_t_s = pyo.Var()  # Nodal voltage magnitude squared
model.SOC_i_t_s = pyo.Var()   # State of charge of energy storage systems
model.rho_D_g_b_t_s = pyo.Var(model.G_D, model.B_g, model.T, model.S) # Accepted bid of block b
model.upsilon_D_g_t_s = pyo.Var() # Carbon allowances absorbed/provided by generator g

# Variáveis relacionadas à transmissão
model.P_T_g_t_s = pyo.Var()    # Generator's total power dispatch
model.PF_l_t_s = pyo.Var()     # Power flow on transmission lines
model.P_SE_i_t_s = pyo.Var()   # Active power trades in the wholesale market
model.Q_SE_i_t_s = pyo.Var()   # Reactive power trades in the wholesale market
model.psi_t_s = pyo.Var(model.T, model.S) # Carbon allowances absorbed/provided by the transmission system
model.rho2_SE_t_s = pyo.Var(model.T, model.S) # ) # Carbon allowance trades in the wholesale market
model.lambda_i_t_s = pyo.Var() # LMP of electricity
model.psi_i_t_s = pyo.Var() # LMP of carbon
model.delta_i_t_s = pyo.Var()  # Voltage angle of transmission nodes
model.nu_T_g_t_s = pyo.Var()   # Generator's carbon emission
model.rho_T_g_b_t_s = pyo.Var() # Accepted bid of block b
model.upsilon_T_g_t_s = pyo.Var() # Carbon allowances absorbed/provided by generator g


# Função Objetivo para o DSO (Upper-level problem)
model.objective = pyo.Objective(
    expr=sum(
        model.p_s_t * (
            # Custos de geração dos geradores conectados à rede de distribuição
            sum(model.C_D_g_b[g, b] * model.rho_D_g_b_t_s[g, b, t, s] for g in model.G_D for b in model.B_g)
            # Receita de troca de allowances de carbono com o ISO
            + model.psi_t_s[t, s] * model.rho2_SE_t_s[t, s]
            # Receita de intercâmbio de energia com o ISO nos nós de interface
            + sum(model.lambda_i_t_s[i, t, s] * model.P_SE_i_t_s[i, t, s] for i in model.N_INF))
        for s in model.S 
        for t in model.T),
    sense=pyo.minimize)

####################################################
#################### RESTRIÇÕES ####################
####################################################


#################### Distribuição ##################

"""
Eq. (2): Garante o balanço de potência ativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência, armazenamento e intercâmbio com o sistema de transmissão.
"""
def active_power_balance_distribution_rule(model, i, t, s):

    pot_geradores = sum(model.P_D_g_t_s[g, t, s] for g in model.G_D if g in model.N_D[i])
    pot_entrando  = sum(model.P_ki_t_s[ki, t, s] for ki in model.D if ki in model.N_D[i])
    pot_saindo    = sum(model.P_ij_t_s[ij, t, s] + dlines['R (pu)'][ij]*model.I[ij,t,s] for ij in model.D if ij in model.N_D[i])
    pot_ess       = model.P_ESS_i_t_s[i, t, s] 
    pot_LS        = model.P_LS_i_t_s[i, t, s]
    pot_L_D       = model.L_D_i_t_s[i, t, s]
    
    return pot_geradores + pot_entrando - pot_saindo + pot_ess == pot_LS + pot_L_D

model.active_power_balance_distribution = pyo.Constraint(model.N_D, model.T, model.S, rule=active_power_balance_distribution_rule)

#solve the model


results = pyo.SolverFactory('gurobi',halt_on_ampl_error=True).solve(model,tee=True) 
model.display()
a=1




'''

"""
Eq. (3): Garante o balanço de potência reativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência e intercâmbio com o sistema de transmissão.
"""
def reactive_power_balance_distribution_rule(model, i, t, s):
    return (
        sum(model.Q_D_g_t_s[g, t, s] for g in model.G_i[i]) -
        sum(model.Q_ij_t_s[ij, t, s] for ij in model.D if ij in model.N_D[i]) +
        model.Q_SE_i_t_s[i, t, s]
        == model.L_i * (model.P_LS_i_t_s[i, t, s] + model.L_D_i_t_s[i, t, s])
    )
model.reactive_power_balance_distribution = pyo.Constraint(model.N_D, model.T, model.S, rule=reactive_power_balance_distribution_rule)


""" Eq. (4): Garante o balanço de potência reativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência reativa e intercâmbio com o sistema de transmissão."""

def reactive_power_balance_rule(model, i, t, s):
    return (
        sum(model.Q_D_g_t_s[g, t, s] for g in model.G_i[i]) +
        sum(model.Q_ki_t_s[ki, t, s] for ki in model.D if ki in model.N_D[i]) -
        sum(model.Q_ij_t_s[ij, t, s] + model.X_ij[ij] * model.ISQ_ij_t_s[ij, t, s] for ij in model.D if ij in model.N_D[i]) +
        model.Q_SE_i_t_s[i, t, s]
        == model.L_i * (model.P_LS_i_t_s[i, t, s] + model.L_D_i_t_s[i, t, s])
    )
model.reactive_power_balance = pyo.Constraint(model.N_D, model.T, model.S, rule=reactive_power_balance_rule)


""" Eq. (5): Modela a queda de tensão nas linhas de distribuição 
usando uma formulação de programação cônica de segunda ordem. """

def voltage_drop_rule(model, ij, t, s):
    return (
        model.V_SQ_i_t_s[ij[0], t, s] - model.V_SQ_i_t_s[ij[1], t, s] ==
        2 * (model.P_ij_t_s[ij, t, s] * model.R_ij[ij] + model.Q_ij_t_s[ij, t, s] * model.X_ij[ij]) +
        (model.R_ij[ij]**2 + model.X_ij[ij]**2) * model.ISQ_ij_t_s[ij, t, s]
    )
model.voltage_drop = pyo.Constraint(model.D, model.T, model.S, rule=voltage_drop_rule)


""" Eq. (6): Modela a queda de tensão considerando resistência, reatância 
e corrente quadrática nas linhas de distribuição. """

def voltage_drop_distribution_rule(model, ij, t, s):
    return (
        model.V_SQ_i_t_s[ij[0], t, s] - model.V_SQ_i_t_s[ij[1], t, s] ==
        2 * (model.P_ij_t_s[ij, t, s] * model.R_ij[ij] + model.Q_ij_t_s[ij, t, s] * model.X_ij[ij]) +
        (model.R_ij[ij]**2 + model.X_ij[ij]**2) * model.ISQ_ij_t_s[ij, t, s]
    )
model.voltage_drop_distribution = pyo.Constraint(model.D, model.T, model.S, rule=voltage_drop_distribution_rule)


"""
Eq. (7): Define os limites de corrente e tensão nas linhas de distribuição, 
garantindo a operação dentro de faixas seguras.
"""
def current_voltage_limit_distribution_rule(model, ij, t, s):
    return (
        model.P_ij_t_s[ij, t, s]**2 + model.Q_ij_t_s[ij, t, s]**2 <=
        model.ISQ_ij_t_s[ij, t, s] * model.V_SQ_i_t_s[ij[1], t, s]
    )
model.current_voltage_limit_distribution = pyo.Constraint(model.D, model.T, model.S, rule=current_voltage_limit_distribution_rule)


#################### Geradores ##################

"""
Eq. (8): Estabelece os limites de geração para geradores convencionais, 
considerando potência ativa e reativa.
"""
def conventional_generator_limit_rule(model, g, t, s):
    return (
        model.P_D_g_t_s[g, t, s]**2 + model.Q_D_g_t_s[g, t, s]**2 <= model.S_D_g[g]**2
    )
model.conventional_generator_limit = pyo.Constraint(
    [g for g in model.G_D if g == 'Conventional'], 
    model.T, model.S, 
    rule=conventional_generator_limit_rule
)


"""
Eq. (9): Define os limites de geração para geradores renováveis, 
considerando a disponibilidade de energia renovável.
"""
def renewable_generator_limit_rule(model, g, t, s):
    return (
        model.P_D_g_t_s[g, t, s]**2 + model.Q_D_g_t_s[g, t, s]**2 <= 
        model.r_g_t_s[g, t, s] * model.S_D_g[g]**2
    )
model.renewable_generator_limit = pyo.Constraint(
    [g for g in model.G_D if g == 'Renewable'], 
    model.T, model.S, 
    rule=renewable_generator_limit_rule
)


"""
Eq. (10): Garante que a soma dos blocos de oferta de geração 
para cada gerador em cada período e cenário seja igual à geração total.
"""
def generation_bid_block_balance_rule(model, g, t, s):
    return (
        sum(model.rho_D_g_b_t_s[g, b, t, s] for b in model.B_g[g]) == 
        model.P_D_g_t_s[g, t, s]
    )
model.generation_bid_block_balance = pyo.Constraint(
    model.G_D, model.T, model.S, 
    rule=generation_bid_block_balance_rule
)


"""
Eq. (11): Define o limite máximo para cada bloco de oferta de geração, 
impedindo que um bloco ultrapasse sua capacidade predefinida.
"""
def generation_bid_block_limit_rule(model, g, b, t, s):
    return (
        model.rho_D_g_b_t_s[g, b, t, s] <= model.P_D_g_b[g, b]
    )
model.generation_bid_block_limit = pyo.Constraint(
    model.G_D, model.B_g, model.T, model.S, 
    rule=generation_bid_block_limit_rule
)


"""
Eq. (12): Garante que a soma dos blocos de oferta aceitos para um gerador 
seja igual à sua geração total de potência ativa.
"""
def bid_block_generation_rule(model, g, t, s):
    return sum(model.rho_D_g_b_t_s[g, b, t, s] for b in b) == model.P_D_g_t_s[g, t, s]
model.bid_block_generation = pyo.Constraint(model.G_D, model.T, model.S, rule=bid_block_generation_rule)


"""
Eq. (13): Limita o bloco de oferta aceito ao tamanho máximo definido 
para aquele bloco específico.
"""
def bid_block_limit_rule(model, g, b, t, s):
    return model.rho_D_g_b_t_s[g, b, t, s] <= model.P_D_g_b[g, b]
model.bid_block_limit = pyo.Constraint(model.G_D, b, model.T, model.S, rule=bid_block_limit_rule)


"""
Eq. (14): Define os limites inferiores e superiores para o deslocamento 
de carga em cada nó.
"""
def load_shifting_limit_rule(model, i, t, s):
    return (-model.P_LS_i[i], model.P_LS_i_t_s[i, t, s], model.P_LS_i[i])
model.load_shifting_limit = pyo.Constraint(model.N_D, model.T, model.S, rule=load_shifting_limit_rule)


"""
Eq. (15): Garante que o deslocamento total de carga para cada nó ao 
longo de todos os períodos seja zero.
"""
def load_shifting_balance_rule(model, i, s):
    return sum(model.P_LS_i_t_s[i, t, s] for t in model.T) == 0
model.load_shifting_balance = pyo.Constraint(model.N_D, model.S, rule=load_shifting_balance_rule)


#################### Armazenamento ##################
"""
Eq. (16): Define os limites de carga e descarga para os sistemas 
de armazenamento de energia.
"""
def storage_power_limit_rule(model, i, t, s):
    return (-model.P_ESS_i[i], model.P_ESS_i_t_s[i, t, s], model.P_ESS_i[i])
model.storage_power_limit = pyo.Constraint(model.N_D, model.T, model.S, rule=storage_power_limit_rule)


"""
Eq. (17): Define os limites mínimo e máximo do Estado de Carga (SOC) 
dos sistemas de armazenamento de energia para cada nó.
"""
def ess_soc_limit_rule(model, i, t, s):
    return (
        model.SOC_ESS_MIN[i] <= model.SOC_i_t_s[i, t, s] <= 
        model.SOC_ESS_MAX[i]
    )
model.ess_soc_limit = pyo.Constraint(model.N_D, model.T, model.S, rule=ess_soc_limit_rule)


"""
Eq. (18-19): Calcula a evolução do Estado de Carga (SOC) 
dos sistemas de armazenamento ao longo do tempo.
"""
def ess_soc_update_rule(model, i, t, s):
    if t == model.T.first():
        return model.SOC_i_t_s[i, t, s] == model.SOC_i_0[i]
    else:
        return (
            model.SOC_i_t_s[i, t, s] == 
            model.SOC_i_t_s[i, t-1, s] + 
            model.PESS_i_t_s[i, t, s] * model.eta_charge - 
            model.PESS_i_t_s[i, t, s] / model.eta_discharge
        )
model.ess_soc_update = pyo.Constraint(model.N_D, model.T, model.S, rule=ess_soc_update_rule)


"""
Eq. (20): Garante que o Estado de Carga (SOC) final 
seja igual ao valor inicial para cada sistema de armazenamento.
"""
def ess_final_soc_balance_rule(model, i, s):
    return (
        model.SOC_i_t_s[i, model.T.last(), s] == model.SOC_i_0[i]
    )
model.ess_final_soc_balance = pyo.Constraint(model.N_D, model.S, rule=ess_final_soc_balance_rule)


#################### Carbono ##################

"""
Eq. (21): Calcula as emissões de carbono dos geradores 
conectados ao sistema de distribuição.
"""
def carbon_emission_distribution_rule(model, g, t, s):
    return (
        model.nu_D_g_t_s[g, t, s] == 
        model.gamma_D_g[g] * model.P_D_g_t_s[g, t, s]
    )
model.carbon_emission_distribution = pyo.Constraint(model.G_D, model.T, model.S, rule=carbon_emission_distribution_rule)


"""
Eq. (22): Estabelece o balanço de carbono para o sistema de distribuição.
"""
def carbon_balance_distribution_rule(model, t, s):
    return (
        sum(model.rho_D_g_t_s[g, t, s] for g in model.G_D) == 
        -model.rho_SE_t_s[t, s]
    )
model.carbon_balance_distribution = pyo.Constraint(model.T, model.S, rule=carbon_balance_distribution_rule)


"""
Eq. (23): Define o limite de emissões de carbono para o sistema de distribuição.
"""
def carbon_emission_limit_distribution_rule(model, t, s):
    return (
        sum(model.nu_D_g_t_s[g, t, s] + model.rho_D_g_t_s[g, t, s] for g in model.G_D) <= 
        model.Gamma_D_t
    )
model.carbon_emission_limit_distribution = pyo.Constraint(model.T, model.S, rule=carbon_emission_limit_distribution_rule)
'''