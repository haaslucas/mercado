import pyomo.environ as pyo
import pyomo.opt as opt
import pandas as pd

pasta = './01_modelo_pyomo/data/'
I_t = pd.read_csv(pasta + 'matriz_incidencia_transmissao.csv')
I_t = I_t.values.tolist()

dlines  = pd.read_csv(pasta + 'distribution_line_data.csv')
tnodes = pd.read_csv(pasta + 'transmission_node_data.csv')
dnodes = pd.read_csv(pasta + 'distribution_node_data.csv')

g_t_nd = pd.read_csv(pasta + 'transmission_connected_non_dispatchable_generators.csv')
g_t_d  = pd.read_csv(pasta + 'transmission_connected_dispatchable_generators.csv')
g_d_nd = pd.read_csv(pasta + 'dso_connected_non_dispatchable_generators.csv')
g_d_d  = pd.read_csv(pasta + 'dso_connected_dispatchable_generators.csv')

SBASE = 100               # MVA
VBASE = 24.9              # kV
IBASE = 100 * 1e3 / 24.90 # A
SE_Capacity = 50          # MWh
U=1                       # Rigidez da regulação de emissões de carbono   
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
s = range(0,5)  # Índice para cenários  
t = range(0,24) # Índice para períodos de tempo  

#################### Conjuntos ####################
um.B_g = pyo.Set()          # Set of bids provided by generator g
um.G_D = pyo.Set()          # Set of distribution-connected generators
um.G_i = pyo.Set()          # Set of generators connected to node i
um.G_T = pyo.Set()          # Set of transmission-connected generators
um.D = pyo.Set()            # Set of distribution lines
um.T = pyo.Set()            # Set of transmission lines

um.d_nodes = pyo.RangeSet(1, 34) # Set of distribution nodes
um.d_lines = pyo.RangeSet(1, 33)          # Set of distribution lines

um.N_T = pyo.Set()          # Set of transmission nodes
um.N_INF = pyo.Set()        # Set of interface nodes
um.S = pyo.Set()            # Set of scenarios
um.T = pyo.Set()            # Set of time periods


#################### Parametros ####################

# Parâmetros relacionados a geradores
um.A_i_l = pyo.Param(i,l,initialize= lambda m,col,lin: I_t[lin][col])# Incidence matrix of transmission lines

# Parâmetros relacionados a lances de geradores
um.C_D_g_b = pyo.Param(g, b, initialize=lambda m, g, b: 1)           # Distribution-connected generators' bids
um.C_T_g_b = pyo.Param(g, b, initialize=lambda m, g, b: 1)           # Transmission-connected generators' bids
um.P_D_g_b = pyo.Param(g, b, initialize=lambda m, g, b: 1)           # Size of the generator's electricity block bid
um.P_T_g_b = pyo.Param(g, b, initialize=lambda m, g, b: 1)           # Size of the generator's electricity block bid

# Parâmetros relacionados a geradores
um.S_D_g = pyo.Param(g, initialize=lambda m, g: 1)                   # Maximum power injection of generators
um.S_T_g = pyo.Param(g, initialize=lambda m, g: 1)                   # Maximum power injection of generators
um.gamma_D_g = pyo.Param(g, initialize=lambda m, g: 1)               # Generator's carbon intensity
um.gamma_T_g = pyo.Param(g, initialize=lambda m, g: 1)               # Generator's carbon intensity
um.nu_D_g_t_s = pyo.Param(g, t, s, initialize=lambda m, g, t, s: 1)  # Generator's carbon emission
um.nu_T_g_t_s = pyo.Param(g, t, s, initialize=lambda m, g, t, s: 1)  # Generator's carbon emission

# Parâmetros relacionados a demanda e disponibilidade
um.L_D_i_t_s = pyo.Param(i, t, s, initialize=lambda m, i, t, s: 1)   # Power demand of a distribution node
um.L_T_i_t_s = pyo.Param(i, t, s, initialize=lambda m, i, t, s: 1)   # Power demand of a transmission node
um.r_g_t_s = pyo.Param(g, t, s, initialize=lambda m, g, t, s: 1)     # Renewable-based power availability

# Parâmetros relacionados a linhas e capacidades
um.I_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Capacity of distribution lines
um.PF_l = pyo.Param(l, initialize=lambda m, l: 1)                    # Power flow capacity of transmission lines
um.R_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Resistance of distribution lines
um.X_l = pyo.Param(l, initialize=lambda m, l: 1)                     # Reactance of transmission lines
um.X_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Reactance of distribution lines


# Parâmetros de mercado e regulação
um.p_s_t = pyo.Param(s, t, initialize=lambda m, s, t: 1)             # Weight of scenario s at period t
um.pi_E_i_t_s = pyo.Param(i, t, s, initialize=lambda m, i, t, s: 1)  # Bids on the electricity market
um.pi_C_t_s = pyo.Param(t, s, initialize=lambda m, t, s: 1)          # Bids on the carbon market
um.U = pyo.Param(initialize=lambda m: 1)                             # Rigidity of the carbon emission regulation
um.Lambda = pyo.Param(initialize=lambda m: 1)                        # Average carbon intensity of the entire power system

# Parâmetros das linhas de distribuição dlines = distribution_line_data.csv
um.R = pyo.Param(i, initialize=lambda m,i: dlines['R (pu)'])                 # branch resistance (ohm)
um.X = pyo.Param(i, initialize=lambda m,i: dlines['X (pu)'])                 # branch reactance (ohm)
#um.Z = pyo.Param(i, initialize=lambda m,i: dlines['X (pu)'])                 # branch reactance (ohm)
um.Imax = pyo.Param(i, initialize=lambda m,i: dlines['Line Ampacity (A)'])   # branch maximum current (kA)

# Parâmetros relacionados a BESS do sistema de distribuição ESS
um.ESS_P_max_charge = pyo.Param(i, initialize=lambda m,i: dr['Max ESS Charge (MW)']) # Max charging power of BESS
um.ESS_P_max_discha = pyo.Param(i, initialize=lambda m,i: dr['Max ESS Discharge (MW)']) # Max discharging power of BESS
um.ESS_SOC_min = pyo.Param(i, initialize=lambda m,i: dr['Min ESS SOC (MWh)']) # Min SOC of BESS
um.ESS_SOC_max = pyo.Param(i, initialize=lambda m,i: dr['Max ESS SOC (MWh)']) # Max SOC of BESS
um.ESS_SOC_0   = pyo.Param(i, initialize=lambda m,i: dr['Ini ESS SOC (MWh)']) # Initial SOC of BESS

um.P_LS_i = pyo.Param(i, initialize=lambda m, i: 1)                  # Maximum nodal load increase/reduction


#################### Variáveis de decisão ####################

# Variáveis relacionadas à distribuição

um.P_D_g_t_s = pyo.Var(um.G_D, um.T, um.S) # Active power dispatch of distribution-connected generators
um.Q_D_g_t_s = pyo.Var()   # Generator reactive power dispatch

um.P_ij_t_s = pyo.Var(um.d_lines, um.d_lines, um.T, um.S) # Active power flow on distribution lines

um.Q_ij_t_s = pyo.Var()    # Reactive power flow on distribution lines
um.P_ESS_i_t_s = pyo.Var() # Power charge/discharge of storage systems
um.P_LS_i_t_s = pyo.Var()  # Nodal load shift
um.V_SQ_i_t_s = pyo.Var()  # Nodal voltage magnitude squared
um.SOC_i_t_s = pyo.Var()   # State of charge of energy storage systems
um.rho_D_g_b_t_s = pyo.Var() # Accepted bid of block b
um.upsilon_D_g_t_s = pyo.Var() # Carbon allowances absorbed/provided by generator g

# Variáveis relacionadas à transmissão
um.P_T_g_t_s = pyo.Var()    # Generator's total power dispatch
um.PF_l_t_s = pyo.Var()     # Power flow on transmission lines
um.P_SE_i_t_s = pyo.Var()   # Active power trades in the wholesale market
um.Q_SE_i_t_s = pyo.Var()   # Reactive power trades in the wholesale market
um.upsilon_SE_t_s = pyo.Var() # Carbon allowance trades in the wholesale market
um.delta_i_t_s = pyo.Var()  # Voltage angle of transmission nodes
um.nu_T_g_t_s = pyo.Var()   # Generator's carbon emission
um.rho_T_g_b_t_s = pyo.Var() # Accepted bid of block b
um.upsilon_T_g_t_s = pyo.Var() # Carbon allowances absorbed/provided by generator g


# Função Objetivo para o DSO (Upper-level problem)
um.objective = pyo.Objective(
    expr=sum(
        um.p_s_t * (
            # Custos de geração dos geradores conectados à rede de distribuição
            sum(
                um.C_D_g_b[g, b] * um.rho_D_g_b_t_s[g, b, t, s] 
                for g in g 
                for b in b
            ) +
            # Receita de troca de allowances de carbono com o ISO
            um.upsilon_SE_t_s[t, s] +
            # Receita de intercâmbio de energia com o ISO nos nós de interface
            sum(
                um.P_SE_i_t_s[i, t, s] 
                for i in um.N_INF
            )
        )
        for s in s 
        for t in t
    ),
    sense=pyo.minimize
)

####################################################
#################### RESTRIÇÕES ####################
####################################################


#################### Distribuição ##################

"""
Eq. (2): Garante o balanço de potência ativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência, armazenamento e intercâmbio com o sistema de transmissão.
"""
def active_power_balance_distribution_rule(um, n, m, t, s):
    pot_geradores = sum(um.P_D_g_t_s[g, t, s] for g in g_d_d['Node'] if g in um.dist_nodes[n])
    pot_entrando = sum(um.P_ij_t_s[n, m, t, s] + dlines['R (pu)'][n,m]*um.I[n, m, t, s])
    
      
    return pot_geradores - pot_entrando
um.active_power_balance_distribution = pyo.Constraint(um.dist_nodes, um.T, um.S, rule=active_power_balance_distribution_rule)


"""
Eq. (2): Garante o balanço de potência ativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência, armazenamento e intercâmbio com o sistema de transmissão.
"""
def active_power_balance_distribution_rule(model, i, t, s):
    return (
        sum(um.P_D_g_t_s[g, t, s] for g in um.G_i[i]) -
        sum(um.P_ij_t_s[ij, t, s] for ij in um.D if ij in um.N_D[i]) +
        um.P_ESS_i_t_s[i, t, s] + um.P_SE_i_t_s[i, t, s]
        == um.P_LS_i_t_s[i, t, s] + um.L_D_i_t_s[i, t, s]
    )
um.active_power_balance_distribution = pyo.Constraint(um.N_D, um.T, um.S, rule=active_power_balance_distribution_rule)


 
"""
Eq. (3): Garante o balanço de potência reativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência e intercâmbio com o sistema de transmissão.
"""
def reactive_power_balance_distribution_rule(model, i, t, s):
    return (
        sum(um.Q_D_g_t_s[g, t, s] for g in um.G_i[i]) -
        sum(um.Q_ij_t_s[ij, t, s] for ij in um.D if ij in um.dist_nodes[i]) +
        um.Q_SE_i_t_s[i, t, s]
        == um.L_i * (um.P_LS_i_t_s[i, t, s] + um.L_D_i_t_s[i, t, s])
    )
um.reactive_power_balance_distribution = pyo.Constraint(um.dist_nodes, um.T, um.S, rule=reactive_power_balance_distribution_rule)


""" Eq. (4): Garante o balanço de potência reativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência reativa e intercâmbio com o sistema de transmissão."""

def reactive_power_balance_rule(model, i, t, s):
    return (
        sum(um.Q_D_g_t_s[g, t, s] for g in um.G_i[i]) +
        sum(um.Q_ki_t_s[ki, t, s] for ki in um.D if ki in um.dist_nodes[i]) -
        sum(um.Q_ij_t_s[ij, t, s] + um.X_ij[ij] * um.ISQ_ij_t_s[ij, t, s] for ij in um.D if ij in um.dist_nodes[i]) +
        um.Q_SE_i_t_s[i, t, s]
        == um.L_i * (um.P_LS_i_t_s[i, t, s] + um.L_D_i_t_s[i, t, s])
    )
um.reactive_power_balance = pyo.Constraint(um.dist_nodes, um.T, um.S, rule=reactive_power_balance_rule)


""" Eq. (5): Modela a queda de tensão nas linhas de distribuição 
usando uma formulação de programação cônica de segunda ordem. """

def voltage_drop_rule(model, ij, t, s):
    return (
        um.V_SQ_i_t_s[ij[0], t, s] - um.V_SQ_i_t_s[ij[1], t, s] ==
        2 * (um.P_ij_t_s[ij, t, s] * um.R_ij[ij] + um.Q_ij_t_s[ij, t, s] * um.X_ij[ij]) +
        (um.R_ij[ij]**2 + um.X_ij[ij]**2) * um.ISQ_ij_t_s[ij, t, s]
    )
um.voltage_drop = pyo.Constraint(um.D, um.T, um.S, rule=voltage_drop_rule)


""" Eq. (6): Modela a queda de tensão considerando resistência, reatância 
e corrente quadrática nas linhas de distribuição. """

def voltage_drop_distribution_rule(model, ij, t, s):
    return (
        um.V_SQ_i_t_s[ij[0], t, s] - um.V_SQ_i_t_s[ij[1], t, s] ==
        2 * (um.P_ij_t_s[ij, t, s] * um.R_ij[ij] + um.Q_ij_t_s[ij, t, s] * um.X_ij[ij]) +
        (um.R_ij[ij]**2 + um.X_ij[ij]**2) * um.ISQ_ij_t_s[ij, t, s]
    )
um.voltage_drop_distribution = pyo.Constraint(um.D, um.T, um.S, rule=voltage_drop_distribution_rule)


"""
Eq. (7): Define os limites de corrente e tensão nas linhas de distribuição, 
garantindo a operação dentro de faixas seguras.
"""
def current_voltage_limit_distribution_rule(model, ij, t, s):
    return (
        um.P_ij_t_s[ij, t, s]**2 + um.Q_ij_t_s[ij, t, s]**2 <=
        um.ISQ_ij_t_s[ij, t, s] * um.V_SQ_i_t_s[ij[1], t, s]
    )
um.current_voltage_limit_distribution = pyo.Constraint(um.D, um.T, um.S, rule=current_voltage_limit_distribution_rule)


#################### Geradores ##################

"""
Eq. (8): Estabelece os limites de geração para geradores convencionais, 
considerando potência ativa e reativa.
"""
def conventional_generator_limit_rule(model, g, t, s):
    return (
        um.P_D_g_t_s[g, t, s]**2 + um.Q_D_g_t_s[g, t, s]**2 <= um.S_D_g[g]**2
    )
um.conventional_generator_limit = pyo.Constraint(
    [g for g in um.G_D if g == 'Conventional'], 
    um.T, um.S, 
    rule=conventional_generator_limit_rule
)


"""
Eq. (9): Define os limites de geração para geradores renováveis, 
considerando a disponibilidade de energia renovável.
"""
def renewable_generator_limit_rule(model, g, t, s):
    return (
        um.P_D_g_t_s[g, t, s]**2 + um.Q_D_g_t_s[g, t, s]**2 <= 
        um.r_g_t_s[g, t, s] * um.S_D_g[g]**2
    )
um.renewable_generator_limit = pyo.Constraint(
    [g for g in um.G_D if g == 'Renewable'], 
    um.T, um.S, 
    rule=renewable_generator_limit_rule
)


"""
Eq. (10): Garante que a soma dos blocos de oferta de geração 
para cada gerador em cada período e cenário seja igual à geração total.
"""
def generation_bid_block_balance_rule(model, g, t, s):
    return (
        sum(um.rho_D_g_b_t_s[g, b, t, s] for b in um.B_g[g]) == 
        um.P_D_g_t_s[g, t, s]
    )
um.generation_bid_block_balance = pyo.Constraint(
    um.G_D, um.T, um.S, 
    rule=generation_bid_block_balance_rule
)


"""
Eq. (11): Define o limite máximo para cada bloco de oferta de geração, 
impedindo que um bloco ultrapasse sua capacidade predefinida.
"""
def generation_bid_block_limit_rule(model, g, b, t, s):
    return (
        um.rho_D_g_b_t_s[g, b, t, s] <= um.P_D_g_b[g, b]
    )
um.generation_bid_block_limit = pyo.Constraint(
    um.G_D, um.B_g, um.T, um.S, 
    rule=generation_bid_block_limit_rule
)


"""
Eq. (12): Garante que a soma dos blocos de oferta aceitos para um gerador 
seja igual à sua geração total de potência ativa.
"""
def bid_block_generation_rule(model, g, t, s):
    return sum(um.rho_D_g_b_t_s[g, b, t, s] for b in b) == um.P_D_g_t_s[g, t, s]
um.bid_block_generation = pyo.Constraint(um.G_D, um.T, um.S, rule=bid_block_generation_rule)


"""
Eq. (13): Limita o bloco de oferta aceito ao tamanho máximo definido 
para aquele bloco específico.
"""
def bid_block_limit_rule(model, g, b, t, s):
    return um.rho_D_g_b_t_s[g, b, t, s] <= um.P_D_g_b[g, b]
um.bid_block_limit = pyo.Constraint(um.G_D, b, um.T, um.S, rule=bid_block_limit_rule)


"""
Eq. (14): Define os limites inferiores e superiores para o deslocamento 
de carga em cada nó.
"""
def load_shifting_limit_rule(model, i, t, s):
    return (-um.P_LS_i[i], um.P_LS_i_t_s[i, t, s], um.P_LS_i[i])
um.load_shifting_limit = pyo.Constraint(um.dist_nodes, um.T, um.S, rule=load_shifting_limit_rule)


"""
Eq. (15): Garante que o deslocamento total de carga para cada nó ao 
longo de todos os períodos seja zero.
"""
def load_shifting_balance_rule(model, i, s):
    return sum(um.P_LS_i_t_s[i, t, s] for t in um.T) == 0
um.load_shifting_balance = pyo.Constraint(um.dist_nodes, um.S, rule=load_shifting_balance_rule)


#################### Armazenamento ##################
"""
Eq. (16): Define os limites de carga e descarga para os sistemas 
de armazenamento de energia.
"""
def storage_power_limit_rule(model, i, t, s):
    return (-um.P_ESS_i[i], um.P_ESS_i_t_s[i, t, s], um.P_ESS_i[i])
um.storage_power_limit = pyo.Constraint(um.dist_nodes, um.T, um.S, rule=storage_power_limit_rule)


"""
Eq. (17): Define os limites mínimo e máximo do Estado de Carga (SOC) 
dos sistemas de armazenamento de energia para cada nó.
"""
def ess_soc_limit_rule(model, i, t, s):
    return (
        um.SOC_ESS_MIN[i] <= um.SOC_i_t_s[i, t, s] <= 
        um.SOC_ESS_MAX[i]
    )
um.ess_soc_limit = pyo.Constraint(um.dist_nodes, um.T, um.S, rule=ess_soc_limit_rule)


"""
Eq. (18-19): Calcula a evolução do Estado de Carga (SOC) 
dos sistemas de armazenamento ao longo do tempo.
"""
def ess_soc_update_rule(model, i, t, s):
    if t == um.T.first():
        return um.SOC_i_t_s[i, t, s] == um.SOC_i_0[i]
    else:
        return (
            um.SOC_i_t_s[i, t, s] == 
            um.SOC_i_t_s[i, t-1, s] + 
            um.PESS_i_t_s[i, t, s] * um.eta_charge - 
            um.PESS_i_t_s[i, t, s] / um.eta_discharge
        )
um.ess_soc_update = pyo.Constraint(um.dist_nodes, um.T, um.S, rule=ess_soc_update_rule)


"""
Eq. (20): Garante que o Estado de Carga (SOC) final 
seja igual ao valor inicial para cada sistema de armazenamento.
"""
def ess_final_soc_balance_rule(model, i, s):
    return (
        um.SOC_i_t_s[i, um.T.last(), s] == um.SOC_i_0[i]
    )
um.ess_final_soc_balance = pyo.Constraint(um.dist_nodes, um.S, rule=ess_final_soc_balance_rule)


#################### Carbono ##################

"""
Eq. (21): Calcula as emissões de carbono dos geradores 
conectados ao sistema de distribuição.
"""
def carbon_emission_distribution_rule(model, g, t, s):
    return (
        um.nu_D_g_t_s[g, t, s] == 
        um.gamma_D_g[g] * um.P_D_g_t_s[g, t, s]
    )
um.carbon_emission_distribution = pyo.Constraint(um.G_D, um.T, um.S, rule=carbon_emission_distribution_rule)


"""
Eq. (22): Estabelece o balanço de carbono para o sistema de distribuição.
"""
def carbon_balance_distribution_rule(model, t, s):
    return (
        sum(um.rho_D_g_t_s[g, t, s] for g in um.G_D) == 
        -um.rho_SE_t_s[t, s]
    )
um.carbon_balance_distribution = pyo.Constraint(um.T, um.S, rule=carbon_balance_distribution_rule)


"""
Eq. (23): Define o limite de emissões de carbono para o sistema de distribuição.
"""
def carbon_emission_limit_distribution_rule(model, t, s):
    return (
        sum(um.nu_D_g_t_s[g, t, s] + um.rho_D_g_t_s[g, t, s] for g in um.G_D) <= 
        um.Gamma_D_t
    )
um.carbon_emission_limit_distribution = pyo.Constraint(um.T, um.S, rule=carbon_emission_limit_distribution_rule)
