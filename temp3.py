import pyomo.environ as pyo
import pyomo.opt as opt
import pandas as pd

I_t = pd.read_csv('data/matriz_incidencia_transmissao.csv')
I_t = I_t.values.tolist()

SBASE = 100 #MVA
VBASE = 24.9 #kV

#nó de conexão entre a rede de distribuição e a rede de transmissão
N_INF_d = 1
N_INF_t = 5

Umodel, LLmodel = pyo.ConcreteModel(), pyo.ConcreteModel()

#################### Índices ####################   
b = range(0,5)  # Índice para bloco de oferta  
g = range(0,8)  # Índice para geradores  
i, j = range(0,14), range(0,14)  # Índice para nós  
ij, ki = range(0,40), range(0,40)  # Índice para linhas de distribuição  
l = range(0,20)  # Índice para linhas de transmissão  
s = range(0,5)  # Índice para cenários  
t = range(0,24) # Índice para períodos de tempo  

#################### Conjuntos ####################
Umodel.B_g = pyo.Set()          # Set of bids provided by generator g
Umodel.G_D = pyo.Set()          # Set of distribution-connected generators
Umodel.G_i = pyo.Set()          # Set of generators connected to node i
Umodel.G_T = pyo.Set()          # Set of transmission-connected generators
Umodel.D = pyo.Set()            # Set of distribution lines
Umodel.T = pyo.Set()            # Set of transmission lines
Umodel.N_D = pyo.Set()          # Set of distribution nodes
Umodel.N_T = pyo.Set()          # Set of transmission nodes
Umodel.N_INF = pyo.Set()        # Set of interface nodes
Umodel.S = pyo.Set()            # Set of scenarios
Umodel.T = pyo.Set()            # Set of time periods


#################### Parametros ####################

# Parâmetros relacionados a geradores
Umodel.A_i_l = pyo.Param(i,l,initialize= lambda m,col,lin: I_t[lin][col])# Incidence matrix of transmission lines
Umodel.C_D_g_b = pyo.Param(g, b, initialize=lambda m, g, b: 1)           # Distribution-connected generators' bids
Umodel.C_T_g_b = pyo.Param(g, b, initialize=lambda m, g, b: 1)           # Transmission-connected generators' bids
Umodel.P_D_g_b = pyo.Param(g, b, initialize=lambda m, g, b: 1)           # Size of the generator's electricity block bid
Umodel.P_T_g_b = pyo.Param(g, b, initialize=lambda m, g, b: 1)           # Size of the generator's electricity block bid
Umodel.S_D_g = pyo.Param(g, initialize=lambda m, g: 1)                   # Maximum power injection of generators
Umodel.S_T_g = pyo.Param(g, initialize=lambda m, g: 1)                   # Maximum power injection of generators
Umodel.gamma_D_g = pyo.Param(g, initialize=lambda m, g: 1)               # Generator's carbon intensity
Umodel.gamma_T_g = pyo.Param(g, initialize=lambda m, g: 1)               # Generator's carbon intensity
Umodel.nu_D_g_t_s = pyo.Param(g, t, s, initialize=lambda m, g, t, s: 1)  # Generator's carbon emission
Umodel.nu_T_g_t_s = pyo.Param(g, t, s, initialize=lambda m, g, t, s: 1)  # Generator's carbon emission

# Parâmetros relacionados a demanda e disponibilidade
Umodel.L_D_i_t_s = pyo.Param(i, t, s, initialize=lambda m, i, t, s: 1)   # Power demand of a distribution node
Umodel.L_T_i_t_s = pyo.Param(i, t, s, initialize=lambda m, i, t, s: 1)   # Power demand of a transmission node
Umodel.r_g_t_s = pyo.Param(g, t, s, initialize=lambda m, g, t, s: 1)     # Renewable-based power availability

# Parâmetros relacionados a linhas e capacidades
Umodel.I_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Capacity of distribution lines
Umodel.PF_l = pyo.Param(l, initialize=lambda m, l: 1)                    # Power flow capacity of transmission lines
Umodel.R_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Resistance of distribution lines
Umodel.X_l = pyo.Param(l, initialize=lambda m, l: 1)                     # Reactance of transmission lines
Umodel.X_ij = pyo.Param(ij, initialize=lambda m, ij: 1)                  # Reactance of distribution lines

# Parâmetros relacionados a armazenamento e flexibilidade
Umodel.P_ESS_i = pyo.Param(i, initialize=lambda m, i: 1)                 # Maximum charging and discharging power of storage systems
Umodel.SOC_ESS_i = pyo.Param(i, initialize=lambda m, i: 1)               # Maximum and minimum state of charge of the storage systems
Umodel.SOC_0_i = pyo.Param(i, initialize=lambda m, i: 1)                 # Initial state of charge of the storage systems
Umodel.P_LS_i = pyo.Param(i, initialize=lambda m, i: 1)                  # Maximum nodal load increase/reduction

# Parâmetros de mercado e regulação
Umodel.p_s_t = pyo.Param(s, t, initialize=lambda m, s, t: 1)             # Weight of scenario s at period t
Umodel.pi_E_i_t_s = pyo.Param(i, t, s, initialize=lambda m, i, t, s: 1)  # Bids on the electricity market
Umodel.pi_C_t_s = pyo.Param(t, s, initialize=lambda m, t, s: 1)          # Bids on the carbon market
Umodel.U = pyo.Param(initialize=lambda m: 1)                             # Rigidity of the carbon emission regulation
Umodel.Lambda = pyo.Param(initialize=lambda m: 1)                        # Average carbon intensity of the entire power system



#################### Variáveis de decisão ####################

# Variáveis relacionadas à distribuição
Umodel.P_D_g_t_s = pyo.Var()   # Generator active power dispatch
Umodel.Q_D_g_t_s = pyo.Var()   # Generator reactive power dispatch
Umodel.P_ij_t_s = pyo.Var()    # Active power flow on distribution lines
Umodel.Q_ij_t_s = pyo.Var()    # Reactive power flow on distribution lines
Umodel.P_ESS_i_t_s = pyo.Var() # Power charge/discharge of storage systems
Umodel.P_LS_i_t_s = pyo.Var()  # Nodal load shift
Umodel.V_SQ_i_t_s = pyo.Var()  # Nodal voltage magnitude squared
Umodel.SOC_i_t_s = pyo.Var()   # State of charge of energy storage systems
Umodel.rho_D_g_b_t_s = pyo.Var() # Accepted bid of block b
Umodel.upsilon_D_g_t_s = pyo.Var() # Carbon allowances absorbed/provided by generator g

# Variáveis relacionadas à transmissão
Umodel.P_T_g_t_s = pyo.Var()    # Generator's total power dispatch
Umodel.PF_l_t_s = pyo.Var()     # Power flow on transmission lines
Umodel.P_SE_i_t_s = pyo.Var()   # Active power trades in the wholesale market
Umodel.Q_SE_i_t_s = pyo.Var()   # Reactive power trades in the wholesale market
Umodel.upsilon_SE_t_s = pyo.Var() # Carbon allowance trades in the wholesale market
Umodel.delta_i_t_s = pyo.Var()  # Voltage angle of transmission nodes
Umodel.nu_T_g_t_s = pyo.Var()   # Generator's carbon emission
Umodel.rho_T_g_b_t_s = pyo.Var() # Accepted bid of block b
Umodel.upsilon_T_g_t_s = pyo.Var() # Carbon allowances absorbed/provided by generator g


# Função Objetivo para o DSO (Upper-level problem)
Umodel.objective = pyo.Objective(
    expr=sum(
        Umodel.p_s_t * (
            # Custos de geração dos geradores conectados à rede de distribuição
            sum(
                Umodel.C_D_g_b[g, b] * Umodel.rho_D_g_b_t_s[g, b, t, s] 
                for g in g 
                for b in b
            ) +
            # Receita de troca de allowances de carbono com o ISO
            Umodel.upsilon_SE_t_s[t, s] +
            # Receita de intercâmbio de energia com o ISO nos nós de interface
            sum(
                Umodel.P_SE_i_t_s[i, t, s] 
                for i in Umodel.N_INF
            )
        )
        for s in s 
        for t in t
    ),
    sense=pyo.minimize
)


#################### Restrições ####################

"""
Eq. (2): Garante o balanço de potência ativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência, armazenamento e intercâmbio com o sistema de transmissão.
"""
def active_power_balance_distribution_rule(model, i, t, s):
    return (
        sum(Umodel.P_D_g_t_s[g, t, s] for g in Umodel.G_i[i]) -
        sum(Umodel.P_ij_t_s[ij, t, s] for ij in Umodel.D if ij in Umodel.N_D[i]) +
        Umodel.P_ESS_i_t_s[i, t, s] + Umodel.P_SE_i_t_s[i, t, s]
        == Umodel.P_LS_i_t_s[i, t, s] + Umodel.L_D_i_t_s[i, t, s]
    )
Umodel.active_power_balance_distribution = pyo.Constraint(Umodel.N_D, Umodel.T, Umodel.S, rule=active_power_balance_distribution_rule)


"""
Eq. (3): Garante o balanço de potência reativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência e intercâmbio com o sistema de transmissão.
"""
def reactive_power_balance_distribution_rule(model, i, t, s):
    return (
        sum(Umodel.Q_D_g_t_s[g, t, s] for g in Umodel.G_i[i]) -
        sum(Umodel.Q_ij_t_s[ij, t, s] for ij in Umodel.D if ij in Umodel.N_D[i]) +
        Umodel.Q_SE_i_t_s[i, t, s]
        == Umodel.L_i * (Umodel.P_LS_i_t_s[i, t, s] + Umodel.L_D_i_t_s[i, t, s])
    )
Umodel.reactive_power_balance_distribution = pyo.Constraint(Umodel.N_D, Umodel.T, Umodel.S, rule=reactive_power_balance_distribution_rule)


""" Eq. (4): Garante o balanço de potência reativa nos nós da rede de distribuição, 
considerando geração, fluxo de potência reativa e intercâmbio com o sistema de transmissão."""

def reactive_power_balance_rule(model, i, t, s):
    return (
        sum(Umodel.Q_D_g_t_s[g, t, s] for g in Umodel.G_i[i]) +
        sum(Umodel.Q_ki_t_s[ki, t, s] for ki in Umodel.D if ki in Umodel.N_D[i]) -
        sum(Umodel.Q_ij_t_s[ij, t, s] + Umodel.X_ij[ij] * Umodel.ISQ_ij_t_s[ij, t, s] for ij in Umodel.D if ij in Umodel.N_D[i]) +
        Umodel.Q_SE_i_t_s[i, t, s]
        == Umodel.L_i * (Umodel.P_LS_i_t_s[i, t, s] + Umodel.L_D_i_t_s[i, t, s])
    )
Umodel.reactive_power_balance = pyo.Constraint(Umodel.N_D, Umodel.T, Umodel.S, rule=reactive_power_balance_rule)


""" Eq. (5): Modela a queda de tensão nas linhas de distribuição 
usando uma formulação de programação cônica de segunda ordem. """

def voltage_drop_rule(model, ij, t, s):
    return (
        Umodel.V_SQ_i_t_s[ij[0], t, s] - Umodel.V_SQ_i_t_s[ij[1], t, s] ==
        2 * (Umodel.P_ij_t_s[ij, t, s] * Umodel.R_ij[ij] + Umodel.Q_ij_t_s[ij, t, s] * Umodel.X_ij[ij]) +
        (Umodel.R_ij[ij]**2 + Umodel.X_ij[ij]**2) * Umodel.ISQ_ij_t_s[ij, t, s]
    )
Umodel.voltage_drop = pyo.Constraint(Umodel.D, Umodel.T, Umodel.S, rule=voltage_drop_rule)


""" Eq. (6): Modela a queda de tensão considerando resistência, reatância 
e corrente quadrática nas linhas de distribuição. """

def voltage_drop_distribution_rule(model, ij, t, s):
    return (
        Umodel.V_SQ_i_t_s[ij[0], t, s] - Umodel.V_SQ_i_t_s[ij[1], t, s] ==
        2 * (Umodel.P_ij_t_s[ij, t, s] * Umodel.R_ij[ij] + Umodel.Q_ij_t_s[ij, t, s] * Umodel.X_ij[ij]) +
        (Umodel.R_ij[ij]**2 + Umodel.X_ij[ij]**2) * Umodel.ISQ_ij_t_s[ij, t, s]
    )
Umodel.voltage_drop_distribution = pyo.Constraint(Umodel.D, Umodel.T, Umodel.S, rule=voltage_drop_distribution_rule)


"""
Eq. (7): Define os limites de corrente e tensão nas linhas de distribuição, 
garantindo a operação dentro de faixas seguras.
"""
def current_voltage_limit_distribution_rule(model, ij, t, s):
    return (
        Umodel.P_ij_t_s[ij, t, s]**2 + Umodel.Q_ij_t_s[ij, t, s]**2 <=
        Umodel.ISQ_ij_t_s[ij, t, s] * Umodel.V_SQ_i_t_s[ij[1], t, s]
    )
Umodel.current_voltage_limit_distribution = pyo.Constraint(Umodel.D, Umodel.T, Umodel.S, rule=current_voltage_limit_distribution_rule)


"""
Eq. (8): Estabelece os limites de geração para geradores convencionais, 
considerando potência ativa e reativa.
"""
def conventional_generator_limit_rule(model, g, t, s):
    return (
        Umodel.P_D_g_t_s[g, t, s]**2 + Umodel.Q_D_g_t_s[g, t, s]**2 <= Umodel.S_D_g[g]**2
    )
Umodel.conventional_generator_limit = pyo.Constraint(
    [g for g in Umodel.G_D if g == 'Conventional'], 
    Umodel.T, Umodel.S, 
    rule=conventional_generator_limit_rule
)


"""
Eq. (9): Define os limites de geração para geradores renováveis, 
considerando a disponibilidade de energia renovável.
"""
def renewable_generator_limit_rule(model, g, t, s):
    return (
        Umodel.P_D_g_t_s[g, t, s]**2 + Umodel.Q_D_g_t_s[g, t, s]**2 <= 
        Umodel.r_g_t_s[g, t, s] * Umodel.S_D_g[g]**2
    )
Umodel.renewable_generator_limit = pyo.Constraint(
    [g for g in Umodel.G_D if g == 'Renewable'], 
    Umodel.T, Umodel.S, 
    rule=renewable_generator_limit_rule
)


"""
Eq. (10): Garante que a soma dos blocos de oferta de geração 
para cada gerador em cada período e cenário seja igual à geração total.
"""
def generation_bid_block_balance_rule(model, g, t, s):
    return (
        sum(Umodel.rho_D_g_b_t_s[g, b, t, s] for b in Umodel.B_g[g]) == 
        Umodel.P_D_g_t_s[g, t, s]
    )
Umodel.generation_bid_block_balance = pyo.Constraint(
    Umodel.G_D, Umodel.T, Umodel.S, 
    rule=generation_bid_block_balance_rule
)


"""
Eq. (11): Define o limite máximo para cada bloco de oferta de geração, 
impedindo que um bloco ultrapasse sua capacidade predefinida.
"""
def generation_bid_block_limit_rule(model, g, b, t, s):
    return (
        Umodel.rho_D_g_b_t_s[g, b, t, s] <= Umodel.P_D_g_b[g, b]
    )
Umodel.generation_bid_block_limit = pyo.Constraint(
    Umodel.G_D, Umodel.B_g, Umodel.T, Umodel.S, 
    rule=generation_bid_block_limit_rule
)


"""
Eq. (12): Garante que a soma dos blocos de oferta aceitos para um gerador 
seja igual à sua geração total de potência ativa.
"""
def bid_block_generation_rule(model, g, t, s):
    return sum(Umodel.rho_D_g_b_t_s[g, b, t, s] for b in b) == Umodel.P_D_g_t_s[g, t, s]
Umodel.bid_block_generation = pyo.Constraint(Umodel.G_D, Umodel.T, Umodel.S, rule=bid_block_generation_rule)


"""
Eq. (13): Limita o bloco de oferta aceito ao tamanho máximo definido 
para aquele bloco específico.
"""
def bid_block_limit_rule(model, g, b, t, s):
    return Umodel.rho_D_g_b_t_s[g, b, t, s] <= Umodel.P_D_g_b[g, b]
Umodel.bid_block_limit = pyo.Constraint(Umodel.G_D, b, Umodel.T, Umodel.S, rule=bid_block_limit_rule)


"""
Eq. (14): Define os limites inferiores e superiores para o deslocamento 
de carga em cada nó.
"""
def load_shifting_limit_rule(model, i, t, s):
    return (-Umodel.P_LS_i[i], Umodel.P_LS_i_t_s[i, t, s], Umodel.P_LS_i[i])
Umodel.load_shifting_limit = pyo.Constraint(Umodel.N_D, Umodel.T, Umodel.S, rule=load_shifting_limit_rule)


"""
Eq. (15): Garante que o deslocamento total de carga para cada nó ao 
longo de todos os períodos seja zero.
"""
def load_shifting_balance_rule(model, i, s):
    return sum(Umodel.P_LS_i_t_s[i, t, s] for t in Umodel.T) == 0
Umodel.load_shifting_balance = pyo.Constraint(Umodel.N_D, Umodel.S, rule=load_shifting_balance_rule)


"""
Eq. (16): Define os limites de carga e descarga para os sistemas 
de armazenamento de energia.
"""
def storage_power_limit_rule(model, i, t, s):
    return (-Umodel.P_ESS_i[i], Umodel.P_ESS_i_t_s[i, t, s], Umodel.P_ESS_i[i])
Umodel.storage_power_limit = pyo.Constraint(Umodel.N_D, Umodel.T, Umodel.S, rule=storage_power_limit_rule)


"""
Eq. (17): Define os limites mínimo e máximo do Estado de Carga (SOC) 
dos sistemas de armazenamento de energia para cada nó.
"""
def ess_soc_limit_rule(model, i, t, s):
    return (
        Umodel.SOC_ESS_MIN[i] <= Umodel.SOC_i_t_s[i, t, s] <= 
        Umodel.SOC_ESS_MAX[i]
    )
Umodel.ess_soc_limit = pyo.Constraint(Umodel.N_D, Umodel.T, Umodel.S, rule=ess_soc_limit_rule)


"""
Eq. (18-19): Calcula a evolução do Estado de Carga (SOC) 
dos sistemas de armazenamento ao longo do tempo.
"""
def ess_soc_update_rule(model, i, t, s):
    if t == Umodel.T.first():
        return Umodel.SOC_i_t_s[i, t, s] == Umodel.SOC_i_0[i]
    else:
        return (
            Umodel.SOC_i_t_s[i, t, s] == 
            Umodel.SOC_i_t_s[i, t-1, s] + 
            Umodel.PESS_i_t_s[i, t, s] * Umodel.eta_charge - 
            Umodel.PESS_i_t_s[i, t, s] / Umodel.eta_discharge
        )
Umodel.ess_soc_update = pyo.Constraint(Umodel.N_D, Umodel.T, Umodel.S, rule=ess_soc_update_rule)


"""
Eq. (20): Garante que o Estado de Carga (SOC) final 
seja igual ao valor inicial para cada sistema de armazenamento.
"""
def ess_final_soc_balance_rule(model, i, s):
    return (
        Umodel.SOC_i_t_s[i, Umodel.T.last(), s] == Umodel.SOC_i_0[i]
    )
Umodel.ess_final_soc_balance = pyo.Constraint(Umodel.N_D, Umodel.S, rule=ess_final_soc_balance_rule)


"""
Eq. (21): Calcula as emissões de carbono dos geradores 
conectados ao sistema de distribuição.
"""
def carbon_emission_distribution_rule(model, g, t, s):
    return (
        Umodel.nu_D_g_t_s[g, t, s] == 
        Umodel.gamma_D_g[g] * Umodel.P_D_g_t_s[g, t, s]
    )
Umodel.carbon_emission_distribution = pyo.Constraint(Umodel.G_D, Umodel.T, Umodel.S, rule=carbon_emission_distribution_rule)


"""
Eq. (22): Estabelece o balanço de carbono para o sistema de distribuição.
"""
def carbon_balance_distribution_rule(model, t, s):
    return (
        sum(Umodel.rho_D_g_t_s[g, t, s] for g in Umodel.G_D) == 
        -Umodel.rho_SE_t_s[t, s]
    )
Umodel.carbon_balance_distribution = pyo.Constraint(Umodel.T, Umodel.S, rule=carbon_balance_distribution_rule)


"""
Eq. (23): Define o limite de emissões de carbono para o sistema de distribuição.
"""
def carbon_emission_limit_distribution_rule(model, t, s):
    return (
        sum(Umodel.nu_D_g_t_s[g, t, s] + Umodel.rho_D_g_t_s[g, t, s] for g in Umodel.G_D) <= 
        Umodel.Gamma_D_t
    )
Umodel.carbon_emission_limit_distribution = pyo.Constraint(Umodel.T, Umodel.S, rule=carbon_emission_limit_distribution_rule)
