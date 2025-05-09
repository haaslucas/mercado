import pyomo.environ as pyo

# Criação do modelo Concreto
model = pyo.ConcreteModel()

#-----------------------------------------------------------------------
#  DECLARE SET OF BUS (Tradução de Parâmetros Escalares)
#-----------------------------------------------------------------------
model.VBASE = pyo.Param(initialize=24.90, doc='base voltage magnitude (kV)')
model.SBASE = pyo.Param(initialize=100, doc='base apparent power (MVA)')

# IBASE é calculado a partir de SBASE e VBASE.
# Em Pyomo, podemos definir isso com uma regra ou calcular após a inicialização.
# Por simplicidade inicial, vamos calcular diretamente e atribuir.
# model.IBASE = pyo.Param(initialize=(model.SBASE.value * 1e3 / model.VBASE.value))
# Uma forma mais robusta é usar uma rule se SBASE ou VBASE pudessem mudar dinamicamente,
# mas para parâmetros fixos, a atribuição direta após a criação do modelo ou na leitura de dados é comum.
# Vamos inicializá-lo como um parâmetro que será preenchido posteriormente ou calculado.
model.IBASE = pyo.Param(initialize= (100 * 1e3 / 24.90), doc='base current (A)')


model.SE_Capacity = pyo.Param(initialize=50, mutable=True) # mutable=True se o valor puder ser alterado como com 'let'

#-----------------------------------------------------------------------
#  DECLARE SETS (Tradução de Conjuntos)
#-----------------------------------------------------------------------
model.N = pyo.Set(doc='set of buses')
model.G_T = pyo.Set(doc='set of transmission-connected generators')
model.G_D = pyo.Set(doc='set of distribution-connected generators')
model.LS = pyo.Set(doc='set of load shifting nodes')
model.ESS = pyo.Set(doc='set of energy storage system nodes')
model.RES_T = pyo.Set(doc='set of transmission-connected renewable energy sources')
model.RES_D = pyo.Set(doc='set of distribution-connected renewable energy sources')

# set L within N cross N;
# L é um conjunto de tuplas (origem, destino) onde origem e destino pertencem a N.
# Isso geralmente é inicializado com dados.
model.L = pyo.Set(within=model.N * model.N, doc='set of branches')

model.T = pyo.Set(ordered=True, doc='set of time periods')
model.Trans_Lines = pyo.Set(doc='set of transmission lines')
model.Trans_Nodes = pyo.Set(doc='set of transmission nodes')
model.LDA = pyo.Set(doc='set of ??? (LDA - verificar significado no contexto)') # TODO: Confirmar o que LDA representa

# set WS within 1..3 cross Trans_Nodes cross (G_T union LDA);
# Este é um conjunto mais complexo.
# 1..3 pode ser representado por um RangeSet.
# A união G_T union LDA pode ser model.G_T | model.LDA
# WS será um conjunto de tuplas de 3 elementos.
model.Range1_3 = pyo.RangeSet(1,3)
model.WS = pyo.Set(within=model.Range1_3 * model.Trans_Nodes * (model.G_T | model.LDA),
                   doc='set of wholesale market participants') # A inicialização real virá dos dados

model.S = pyo.Set(doc='set of scenarios')

#-----------------------------------------------------------------------
# DECLARE ITEMS TO BE READ FROM THE .bus DATA FILE (Tradução de Parâmetros Indexados)
#-----------------------------------------------------------------------

model.NBUS = pyo.Param(model.N, within=pyo.Any, doc='bus name (symbolic)')
model.TBUS = pyo.Param(model.N, within=pyo.Any, doc='bus type')
model.Vmax = pyo.Param(model.N, doc='maximum voltage magnitude (kV)')
model.Vmin = pyo.Param(model.N, doc='minimum voltage magnitude (kV)')
model.Vnom = pyo.Param(model.N, doc='nominal voltage magnitude (kV)')
model.Ve = pyo.Param(model.N, doc='voltage magnitude (kV) - verificar uso, pode ser uma variável ou valor inicial')
model.Pgmx = pyo.Param(model.N, doc='maximum real power (kW)')
model.Pgmn = pyo.Param(model.N, doc='minimum real power (kW)')

model.R = pyo.Param(model.L, doc='branch resistance (ohm)')
model.X = pyo.Param(model.L, doc='branch reactance (ohm)')
model.Z = pyo.Param(model.L, mutable=True, doc='branch impedance (ohm)') # mutable=True pois é calculado no execute.run
model.Imax = pyo.Param(model.L, doc='branch maximum current (kA)')


model.load1 = pyo.Param(model.S, model.T, doc='load profile 1 factor')
model.load2 = pyo.Param(model.S, model.T, doc='load profile 2 factor')
model.wind1 = pyo.Param(model.S, model.T, doc='wind profile 1 factor')
model.wind2 = pyo.Param(model.S, model.T, doc='wind profile 2 factor')
model.solar = pyo.Param(model.S, model.T, doc='solar profile factor')
model.probability = pyo.Param(model.S, model.T, doc='probability of scenario s at time t')

# Parâmetros de perfil de carga (load shape parameters)
model.ind1 = pyo.Param(model.T, doc='industrial load profile 1 factor')
model.ind2 = pyo.Param(model.T, doc='industrial load profile 2 factor')
model.ind3 = pyo.Param(model.T, doc='industrial load profile 3 factor')

model.res1 = pyo.Param(model.T, doc='residential load profile 1 factor')
model.res4 = pyo.Param(model.T, doc='residential load profile 4 factor')
model.res5 = pyo.Param(model.T, doc='residential load profile 5 factor')

model.off1 = pyo.Param(model.T, doc='office load profile 1 factor')
model.off4 = pyo.Param(model.T, doc='office load profile 4 factor')

model.L2 = pyo.Param(model.T, doc='load profile L2 factor')
model.L3 = pyo.Param(model.T, doc='load profile L3 factor')
model.L4 = pyo.Param(model.T, doc='load profile L4 factor')
model.L6 = pyo.Param(model.T, doc='load profile L6 factor')
model.L8 = pyo.Param(model.T, doc='load profile L8 factor')

model.PCC = pyo.Param(model.LDA, doc='Point of Common Coupling for LDA entities')

# Parâmetros relacionados à transmissão
model.Trans_Load = pyo.Param(model.Trans_Nodes, model.T, model.S, doc='load at transmission nodes')
model.Trans_Shift_Max = pyo.Param(model.Trans_Nodes, doc='maximum load shifting at transmission nodes')
model.Trans_Shift_Min = pyo.Param(model.Trans_Nodes, doc='minimum load shifting at transmission nodes')
model.Trans_P_ESS_min = pyo.Param(model.Trans_Nodes, doc='minimum power from ESS at transmission nodes')
model.Trans_P_ESS_max = pyo.Param(model.Trans_Nodes, doc='maximum power from ESS at transmission nodes')
model.Trans_SOC_max = pyo.Param(model.Trans_Nodes, doc='maximum State Of Charge for ESS at transmission nodes')
model.Trans_SOC_min = pyo.Param(model.Trans_Nodes, doc='minimum State Of Charge for ESS at transmission nodes')
model.Trans_SOC_initial = pyo.Param(model.Trans_Nodes, doc='initial State Of Charge for ESS at transmission nodes')
model.Trans_Type = pyo.Param(model.Trans_Nodes, within=pyo.Any, doc='type of transmission node (e.g., load type)')
model.Trans_Inst_Cap = pyo.Param(model.Trans_Nodes, doc='installed capacity at transmission nodes')

# Parâmetros de fontes de energia renovável conectadas à transmissão (RES-T)
model.Trans_Wind_Inj = pyo.Param(model.RES_T, model.T, model.S, doc='wind injection at transmission RES')
model.Trans_Solar_Inj = pyo.Param(model.RES_T, model.T, model.S, doc='solar injection at transmission RES')
model.RES_type_t = pyo.Param(model.RES_T, within=pyo.Any, doc='type of transmission RES (e.g., SOLAR, WIND)')
model.RES_Node_t = pyo.Param(model.RES_T, doc='node of transmission RES') # Assumindo que mapeia para um nó em Trans_Nodes
model.RES_Cap_t = pyo.Param(model.RES_T, doc='capacity of transmission RES')
model.Wind_min_t = pyo.Param(model.RES_T, doc='minimum wind power output for transmission RES (pu or MW)') # Verificar unidade
model.Wind_max_t = pyo.Param(model.RES_T, doc='maximum wind power output for transmission RES (pu or MW)') # Verificar unidade
model.Wind_nom_t = pyo.Param(model.RES_T, doc='nominal wind power output for transmission RES (pu or MW)') # Verificar unidade

# Parâmetros do sistema de distribuição (Cargas, Load Shifting, ESS)
model.Dist_Load = pyo.Param(model.N, model.T, model.S, doc='load at distribution nodes')
model.LS_Node = pyo.Param(model.LS, doc='mapping of load shifting entity to distribution node N') # No .dat, LS é 1..8, LS_Node[ls] é um nó em N
model.Dist_Shift_Max = pyo.Param(model.LS, doc='maximum load shifting at distribution nodes')
model.Dist_Shift_Min = pyo.Param(model.LS, doc='minimum load shifting at distribution nodes')
model.ESS_Node = pyo.Param(model.ESS, doc='mapping of ESS entity to distribution node N') # No .dat, ESS é 1..2, ESS_Node[ess] é um nó em N
model.Dist_P_ESS_min = pyo.Param(model.ESS, doc='minimum power from ESS at distribution nodes')
model.Dist_P_ESS_max = pyo.Param(model.ESS, doc='maximum power from ESS at distribution nodes')
model.Dist_SOC_max = pyo.Param(model.ESS, doc='maximum State Of Charge for ESS at distribution nodes')
model.Dist_SOC_min = pyo.Param(model.ESS, doc='minimum State Of Charge for ESS at distribution nodes')
model.Dist_SOC_initial = pyo.Param(model.ESS, doc='initial State Of Charge for ESS at distribution nodes')
model.Dist_Type = pyo.Param(model.N, within=pyo.Any, doc='type of distribution node (e.g., load type like IND1, RES1)')
model.Dist_Inst_Cap = pyo.Param(model.N, doc='installed capacity at distribution nodes')

# Parâmetros de fontes de energia renovável conectadas à distribuição (RES-D)
model.Dist_Wind_Inj = pyo.Param(model.RES_D, model.T, model.S, doc='wind injection at distribution RES')
model.Dist_Solar_Inj = pyo.Param(model.RES_D, model.T, model.S, doc='solar injection at distribution RES')
model.RES_type_d = pyo.Param(model.RES_D, within=pyo.Any, doc='type of distribution RES (e.g., SOLAR, WIND)')
model.RES_Node_d = pyo.Param(model.RES_D, doc='node of distribution RES') # Assumindo que mapeia para um nó em N
model.RES_Cap_d = pyo.Param(model.RES_D, doc='capacity of distribution RES')
model.Wind_min_d = pyo.Param(model.RES_D, doc='minimum wind power output for distribution RES (pu or MW)') # Verificar unidade
model.Wind_max_d = pyo.Param(model.RES_D, doc='maximum wind power output for distribution RES (pu or MW)') # Verificar unidade
model.Wind_nom_d = pyo.Param(model.RES_D, doc='nominal wind power output for distribution RES (pu or MW)') # Verificar unidade

# Parâmetros das linhas de transmissão
model.Trans_Incidencia = pyo.Param(model.Trans_Nodes, model.Trans_Lines, mutable=True, doc='incidence matrix for transmission lines') # Mutable pois é calculada no execute.run
model.Trans_Reactance = pyo.Param(model.Trans_Lines, doc='reactance of transmission lines')
model.Trans_Capacity = pyo.Param(model.Trans_Lines, mutable=True, doc='capacity of transmission lines') # Mutable pois é ajustada no execute.run
model.Trans_Status = pyo.Param(model.Trans_Lines, doc='status of transmission lines (on/off)')
model.From = pyo.Param(model.Trans_Lines, doc='origin node of transmission lines') # Usado para construir Trans_Incidencia
model.To = pyo.Param(model.Trans_Lines, doc='destination node of transmission lines') # Usado para construir Trans_Incidencia

# Parâmetros das linhas de distribuição
# No input.dat, R e X são fornecidos para o conjunto L (distribuição).
# Z é calculado no execute.run.
# Dist_X e Dist_R não são explicitamente declarados no Novo.mod, mas R{L} e X{L} são.
# Vou assumir que R{L} e X{L} são os parâmetros para resistência e reatância da distribuição.
# Dist_Z e Dist_Status são declarados no Novo.mod.
# model.Dist_X = pyo.Param(model.L, doc='reactance of distribution lines') # Já coberto por model.X
# model.Dist_R = pyo.Param(model.L, doc='resistance of distribution lines') # Já coberto por model.R
model.Dist_Z = pyo.Param(model.L, mutable=True, doc='impedance of distribution lines') # Calculado no execute.run, igual a model.Z
model.Dist_Status = pyo.Param(model.L, doc='status of distribution lines (on/off)')

# Parâmetros dos geradores conectados à transmissão (G_T)
model.Pmax_t = pyo.Param(model.G_T, mutable=True, doc='maximum active power for transmission generator')
model.Pmin_t = pyo.Param(model.G_T, mutable=True, doc='minimum active power for transmission generator')
model.a_t = pyo.Param(model.G_T, mutable=True, doc='cost coefficient a for transmission generator')
model.b_t = pyo.Param(model.G_T, mutable=True, doc='cost coefficient b for transmission generator')
model.c_t = pyo.Param(model.G_T, doc='cost coefficient c for transmission generator') # Não modificado no execute.run
model.G_Node = pyo.Param(model.G_T, doc='node of transmission generator')
model.G_Owner_t = pyo.Param(model.G_T, doc='owner of transmission generator')
model.Carbon_Cost_t = pyo.Param(model.G_T, mutable=True, doc='carbon cost for transmission generator')

# Parâmetros dos geradores conectados à distribuição (G_D)
model.Pmax_d = pyo.Param(model.G_D, mutable=True, doc='maximum active power for distribution generator')
model.Pmin_d = pyo.Param(model.G_D, mutable=True, doc='minimum active power for distribution generator')
model.a_d = pyo.Param(model.G_D, mutable=True, doc='cost coefficient a for distribution generator')
model.b_d = pyo.Param(model.G_D, mutable=True, doc='cost coefficient b for distribution generator')
model.c_d = pyo.Param(model.G_D, doc='cost coefficient c for distribution generator') # Não modificado no execute.run
model.G_Owner_d = pyo.Param(model.G_D, doc='owner of distribution generator')
model.G_LDA = pyo.Param(model.G_D, doc='LDA mapping for distribution generator') # Verificar o significado exato de G_LDA
model.G_LDA_Node = pyo.Param(model.G_D, doc='node associated with G_LDA for distribution generator')
model.Carbon_Cost_d = pyo.Param(model.G_D, mutable=True, doc='carbon cost for distribution generator')

# Parâmetros de carbono
model.Carbon_Limit_Trans = pyo.Param(model.T, mutable=True, doc='carbon limit for transmission system')
model.Carbon_Limit_Dist = pyo.Param(model.T, mutable=True, doc='carbon limit for distribution system')
model.Emission_Weighted_Average = pyo.Param(mutable=True, doc='emission weighted average cost/factor')

#-----------------------------------------------------------------------
# DECLARE VARIABLES (Tradução de Variáveis)
#-----------------------------------------------------------------------

# Variáveis Primárias
model.VM = pyo.Var(model.N, model.T, model.S, doc='voltage magnitude at bus n, time t, scenario s')
model.I = pyo.Var(model.L, model.T, model.S, doc='current in branch l, time t, scenario s')

model.Trans_Shift = pyo.Var(model.Trans_Nodes, model.T, model.S, doc='load shifting at transmission node, time t, scenario s')
model.Trans_P_ESS = pyo.Var(model.Trans_Nodes, model.T, model.S, doc='power from ESS at transmission node, time t, scenario s')
model.Trans_Flow = pyo.Var(model.Trans_Lines, model.T, model.S, doc='power flow in transmission line, time t, scenario s')
model.Trans_Theta = pyo.Var(model.Trans_Nodes, model.T, model.S, doc='voltage angle at transmission node, time t, scenario s')

model.Dist_Shift = pyo.Var(model.LS, model.T, model.S, doc='load shifting at distribution load, time t, scenario s')
model.Dist_P_ESS = pyo.Var(model.ESS, model.T, model.S, doc='power from ESS at distribution node, time t, scenario s')
model.Dist_SOC = pyo.Var(model.ESS, model.T, model.S, doc='state of charge of ESS at distribution node, time t, scenario s')
# Dist_Flow não é usado nas restrições do .mod, mas P, Q, Dist_Pfm, Dist_Pto são.
# model.Dist_Flow = pyo.Var(model.L, model.T, model.S, doc='flow in distribution line l, time t, scenario s')
model.Dist_Pfm = pyo.Var(model.L, model.T, model.S, doc='active power flow "from" in distribution line l, time t, scenario s')
model.Dist_Pto = pyo.Var(model.L, model.T, model.S, doc='active power flow "to" in distribution line l, time t, scenario s')
model.Dist_Qfm = pyo.Var(model.L, model.T, model.S, doc='reactive power flow "from" in distribution line l, time t, scenario s')
model.Dist_Qto = pyo.Var(model.L, model.T, model.S, doc='reactive power flow "to" in distribution line l, time t, scenario s')
model.P = pyo.Var(model.L, model.T, model.S, doc='active power flow in distribution line l, time t, scenario s')
model.Q = pyo.Var(model.L, model.T, model.S, doc='reactive power flow in distribution line l, time t, scenario s')

model.P_thermal_dist = pyo.Var(model.G_D, model.T, model.S, domain=pyo.NonNegativeReals, doc='active power from distribution generator, time t, scenario s')
model.Q_thermal_dist = pyo.Var(model.G_D, model.T, model.S, doc='reactive power from distribution generator, time t, scenario s')
model.P_thermal_trans = pyo.Var(model.G_T, model.T, model.S, domain=pyo.NonNegativeReals, doc='active power from transmission generator, time t, scenario s')
model.P_DSO = pyo.Var(model.T, model.S, doc='power exchanged by DSO with transmission, time t, scenario s')
model.Q_DSO = pyo.Var(model.T, model.S, doc='reactive power exchanged by DSO with transmission, time t, scenario s')

model.Footprint_trans = pyo.Var(model.G_T, model.T, model.S, doc='carbon footprint of transmission generator, time t, scenario s')
model.Footprint_dist = pyo.Var(model.G_D, model.T, model.S, doc='carbon footprint of distribution generator, time t, scenario s')
model.Carbon_T = pyo.Var(model.G_T, model.T, model.S, doc='carbon traded by transmission generator, time t, scenario s') # Ou emissão de carbono
model.Carbon_SE = pyo.Var(model.T, model.S, doc='carbon exchanged at substation/SE, time t, scenario s')

# Outras Variáveis (algumas com limites)
# A variável Bid é usada na função objetivo e em Deriv_PDSO.
# A variável Carbon_Price é usada na função objetivo e em Deriv_Carbon_SE.
# No contexto de um problema de otimização padrão, se Bid e Carbon_Price são decisões do modelo,
# elas devem ser variáveis. Se são preços dados (parâmetros), devem ser Params.
# O arquivo .mod as declara como 'var' com limites, sugerindo que são variáveis de decisão.
# No entanto, as restrições "Deriv_PDSO" e "Deriv_Carbon_SE" parecem ser condições de KKT
# de um problema de nível inferior ou de um mercado, onde Bid e Carbon_Price poderiam ser
# multiplicadores de Lagrange (variáveis duplas) desse problema inferior.
# Para uma tradução inicial do problema do DSO, vou assumi-las como variáveis de decisão do DSO,
# conforme declarado no .mod. Se este for um modelo bi-nível, a interpretação mudaria.

def bid_bounds_rule(model, t, s):
    return (None, 500.0) # (lower_bound, upper_bound)
model.Bid = pyo.Var(model.T, model.S, bounds=bid_bounds_rule, doc='Bid price by DSO, time t, scenario s')

def carbon_price_bounds_rule(model, t, s):
    return (None, 5.0)
model.Carbon_Price = pyo.Var(model.T, model.S, bounds=carbon_price_bounds_rule, doc='Carbon price, time t, scenario s')

model.DSO_Revenue = pyo.Var(model.T, model.S, doc='DSO revenue, time t, scenario s') # Usada em restrição comentada
model.GenCosts_dist = pyo.Var(model.G_D, model.T, model.S, domain=pyo.NonNegativeReals, doc='generation costs for distribution generators')
model.GenCosts_trans = pyo.Var(model.G_T, model.T, model.S, domain=pyo.NonNegativeReals, doc='generation costs for transmission generators')

#-----------------------------------------------------------------------
# DUAL VARIABLES (Tradução das Variáveis Duplas para KKT)
#-----------------------------------------------------------------------
# Estas são as variáveis duplas associadas às restrições do problema de nível inferior (ISO/Transmissão)
model.lambda_kkt = pyo.Var(model.Trans_Nodes, model.T, model.S, doc='Dual variable for Trans_Power_Balance') # Renomeado para evitar conflito com a função lambda do Python
model.omega_U = pyo.Var(model.Trans_Lines, model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for Trans_Flow_Upper')
model.omega_L = pyo.Var(model.Trans_Lines, model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for Trans_Flow_Lower')

model.kappa_U = pyo.Var(model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for MAXIMUM_SE_LIMIT')
model.kappa_L = pyo.Var(model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for MINIMUM_SE_LIMIT')
model.alpha_U = pyo.Var(model.G_T, model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for Trans_Gen_Upper')
model.alpha_L = pyo.Var(model.G_T, model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for Trans_Gen_Lower')

# Duals para restrições de rampa. Definidas como NonNegativeReals com base nas constraints Omega_Lower/Upper etc.
model.ramp_up = pyo.Var(model.G_T, model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for Trans_Gen_Upper_Ramp (Ramp Down Limit)')
model.ramp_low = pyo.Var(model.G_T, model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for Trans_Gen_Lower_Ramp (Ramp Up Limit)')

model.mu_kkt = pyo.Var(model.Trans_Lines, model.T, model.S, doc='Dual variable for Calculate_Trans_PowerFlow') # Renomeado para evitar conflito com outros usos de mu
model.phi_kkt = pyo.Var(model.G_T, model.T, model.S, doc='Dual variable for Calculate_Footprint_Trans') # Renomeado
model.psi_kkt = pyo.Var(model.T, model.S, doc='Dual variable for Carbon_Balance_Trans') # Renomeado
model.psi_L_kkt = pyo.Var(model.T, model.S, doc='Dual variable (psi_L), não utilizada nas restrições KKT do Novo.mod') # Renomeado

# Epsilon é usado em restrições comentadas no AMPL, declarado como free aqui.
model.epsilon_kkt = pyo.Var(model.G_T, model.T, model.S, doc='Dual variable (epsilon), usada em KKTs comentadas') # Renomeado

model.Xi_kkt = pyo.Var(model.T, model.S, domain=pyo.NonNegativeReals, doc='Dual variable for Carbon_Limits_Trans') # Renomeado

# Eta_U e Eta_L são usados em restrições comentadas no AMPL, declarados como free aqui.
model.eta_U_kkt = pyo.Var(model.T, model.S, doc='Dual variable (eta_U), usada em KKTs comentadas') # Renomeado
model.eta_L_kkt = pyo.Var(model.T, model.S, doc='Dual variable (eta_L), usada em KKTs comentadas') # Renomeado


#-----------------------------------------------------------------------
# OBJECTIVE FUNCTION (Tradução da Função Objetivo)
#-----------------------------------------------------------------------
def dso_costs_rule(model):
    cost = 0
    for s in model.S:
        for t in model.T:
            # Custo da energia comprada/vendida do DSO para o sistema de transmissão
            energy_cost_dso_trans = model.Bid[t,s] * model.P_DSO[t,s]
            
            # Custos de geração dos geradores de distribuição
            dist_gen_costs = sum(model.GenCosts_dist[g,t,s] for g in model.G_D)
            
            # Custos/Receitas de carbono na subestação
            carbon_exchange_cost = model.Carbon_Price[t,s] * model.Carbon_SE[t,s]
            
            cost += model.probability[s,t] * (energy_cost_dso_trans + dist_gen_costs + carbon_exchange_cost)
    return cost

model.DSO_Costs = pyo.Objective(rule=dso_costs_rule, sense=pyo.minimize)

#-----------------------------------------------------------------------
# CONSTRAINTS (Tradução das Restrições)
#-----------------------------------------------------------------------

# s.t. Calculate_GenCosts_Dist {t in T, g in G_D, s in S}:
# 	GenCosts_dist[g,t,s] = a_d[g]*P_thermal_dist[g,t,s]^2 + b_d[g]*P_thermal_dist[g,t,s] + c_d[g];
def calculate_gen_costs_dist_rule(model, t, g, s):
    return model.GenCosts_dist[g,t,s] == (model.a_d[g] * model.P_thermal_dist[g,t,s]**2 +
                                           model.b_d[g] * model.P_thermal_dist[g,t,s] +
                                           model.c_d[g])
model.Calculate_GenCosts_Dist = pyo.Constraint(model.T, model.G_D, model.S, rule=calculate_gen_costs_dist_rule)

# s.t. Calculate_GenCosts_Trans {t in T, g in G_T, s in S}:
# 	GenCosts_trans[g,t,s] = a_t[g]*P_thermal_trans[g,t,s]^2 + b_t[g]*P_thermal_trans[g,t,s] + c_t[g];
def calculate_gen_costs_trans_rule(model, t, g, s):
    return model.GenCosts_trans[g,t,s] == (model.a_t[g] * model.P_thermal_trans[g,t,s]**2 +
                                            model.b_t[g] * model.P_thermal_trans[g,t,s] +
                                            model.c_t[g])
model.Calculate_GenCosts_Trans = pyo.Constraint(model.T, model.G_T, model.S, rule=calculate_gen_costs_trans_rule)

# s.t. ACTIVE_POWER_BALANCE{n_node in N, t in T, s in S}:
# 	sum{g in G_D:G_LDA_Node[g] == n_node}(P_thermal_dist[g,t,s])
# 	-sum{(n_node,m) in L}(P[n_node,m,t,s] + R[n_node,m]*I[n_node,m,t,s])
# 	+sum{(l,n_node) in L}(P[l,n_node,t,s])
# 	+(if 1 == n_node then P_DSO[t,s])
# 	-sum{i in LS:LS_Node[i] == n_node} (Dist_Shift[i,t,s])
# 	+sum{i in ESS:ESS_Node[i] == n_node} (Dist_P_ESS[i,t,s])
# 	= Dist_Load[n_node,t,s];
def active_power_balance_rule(model, n_node, t, s):
    generation_at_node = sum(model.P_thermal_dist[g,t,s] for g in model.G_D if model.G_LDA_Node[g] == n_node)
    
    power_exchange_dso = 0
    # Assuming node 1 is the point of common coupling with transmission.
    # The set N is 1-indexed based on input.dat.
    if n_node == 1: # Or use a parameter for PCC node if it's not always 1
        power_exchange_dso = model.P_DSO[t,s]
        
    load_shifting_at_node = sum(model.Dist_Shift[i,t,s] for i in model.LS if model.LS_Node[i] == n_node)
    ess_power_at_node = sum(model.Dist_P_ESS[i,t,s] for i in model.ESS if model.ESS_Node[i] == n_node)

    # Sum of power flowing out of n_node
    # P[line,t,s] is power on line 'line' at time 't', scenario 's'.
    # R[line]*I[line,t,s] - this term is dimensionally unusual (V not W). Translating literally.
    power_flow_out = sum(model.P[line,t,s] + model.R[line]*model.I[line,t,s] 
                         for line in model.L if line[0] == n_node)
    
    # Sum of power flowing into n_node
    power_flow_in = sum(model.P[line,t,s] 
                        for line in model.L if line[1] == n_node)

    return (generation_at_node + 
            power_exchange_dso + 
            ess_power_at_node +
            power_flow_in - # Flows into n_node are positive contributions
            power_flow_out - # Flows out of n_node are negative contributions
            load_shifting_at_node # Load shift can be positive (load reduction) or negative (load increase)
            == model.Dist_Load[n_node,t,s])

model.ACTIVE_POWER_BALANCE = pyo.Constraint(model.N, model.T, model.S, rule=active_power_balance_rule)

# s.t. LOSS_REAL_POWER {(i,j) in L, t in T, s in S}: 
# 	(Dist_Pfm[i,j,t,s] + Dist_Pto[i,j,t,s]) = R[i,j] *I[i,j,t,s]^2;
def loss_real_power_rule(model, i_node, j_node, t, s): # Iterating over L directly is better
    line = (i_node, j_node)
    if line not in model.L: # Should not happen if rule is indexed by model.L
        return pyo.Constraint.Skip
    return model.Dist_Pfm[line,t,s] + model.Dist_Pto[line,t,s] == model.R[line] * model.I[line,t,s]**2
# Correct indexing for the constraint definition:
model.LOSS_REAL_POWER = pyo.Constraint(model.L, model.T, model.S, rule=lambda model, i,j,t,s: \
    model.Dist_Pfm[(i,j),t,s] + model.Dist_Pto[(i,j),t,s] == model.R[(i,j)] * model.I[(i,j),t,s]**2)


# s.t. REAL_POWER_FLOW {(i,j) in L, t in T, s in S}:
# 	(Dist_Pfm[i,j,t,s] - Dist_Pto[i,j,t,s]) = R[i,j]/Z[i,j]^2 * (VM[i,t,s]^2 - VM[j,t,s]^2);
def real_power_flow_rule(model, i_node, j_node, t, s):
    line = (i_node, j_node)
    if line not in model.L:
        return pyo.Constraint.Skip
    
    # Add epsilon for numerical stability if Z might be zero
    # If Z is zero, R must also be zero, so R / (Z^2+eps) = 0 / eps = 0.
    denominator = model.Z[line]**2 + 1e-9 
    
    return (model.Dist_Pfm[line,t,s] - model.Dist_Pto[line,t,s] == 
            model.R[line] / denominator * 
            (model.VM[i_node,t,s]**2 - model.VM[j_node,t,s]**2))
# Correct indexing for the constraint definition:
model.REAL_POWER_FLOW = pyo.Constraint(model.L, model.T, model.S, rule=lambda model, i,j,t,s: \
    model.Dist_Pfm[(i,j),t,s] - model.Dist_Pto[(i,j),t,s] == \
    model.R[(i,j)] / (model.Z[(i,j)]**2 + 1e-9) * \
    (model.VM[i,t,s]**2 - model.VM[j,t,s]**2))


# s.t. CURRENT_FLOW {(i,j) in L, t in T, s in S}:
# 	I[i,j,t,s] = (VM[i,t,s]-VM[j,t,s])/Z[i,j];
def current_flow_rule(model, i_node, j_node, t, s):
    line = (i_node, j_node)
    if line not in model.L:
        return pyo.Constraint.Skip
    # Add epsilon for numerical stability if Z might be zero
    denominator = model.Z[line] + 1e-9
    return model.I[line,t,s] == (model.VM[i_node,t,s] - model.VM[j_node,t,s]) / denominator
# Correct indexing for the constraint definition:
model.CURRENT_FLOW = pyo.Constraint(model.L, model.T, model.S, rule=lambda model, i,j,t,s: \
    model.I[(i,j),t,s] == (model.VM[i,t,s] - model.VM[j,t,s]) / (model.Z[(i,j)] + 1e-9))

#---- LINE RATING LIMITS
# s.t. MINIMUM_CURRENT_FLOW {(n,m) in L, t in T, s in S}: 
# 	I[n,m,t,s] >= -Imax[n,m];
def minimum_current_flow_rule(model, n, m, t, s):
    line = (n,m)
    return model.I[line,t,s] >= -model.Imax[line]
model.MINIMUM_CURRENT_FLOW = pyo.Constraint(model.L, model.T, model.S, rule=minimum_current_flow_rule)
  
# s.t. MAXIMUM_CURRENT_FLOW {(n,m) in L, t in T, s in S}: 
#     I[n,m,t,s] <= Imax[n,m];
def maximum_current_flow_rule(model, n, m, t, s):
    line = (n,m)
    return model.I[line,t,s] <= model.Imax[line]
model.MAXIMUM_CURRENT_FLOW = pyo.Constraint(model.L, model.T, model.S, rule=maximum_current_flow_rule)											

#---- VOLTAGE MAGNITUDE LIMITS
# s.t. MAXIMUM_VOLTAGE_MAGNITUDE {t in T,n in N, s in S}:
#     VM[n,t,s] <= Vmax[n];
def maximum_voltage_magnitude_rule(model, t, n, s):
    return model.VM[n,t,s] <= model.Vmax[n]
model.MAXIMUM_VOLTAGE_MAGNITUDE = pyo.Constraint(model.T, model.N, model.S, rule=maximum_voltage_magnitude_rule)												
  
# s.t. MINIMUM_VOLTAGE_MAGNITUDE {t in T,n in N, s in S}:
#   	VM[n,t,s] >= Vmin[n];
def minimum_voltage_magnitude_rule(model, t, n, s):
    return model.VM[n,t,s] >= model.Vmin[n]
model.MINIMUM_VOLTAGE_MAGNITUDE = pyo.Constraint(model.T, model.N, model.S, rule=minimum_voltage_magnitude_rule)												
  
#---- ACTIVE GENERATION LIMITS (Distribution Generators)
# s.t. MAXIMUM_DG_ACTIVE_POWER {t in T, g in G_D, s in S}:
#     P_thermal_dist[g,t,s] <=   Pmax_d[g];
def maximum_dg_active_power_rule(model, t, g, s):
    return model.P_thermal_dist[g,t,s] <= model.Pmax_d[g]
model.MAXIMUM_DG_ACTIVE_POWER = pyo.Constraint(model.T, model.G_D, model.S, rule=maximum_dg_active_power_rule)

# s.t. MINIMUM_DG_ACTIVE_POWER {t in T, g in G_D, s in S}:
#     P_thermal_dist[g,t,s] >=   Pmin_d[g];
def minimum_dg_active_power_rule(model, t, g, s):
    return model.P_thermal_dist[g,t,s] >= model.Pmin_d[g]
model.MINIMUM_DG_ACTIVE_POWER = pyo.Constraint(model.T, model.G_D, model.S, rule=minimum_dg_active_power_rule)

#---- LOAD SHIFTING LIMITS
# s.t. MAXIMUM_LS_LIMIT {t in T, n_ls in LS, s in S}:
# 	Dist_Shift[n_ls,t,s] <= Dist_Shift_Max[n_ls] *Dist_Load[LS_Node[n_ls],t,s];
def maximum_ls_limit_rule(model, t, n_ls, s):
    # LS_Node[n_ls] gives the actual node in N for this load shifting entity
    actual_node = model.LS_Node[n_ls]
    # It's possible Dist_Load might not be defined for all (actual_node, t, s) if load is zero.
    # However, AMPL would treat missing param values as an error during constraint generation
    # unless 'default' is used. Pyomo would error if model.Dist_Load[actual_node,t,s] is not populated.
    # Assuming Dist_Load is fully populated or 0 for relevant nodes.
    return model.Dist_Shift[n_ls,t,s] <= model.Dist_Shift_Max[n_ls] * model.Dist_Load[actual_node,t,s]
model.MAXIMUM_LS_LIMIT = pyo.Constraint(model.T, model.LS, model.S, rule=maximum_ls_limit_rule)

# s.t. MINIMUM_LS_LIMIT {t in T, n_ls in LS, s in S}:
# 	-Dist_Shift[n_ls,t,s] <= -Dist_Shift_Min[n_ls] *Dist_Load[LS_Node[n_ls],t,s];
# This is equivalent to: Dist_Shift[n_ls,t,s] >= Dist_Shift_Min[n_ls] * Dist_Load[LS_Node[n_ls],t,s]
def minimum_ls_limit_rule(model, t, n_ls, s):
    actual_node = model.LS_Node[n_ls]
    return model.Dist_Shift[n_ls,t,s] >= model.Dist_Shift_Min[n_ls] * model.Dist_Load[actual_node,t,s]
model.MINIMUM_LS_LIMIT = pyo.Constraint(model.T, model.LS, model.S, rule=minimum_ls_limit_rule)

# s.t. ZERO_SUM_LS {n_ls in LS, s in S}:
# 	sum{t_period in T} Dist_Shift[n_ls,t_period,s] = 0;
def zero_sum_ls_rule(model, n_ls, s):
    return sum(model.Dist_Shift[n_ls,t_period,s] for t_period in model.T) == 0
model.ZERO_SUM_LS = pyo.Constraint(model.LS, model.S, rule=zero_sum_ls_rule)

#---- ENERGY STORAGE SYSTEMS LIMITS (Distribution)
# s.t. Dist_ESS_Limits  {n in ESS, t in T, s in S}: Dist_P_ESS[n,t,s] >= Dist_P_ESS_min[n];
def dist_ess_limits_rule(model, n_ess, t, s):
    return model.Dist_P_ESS[n_ess,t,s] >= model.Dist_P_ESS_min[n_ess]
model.Dist_ESS_Limits = pyo.Constraint(model.ESS, model.T, model.S, rule=dist_ess_limits_rule)

# s.t. Dist_ESS_Limits2 {n in ESS, t in T, s in S}: Dist_P_ESS[n,t,s] <= Dist_P_ESS_max[n];
def dist_ess_limits2_rule(model, n_ess, t, s):
    return model.Dist_P_ESS[n_ess,t,s] <= model.Dist_P_ESS_max[n_ess]
model.Dist_ESS_Limits2 = pyo.Constraint(model.ESS, model.T, model.S, rule=dist_ess_limits2_rule)

# s.t. Dispatch_Dist_SOCLimits  {n in ESS, t in T, s in S}: Dist_SOC[n,t,s] >= Dist_SOC_min[n];
def dispatch_dist_soc_limits_rule(model, n_ess, t, s):
    return model.Dist_SOC[n_ess,t,s] >= model.Dist_SOC_min[n_ess]
model.Dispatch_Dist_SOCLimits = pyo.Constraint(model.ESS, model.T, model.S, rule=dispatch_dist_soc_limits_rule)

# s.t. Dispatch_Dist_SOCLimits2 {n in ESS, t in T, s in S}: Dist_SOC[n,t,s] <= Dist_SOC_max[n];
def dispatch_dist_soc_limits2_rule(model, n_ess, t, s):
    return model.Dist_SOC[n_ess,t,s] <= model.Dist_SOC_max[n_ess]
model.Dispatch_Dist_SOCLimits2 = pyo.Constraint(model.ESS, model.T, model.S, rule=dispatch_dist_soc_limits2_rule)

# s.t. Dispatch_Dist_CalculateSOC {n in ESS, t in T, s in S}:
# (if t = first(T) then Dist_SOC_initial[n] else Dist_SOC[n,prev(t),s]) - Dist_P_ESS[n,t,s] = Dist_SOC[n,t,s];
def dispatch_dist_calculate_soc_rule(model, n_ess, t, s):
    if t == model.T.first(): # model.T is 1-indexed and ordered
        soc_previous = model.Dist_SOC_initial[n_ess]
    else:
        soc_previous = model.Dist_SOC[n_ess, model.T.prev(t), s]
    # Assuming Dist_P_ESS is positive for discharge and negative for charge.
    # The equation implies: SOC_current = SOC_previous - Power_Discharged (or + Power_Charged if P_ESS is negative for charge)
    # If Dist_P_ESS > 0 means discharging, then SOC decreases.
    # If Dist_P_ESS < 0 means charging, then SOC increases.
    # The formulation SOC_prev - P_ESS = SOC_current is consistent with P_ESS > 0 for discharge.
    return soc_previous - model.Dist_P_ESS[n_ess,t,s] == model.Dist_SOC[n_ess,t,s]
model.Dispatch_Dist_CalculateSOC = pyo.Constraint(model.ESS, model.T, model.S, rule=dispatch_dist_calculate_soc_rule)

# s.t. Dispatch_Dist_SameSOC {n in ESS, s in S}: Dist_SOC_initial[n] = Dist_SOC[n, last(T),s];
def dispatch_dist_same_soc_rule(model, n_ess, s):
    return model.Dist_SOC_initial[n_ess] == model.Dist_SOC[n_ess, model.T.last(), s]
model.Dispatch_Dist_SameSOC = pyo.Constraint(model.ESS, model.S, rule=dispatch_dist_same_soc_rule)

#---- CARBON RELATED CONSTRAINTS (Distribution)
# s.t. Calculate_Footprint_Dist {t in T, g in G_D, s in S}:
#	Footprint_dist[g,t,s] = Carbon_Cost_d[g] * P_thermal_dist[g,t,s] ;
def calculate_footprint_dist_rule(model, t, g, s):
    return model.Footprint_dist[g,t,s] == model.Carbon_Cost_d[g] * model.P_thermal_dist[g,t,s]
model.Calculate_Footprint_Dist = pyo.Constraint(model.T, model.G_D, model.S, rule=calculate_footprint_dist_rule)
	
# s.t. Dist_Carbon_Trade_Limit2 {t in T, s in S}:
#	sum{g in G_D} (Footprint_dist[g,t,s]) - Carbon_SE[t,s] <= Carbon_Limit_Dist[t] ;	
def dist_carbon_trade_limit2_rule(model, t, s):
    sum_footprint_dist = sum(model.Footprint_dist[g,t,s] for g in model.G_D)
    return sum_footprint_dist - model.Carbon_SE[t,s] <= model.Carbon_Limit_Dist[t]
model.Dist_Carbon_Trade_Limit2 = pyo.Constraint(model.T, model.S, rule=dist_carbon_trade_limit2_rule)

#-----------------------------------------------------------------------
# TRANSMISSION PRIMAL CONSTRAINTS
#-----------------------------------------------------------------------

# s.t. Calculate_Trans_PowerFlow {l in Trans_Lines, t in T, s in S}: 
# 	Trans_Flow[l,t,s] - Trans_Status[l] * 1/Trans_Reactance[l] * 
# 	sum{n in Trans_Nodes}(Trans_Incidencia[n,l]*Trans_Theta[n,t,s]) = 0;
def calculate_trans_power_flow_rule(model, l, t, s):
    if model.Trans_Reactance[l] == 0: # Avoid division by zero if reactance is zero
        # If reactance is zero, for a non-zero flow, angle difference must be zero.
        # Or, if status is 0, flow must be 0.
        # This case might need specific handling based on model assumptions for zero reactance lines.
        # Assuming active lines have non-zero reactance. If status is 0, flow should be 0.
        if model.Trans_Status[l] == 0:
            return model.Trans_Flow[l,t,s] == 0
        else: # Status is 1, Reactance is 0. This implies theta_from == theta_to for finite flow.
              # The sum term should be zero if flow is not infinite.
              # For DC flow, Trans_Flow = (Theta_from - Theta_to) / X. If X=0, Theta_from must equal Theta_to unless flow is infinite.
              # A more robust way: X * Trans_Flow = Theta_from - Theta_to
            theta_difference = sum(model.Trans_Incidencia[n,l] * model.Trans_Theta[n,t,s] for n in model.Trans_Nodes)
            return model.Trans_Reactance[l] * model.Trans_Flow[l,t,s] == model.Trans_Status[l] * theta_difference
    
    # Original formulation:
    angle_difference_term = sum(model.Trans_Incidencia[n,l] * model.Trans_Theta[n,t,s] for n in model.Trans_Nodes)
    return model.Trans_Flow[l,t,s] - model.Trans_Status[l] * (1/model.Trans_Reactance[l]) * angle_difference_term == 0
model.Calculate_Trans_PowerFlow = pyo.Constraint(model.Trans_Lines, model.T, model.S, rule=calculate_trans_power_flow_rule)

# s.t. Trans_Power_Balance{n_node in Trans_Nodes, t in T, s in S}: 
# 	sum{g in G_T : G_Node[g] == n_node}(P_thermal_trans[g,t,s])
# 	-sum{l in Trans_Lines}(Trans_Incidencia[n_node,l]*Trans_Flow[l,t,s])
# 	-(if 5 == n_node then P_DSO[t,s]) 
# 	= Trans_Load[n_node,t,s];
def trans_power_balance_rule(model, n_node, t, s):
    gen_at_node = sum(model.P_thermal_trans[g,t,s] for g in model.G_T if model.G_Node[g] == n_node)
    
    net_flow_out_node = sum(model.Trans_Incidencia[n_node,l] * model.Trans_Flow[l,t,s] for l in model.Trans_Lines)
    
    p_dso_exchange = 0
    # Node 5 is hardcoded as the PCC with DSO in the AMPL model.
    # Trans_Nodes in input.dat are 1 to 14.
    if n_node == 5: # This should ideally be a parameter if it can change
        p_dso_exchange = model.P_DSO[t,s]
        
    return gen_at_node - net_flow_out_node - p_dso_exchange == model.Trans_Load[n_node,t,s]
model.Trans_Power_Balance = pyo.Constraint(model.Trans_Nodes, model.T, model.S, rule=trans_power_balance_rule)

# s.t. Trans_Flow_Upper{l in Trans_Lines, t in T, s in S}:  Trans_Flow[l,t,s] <= Trans_Capacity[l];
def trans_flow_upper_rule(model, l, t, s):
    return model.Trans_Flow[l,t,s] <= model.Trans_Capacity[l]
model.Trans_Flow_Upper = pyo.Constraint(model.Trans_Lines, model.T, model.S, rule=trans_flow_upper_rule)

# s.t. Trans_Flow_Lower{l in Trans_Lines, t in T, s in S}: -Trans_Flow[l,t,s] <= Trans_Capacity[l];
def trans_flow_lower_rule(model, l, t, s):
    return -model.Trans_Flow[l,t,s] <= model.Trans_Capacity[l] # or model.Trans_Flow[l,t,s] >= -model.Trans_Capacity[l]
model.Trans_Flow_Lower = pyo.Constraint(model.Trans_Lines, model.T, model.S, rule=trans_flow_lower_rule)

# s.t. Trans_Gen_Upper{g in G_T, t in T, s in S}:  P_thermal_trans[g,t,s] <=  Pmax_t[g];
def trans_gen_upper_rule(model, g, t, s):
    return model.P_thermal_trans[g,t,s] <= model.Pmax_t[g]
model.Trans_Gen_Upper = pyo.Constraint(model.G_T, model.T, model.S, rule=trans_gen_upper_rule)

# s.t. Trans_Gen_Lower{g in G_T, t in T, s in S}: -P_thermal_trans[g,t,s] <= -Pmin_t[g];
def trans_gen_lower_rule(model, g, t, s):
    return -model.P_thermal_trans[g,t,s] <= -model.Pmin_t[g] # or model.P_thermal_trans[g,t,s] >= model.Pmin_t[g]
model.Trans_Gen_Lower = pyo.Constraint(model.G_T, model.T, model.S, rule=trans_gen_lower_rule)

# s.t. Trans_Gen_Upper_Ramp{g in G_T, t in T, s in S: t<>first(T)}:  P_thermal_trans[g,t-1,s] - P_thermal_trans[g,t,s] <=  0.25 * Pmax_t[g];
def trans_gen_upper_ramp_rule(model, g, t, s):
    if t == model.T.first():
        return pyo.Constraint.Skip # No ramp constraint for the first period
    # This is Ramp Down: P_prev - P_curr <= MaxRampDown
    return model.P_thermal_trans[g,model.T.prev(t),s] - model.P_thermal_trans[g,t,s] <= 0.25 * model.Pmax_t[g]
model.Trans_Gen_Upper_Ramp = pyo.Constraint(model.G_T, model.T, model.S, rule=trans_gen_upper_ramp_rule)

# s.t. Trans_Gen_Lower_Ramp{g in G_T, t in T, s in S: t<>first(T)}: -P_thermal_trans[g,t-1,s] + P_thermal_trans[g,t,s] <=  0.25 * Pmax_t[g];
def trans_gen_lower_ramp_rule(model, g, t, s):
    if t == model.T.first():
        return pyo.Constraint.Skip # No ramp constraint for the first period
    # This is Ramp Up: P_curr - P_prev <= MaxRampUp
    return model.P_thermal_trans[g,t,s] - model.P_thermal_trans[g,model.T.prev(t),s] <= 0.25 * model.Pmax_t[g]
model.Trans_Gen_Lower_Ramp = pyo.Constraint(model.G_T, model.T, model.S, rule=trans_gen_lower_ramp_rule)

# s.t. MAXIMUM_SE_LIMIT {t in T, s in S}:    P_DSO[t,s] <=   SE_Capacity;
def maximum_se_limit_rule(model, t, s):
    return model.P_DSO[t,s] <= model.SE_Capacity
model.MAXIMUM_SE_LIMIT = pyo.Constraint(model.T, model.S, rule=maximum_se_limit_rule)
											
# s.t. MINIMUM_SE_LIMIT {t in T, s in S}:   -P_DSO[t,s] <=  SE_Capacity;
def minimum_se_limit_rule(model, t, s):
    return -model.P_DSO[t,s] <= model.SE_Capacity # or model.P_DSO[t,s] >= -model.SE_Capacity
model.MINIMUM_SE_LIMIT = pyo.Constraint(model.T, model.S, rule=minimum_se_limit_rule)

#---- CARBON RELATED CONSTRAINTS (Transmission)
# subject to Calculate_Footprint_Trans {g in G_T, t in T, s in S}:
#	Footprint_trans[g,t,s] = Carbon_Cost_t[g] * P_thermal_trans[g,t,s];		#phi
def calculate_footprint_trans_rule(model, g, t, s):
    return model.Footprint_trans[g,t,s] == model.Carbon_Cost_t[g] * model.P_thermal_trans[g,t,s]
model.Calculate_Footprint_Trans = pyo.Constraint(model.G_T, model.T, model.S, rule=calculate_footprint_trans_rule)

# s.t. Carbon_Balance_Trans {t in T, s in S}:
#	sum{g in G_T} (Carbon_T[g,t,s]) - Carbon_SE[t,s] = 0;					#psi
def carbon_balance_trans_rule(model, t, s):
    sum_carbon_t = sum(model.Carbon_T[g,t,s] for g in model.G_T)
    return sum_carbon_t - model.Carbon_SE[t,s] == 0
model.Carbon_Balance_Trans = pyo.Constraint(model.T, model.S, rule=carbon_balance_trans_rule)

# s.t. Carbon_Limits_Trans {t in T, s in S}:									
#	sum{g in G_T} (Footprint_trans[g,t,s] + Carbon_T[g,t,s]) <= Carbon_Limit_Trans[t] ;		#Xi
def carbon_limits_trans_rule(model, t, s):
    sum_total_carbon_footprint = sum(model.Footprint_trans[g,t,s] + model.Carbon_T[g,t,s] for g in model.G_T)
    return sum_total_carbon_footprint <= model.Carbon_Limit_Trans[t]
model.Carbon_Limits_Trans = pyo.Constraint(model.T, model.S, rule=carbon_limits_trans_rule)

#-----------------------------------------------------------------------
# TRANSMISSION DUAL CONSTRAINTS (KKT Conditions)
#-----------------------------------------------------------------------

# s.t. Deriv_Potencia_with_Ramp {g in G_T, t in T, s in S}:
#	(2*a_t[g] * P_thermal_trans[g,t,s] + b_t[g]) - lambda[G_Node[g],t,s] 
#	-phi[g,t,s]*Carbon_Cost_t[g] 
#	- alpha_L[g,t,s] + alpha_U[g,t,s] 
#	+(if t<>last(T) then (ramp_up[g,t+1,s] - ramp_low[g,t+1,s]))
#	+(if t<>first(T) then (-ramp_up[g,t,s] + ramp_low[g,t,s]))
#	= 0;
def deriv_potencia_with_ramp_rule(model, g, t, s):
    expr = (2 * model.a_t[g] * model.P_thermal_trans[g,t,s] + model.b_t[g]) - \
           model.lambda_kkt[model.G_Node[g],t,s] - \
           model.phi_kkt[g,t,s] * model.Carbon_Cost_t[g] - \
           model.alpha_L[g,t,s] + model.alpha_U[g,t,s]
    
    if t != model.T.last():
        expr += (model.ramp_up[g,model.T.next(t),s] - model.ramp_low[g,model.T.next(t),s])
    if t != model.T.first():
        expr += (-model.ramp_up[g,t,s] + model.ramp_low[g,t,s])
    return expr == 0
model.Deriv_Potencia_with_Ramp = pyo.Constraint(model.G_T, model.T, model.S, rule=deriv_potencia_with_ramp_rule)

# s.t. Deriv_Potencia_without_Ramp {g in G_T, t in T, s in S}:
#	(2*a_t[g] * P_thermal_trans[g,t,s] + b_t[g]) - lambda[G_Node[g],t,s] 
#	-phi[g,t,s]*Carbon_Cost_t[g] 
#	- alpha_L[g,t,s] + alpha_U[g,t,s] 
#	= 0;
def deriv_potencia_without_ramp_rule(model, g, t, s):
    return (2 * model.a_t[g] * model.P_thermal_trans[g,t,s] + model.b_t[g]) - \
           model.lambda_kkt[model.G_Node[g],t,s] - \
           model.phi_kkt[g,t,s] * model.Carbon_Cost_t[g] - \
           model.alpha_L[g,t,s] + model.alpha_U[g,t,s] == 0
model.Deriv_Potencia_without_Ramp = pyo.Constraint(model.G_T, model.T, model.S, rule=deriv_potencia_without_ramp_rule)

# s.t. Deriv_Teta{n in Trans_Nodes, t in T, s in S}: -sum{l in Trans_Lines}(
#	-Trans_Status[l] * (1/Trans_Reactance[l]) * Trans_Incidencia[n,l]*mu[l,t,s]) = 0;
def deriv_teta_rule(model, n_node, t, s):
    # Handle Trans_Reactance[l] == 0 carefully.
    # Original term: Trans_Status[l] * (1/Trans_Reactance[l]) * Trans_Incidencia[n,l] * mu[l,t,s]
    # If Trans_Reactance[l] is 0, this term is problematic.
    # The primal constraint was: Trans_Reactance[l] * Trans_Flow[l,t,s] - Trans_Status[l] * sum_incid_theta = 0
    # Derivative of L w.r.t Trans_Theta[n,t,s] involves:
    # - mu[l,t,s] * Trans_Status[l] * Trans_Incidencia[n,l]
    # This seems to be the coefficient of Trans_Theta[n,t,s] in the power flow equation, multiplied by mu.
    # The AMPL KKT is: sum{l in Trans_Lines} (Trans_Status[l] * (1/Trans_Reactance[l]) * Trans_Incidencia[n_node,l] * model.mu_kkt[l,t,s]) == 0
    # This assumes Trans_Reactance is not zero.
    sum_terms = 0
    for l_line in model.Trans_Lines:
        if model.Trans_Reactance[l_line] != 0:
            sum_terms += (model.Trans_Status[l_line] * (1/model.Trans_Reactance[l_line]) * 
                          model.Trans_Incidencia[n_node,l_line] * model.mu_kkt[l_line,t,s])
        # If Trans_Reactance is 0, and Trans_Status is 1, the original constraint implies angle equality.
        # The derivative term might be handled differently or assumed that such lines don't exist if active.
        # For now, skip if reactance is 0 to avoid division by zero, assuming such lines are inactive or handled.
    return sum_terms == 0
model.Deriv_Teta = pyo.Constraint(model.Trans_Nodes, model.T, model.S, rule=deriv_teta_rule)

# s.t. Deriv_Fluxo{l in Trans_Lines, t in T, s in S}: mu[l,t,s] - omega_L[l,t,s] + omega_U[l,t,s] 
#	-sum{n in Trans_Nodes}(Trans_Incidencia[n,l] * lambda[n,t,s]) = 0;
def deriv_fluxo_rule(model, l_line, t, s):
    # The term sum{n in Trans_Nodes}(Trans_Incidencia[n,l] * lambda[n,t,s]) is:
    # Trans_Incidencia[From[l],l]*lambda[From[l],t,s] + Trans_Incidencia[To[l],l]*lambda[To[l],t,s]
    # = 1 * lambda[From[l],t,s] - 1 * lambda[To[l],t,s]
    # = lambda[model.From[l_line],t,s] - lambda[model.To[l_line],t,s]
    # This is the derivative of lambda_from*Flow - lambda_to*Flow if Flow is defined from->to.
    # The Trans_Power_Balance has -sum{l}(Trans_Incidencia[n,l]*Trans_Flow[l,t,s]).
    # For a given l, Trans_Flow[l,t,s] is multiplied by -Trans_Incidencia[n,l].
    # So, for node From[l], it's -1*Trans_Flow[l,t,s]. For node To[l], it's -(-1)*Trans_Flow[l,t,s] = +Trans_Flow[l,t,s].
    # Derivative of Lagrangian w.r.t. Trans_Flow[l,t,s]:
    # from mu_kkt: +mu_kkt[l,t,s] (from Calculate_Trans_PowerFlow)
    # from omega_L: -omega_L[l,t,s] (from -Trans_Flow <= Capacity => Trans_Flow >= -Capacity)
    # from omega_U: +omega_U[l,t,s] (from Trans_Flow <= Capacity)
    # from lambda_kkt: -lambda_kkt[model.From[l_line],t,s] + lambda_kkt[model.To[l_line],t,s]
    # The AMPL is: mu - omega_L + omega_U - (lambda_from - lambda_to) = 0
    # This matches: mu - omega_L + omega_U - lambda_from + lambda_to = 0
    
    lambda_term = model.lambda_kkt[model.From[l_line],t,s] - model.lambda_kkt[model.To[l_line],t,s]
    
    return model.mu_kkt[l_line,t,s] - model.omega_L[l_line,t,s] + model.omega_U[l_line,t,s] - lambda_term == 0
model.Deriv_Fluxo = pyo.Constraint(model.Trans_Lines, model.T, model.S, rule=deriv_fluxo_rule)

# s.t. Deriv_PDSO{t in T, s in S}:  
#	-Bid[t,s] + lambda[5,t,s] + kappa_U[t,s] - kappa_L[t,s] = 0;
# Node 5 is PCC for P_DSO in Trans_Power_Balance: -(if 5 == n then P_DSO[t,s])
# So derivative of -lambda[5,t,s]*P_DSO[t,s] is -lambda[5,t,s].
# The objective has Bid[t,s]*P_DSO[t,s]. Derivative is Bid[t,s].
# kappa_U for P_DSO <= SE_Capacity. Derivative +kappa_U.
# kappa_L for -P_DSO <= SE_Capacity (P_DSO >= -SE_Capacity). Derivative -kappa_L.
# So: Bid[t,s] - lambda[5,t,s] + kappa_U[t,s] - kappa_L[t,s] = 0.
# AMPL has: -Bid[t,s] + lambda[5,t,s] + kappa_U[t,s] - kappa_L[t,s] = 0.
# This implies the objective for the lower level problem (ISO) is to MINIMIZE -Bid*P_DSO + GenCosts_trans ...
# Or the DSO's objective term Bid*P_DSO is taken as a negative cost for ISO.
# Let's assume the AMPL KKT is derived correctly from some formulation.
# The node 5 should be parameterized if possible. For now, hardcoding.
PCC_Trans_Node = 5 # As per AMPL model
def deriv_pdso_rule(model, t, s):
    return -model.Bid[t,s] + model.lambda_kkt[PCC_Trans_Node,t,s] + model.kappa_U[t,s] - model.kappa_L[t,s] == 0
model.Deriv_PDSO = pyo.Constraint(model.T, model.S, rule=deriv_pdso_rule)

# s.t. Deriv_Footprint{g in G_T, t in T, s in S}:
#	+phi[g,t,s]  #- epsilon[g,t,s] 
#	+ Xi[t,s] = 0;
def deriv_footprint_rule(model, g, t, s):
    expr = model.phi_kkt[g,t,s] + model.Xi_kkt[t,s]
    # if model.epsilon_kkt is part of the model and active:
    # expr -= model.epsilon_kkt[g,t,s] 
    return expr == 0
model.Deriv_Footprint = pyo.Constraint(model.G_T, model.T, model.S, rule=deriv_footprint_rule)

# s.t. Deriv_Carbon_T {g in G_T, t in T, s in S}:
#	#-epsilon[g,t,s] 
#	+ psi[t,s] + Xi[t,s] = 0;
def deriv_carbon_t_rule(model, g, t, s):
    expr = model.psi_kkt[t,s] + model.Xi_kkt[t,s]
    # if model.epsilon_kkt is part of the model and active:
    # expr -= model.epsilon_kkt[g,t,s]
    return expr == 0
model.Deriv_Carbon_T = pyo.Constraint(model.G_T, model.T, model.S, rule=deriv_carbon_t_rule)
	
# s.t. Deriv_Carbon_SE {t in T, s in S}:
#	 -Carbon_Price[t,s] -psi[t,s] 
#	 # + eta_U[t,s] - eta_L[t,s] 
#	 = 0;
# Objective has Carbon_Price[t,s] * Carbon_SE[t,s].
# Carbon_Balance_Trans: sum(Carbon_T) - Carbon_SE = 0. Dual psi_kkt. Term: -psi_kkt * Carbon_SE.
# So derivative w.r.t Carbon_SE is: Carbon_Price - psi_kkt = 0.
# AMPL has: -Carbon_Price - psi = 0. This implies the objective term for ISO was -Carbon_Price*Carbon_SE or psi_kkt definition is negative.
# Assuming AMPL KKT is correct.
def deriv_carbon_se_rule(model, t, s):
    expr = -model.Carbon_Price[t,s] - model.psi_kkt[t,s]
    # if model.eta_U_kkt and model.eta_L_kkt are active:
    # expr += model.eta_U_kkt[t,s] - model.eta_L_kkt[t,s]
    return expr == 0
model.Deriv_Carbon_SE = pyo.Constraint(model.T, model.S, rule=deriv_carbon_se_rule)

# Non-negativity of duals already handled by domain=pyo.NonNegativeReals in Var declaration.
# Explicit constraints for these are not strictly needed in Pyomo unless for conditional dropping.
# e.g. model.Omega_Lower = pyo.Constraint(model.Trans_Lines, model.T, model.S, rule=lambda model,l,t,s: model.omega_L[l,t,s] >=0)

#-----------------------------------------------------------------------
# LINEAR KKT COMPLEMENTARY SLACKNESS CONDITIONS
#-----------------------------------------------------------------------
# s.t. slack_NL2 {l in Trans_Lines, t in T, s in S}: omega_U[l,t,s]*(-Trans_Capacity[l] + Trans_Flow[l,t,s]) = 0;
def slack_nl2_rule(model, l, t, s):
    return model.omega_U[l,t,s] * (-model.Trans_Capacity[l] + model.Trans_Flow[l,t,s]) == 0
model.slack_NL2 = pyo.Constraint(model.Trans_Lines, model.T, model.S, rule=slack_nl2_rule)

# s.t. slack_NL1 {l in Trans_Lines, t in T, s in S}: omega_L[l,t,s]*(-Trans_Flow[l,t,s] - Trans_Capacity[l]) = 0;
def slack_nl1_rule(model, l, t, s):
    return model.omega_L[l,t,s] * (-model.Trans_Flow[l,t,s] - model.Trans_Capacity[l]) == 0
model.slack_NL1 = pyo.Constraint(model.Trans_Lines, model.T, model.S, rule=slack_nl1_rule)

# s.t. slack_NL4 {g in G_T, t in T, s in S}: alpha_U[g,t,s]*(P_thermal_trans[g,t,s] - Pmax_t[g]) = 0;
def slack_nl4_rule(model, g, t, s):
    return model.alpha_U[g,t,s] * (model.P_thermal_trans[g,t,s] - model.Pmax_t[g]) == 0
model.slack_NL4 = pyo.Constraint(model.G_T, model.T, model.S, rule=slack_nl4_rule)

# s.t. slack_NL3 {g in G_T, t in T, s in S}: alpha_L[g,t,s]*(-P_thermal_trans[g,t,s] + Pmin_t[g]) = 0;
def slack_nl3_rule(model, g, t, s):
    return model.alpha_L[g,t,s] * (-model.P_thermal_trans[g,t,s] + model.Pmin_t[g]) == 0
model.slack_NL3 = pyo.Constraint(model.G_T, model.T, model.S, rule=slack_nl3_rule)

# s.t. slack_NL5 {t in T, s in S}: kappa_U[t,s]*(-SE_Capacity + P_DSO[t,s]) = 0;
def slack_nl5_rule(model, t, s):
    return model.kappa_U[t,s] * (-model.SE_Capacity + model.P_DSO[t,s]) == 0
model.slack_NL5 = pyo.Constraint(model.T, model.S, rule=slack_nl5_rule)

# s.t. slack_NL6 {t in T, s in S}: kappa_L[t,s]*(-P_DSO[t,s] - SE_Capacity) = 0;
def slack_nl6_rule(model, t, s):
    return model.kappa_L[t,s] * (-model.P_DSO[t,s] - model.SE_Capacity) == 0
model.slack_NL6 = pyo.Constraint(model.T, model.S, rule=slack_nl6_rule)

# s.t. slack_NL7 {t in T, s in S} :Xi[t,s]*(sum{g in G_T} (Footprint_trans[g,t,s] + Carbon_T[g,t,s]) - Carbon_Limit_Trans[t] ) = 0;
def slack_nl7_rule(model, t, s):
    sum_term = sum(model.Footprint_trans[g,t,s] + model.Carbon_T[g,t,s] for g in model.G_T)
    return model.Xi_kkt[t,s] * (sum_term - model.Carbon_Limit_Trans[t]) == 0
model.slack_NL7 = pyo.Constraint(model.T, model.S, rule=slack_nl7_rule)

# slack_NL8 for epsilon is commented in AMPL

# slack_NL9, slack_NL10 for eta are commented in AMPL

# s.t. slack_NL11{g in G_T, t in T, s in S : t<>first(T)}: ramp_up[g,t,s]*(P_thermal_trans[g,t-1,s] - P_thermal_trans[g,t,s] -0.25 * Pmax_t[g]) = 0;
def slack_nl11_rule(model, g, t, s):
    if t == model.T.first():
        return pyo.Constraint.Skip
    return model.ramp_up[g,t,s] * (model.P_thermal_trans[g,model.T.prev(t),s] - model.P_thermal_trans[g,t,s] - 0.25 * model.Pmax_t[g]) == 0
model.slack_NL11 = pyo.Constraint(model.G_T, model.T, model.S, rule=slack_nl11_rule)

# s.t. slack_NL12{g in G_T, t in T, s in S : t<>first(T)}: ramp_low[g,t,s]*(-P_thermal_trans[g,t-1,s] + P_thermal_trans[g,t,s] -0.25 * Pmax_t[g]) = 0;
def slack_nl12_rule(model, g, t, s):
    if t == model.T.first():
        return pyo.Constraint.Skip
    return model.ramp_low[g,t,s] * (-model.P_thermal_trans[g,model.T.prev(t),s] + model.P_thermal_trans[g,t,s] - 0.25 * model.Pmax_t[g]) == 0
model.slack_NL12 = pyo.Constraint(model.G_T, model.T, model.S, rule=slack_nl12_rule)


# TODO: Implement data loading from input.dat and Scenarios.dat
# TODO: Implement logic from execute.run (parameter adjustments, conditional constraint dropping)
# A leitura dos dados (equivalente ao input.dat) será tratada posteriormente.
# O arquivo execute.run também contém lógica que precisará ser traduzida para Python.

# Exemplo de como os dados seriam carregados (placeholder):
# data_input_dat = {
#    None: {
#        'N': {None: [1,2,3,...]}, # Lista de nós
#        'G_T': {None: [101,102,...]}, # Lista de geradores de transmissão
#        # ... outros conjuntos e parâmetros
#    }
#}
# instance = model.create_instance(data_input_dat)
# No entanto, para modelos grandes, é mais comum carregar dados de arquivos CSV ou Excel usando Pandas.

print("Modelo Pyomo com parâmetros, variáveis e objetivo definidos. Próximo passo: traduzir restrições.")
