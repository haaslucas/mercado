from pyomo.environ import *

# Criar um modelo
model = ConcreteModel()

## Definições de Parâmetros
model.VBASE = Param(initialize=24.90)  # tensão base (kV)
model.SBASE = Param(initialize=100)     # potência aparente base (MVA)
model.IBASE = model.SBASE * 1e3 / model.VBASE
model.SE_Capacity = Param(initialize=50)

# Conjuntos
model.N = Set()
model.G_T = Set()
model.G_D = Set()
model.LS = Set()
model.ESS = Set()
model.RES_T = Set()
model.RES_D = Set()
model.L = Set(initialize=lambda model: [(i,j) for i in model.N for j in model.N])  # Defina as conexões apropriadas
model.T = Set(ordered=True)
model.Trans_Lines = Set()
model.Trans_Nodes = Set()
model.LDA = Set()
model.WS = Set(initialize=lambda model: [(1, n, g) for n in model.Trans_Nodes for g in model.G_T.union(model.LDA)])  # Ajuste conforme a necessidade

## Parâmetros lidos do arquivo de dados
model.NBUS = Param(model.N, initialize={})  # Dicionário com os nomes das buses
model.TBUS = Param(model.N, initialize={})   # Tipo de bus
model.Vmax = Param(model.N, initialize={})
model.Vmin = Param(model.N, initialize={})
model.Vnom = Param(model.N, initialize={})
model.Ve = Param(model.N, initialize={})
model.Pgmx = Param(model.N, initialize={})
model.Pgmn = Param(model.N, initialize={})

# Adicione os parâmetros restantes conforme necessário, seguindo a mesma estrutura

# Variáveis
model.VM = Var(model.N, model.T, within=NonNegativeReals)        # Tensão
model.I = Var(model.L, model.T, within=NonNegativeReals)         # Corrente
model.Trans_Shift = Var(model.Trans_Nodes, model.T, within=Reals)  # Carregamento da transmissão
# Defina as demais variáveis da mesma forma


# Função objetivo
def objective_rule(model):
    return sum((model.Bid[t, s] * model.P_DSO[t, s] +
                sum(model.GenCosts_dist[g, t, s] for g in model.G_D) +
                (model.Carbon_Price[t, s] * model.Carbon_SE[t, s]) 
                for t in model.T for s in model.S))

model.Objective = Objective(rule=objective_rule, sense=minimize)

# Restrições
def active_power_balance_rule(model, n, t, s):
    return (sum(model.P_thermal_dist[g, t, s] for g in model.G_D if model.G_LDA_Node[g] == n) -
            sum(model.P[i, j, t, s] + model.R[i, j] * model.I[i, j, t, s] for (i, j) in model.L) +
            sum(model.P[l, n, t, s] for l in model.L) +
            (model.P_DSO[t, s] if n == 1 else 0) -
            sum(model.Dist_Shift[i, t, s] for i in model.LS if model.LS_Node[i] == n) +
            sum(model.Trans_P_ESS[i, t, s] for i in model.ESS if model.ESS_Node[i] == n) 
            == model.Dist_Load[n, t, s])

model.Active_Power_Balance = Constraint(model.N, model.T, model.S, rule=active_power_balance_rule)

# Adicione as demais restrições da mesma forma

# Resolve o modelo
solver = SolverFactory('glpk')  # ou outro solver disponível
results = solver.solve(model, tee=True)

# Exibir os resultados
model.display()