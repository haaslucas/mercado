from pyomo.environ import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Criação do modelo
model = ConcreteModel()

# Conjuntos
model.i = RangeSet(1, 24)  # network buses
model.slack = 13 # slack bus
model.t = RangeSet(1, 24)  # time periods
model.GB = Set(initialize=[1, 2, 7, 13, 15, 16, 18, 21, 22, 23])  # generating buses

# Parâmetros
model.Sbase = Param(initialize=100)
model.VOLL = Param(initialize=10000)
model.VOLW = Param(initialize=50)

# Características das unidades geradoras
GenD = {
    1: {'pmax': 152, 'pmin': 30.4, 'b': 13.32, 'Qmax': 192, 'Qmin': -50, 'Vg': 1.035, 'RU': 21, 'RD': 21},
    2: {'pmax': 152, 'pmin': 30.4, 'b': 13.32, 'Qmax': 192, 'Qmin': -50, 'Vg': 1.035, 'RU': 21, 'RD': 21},
    7: {'pmax': 350, 'pmin': 75, 'b': 20.7, 'Qmax': 300, 'Qmin': 0, 'Vg': 1.025, 'RU': 43, 'RD': 43},
    13: {'pmax': 591, 'pmin': 206.85, 'b': 20.93, 'Qmax': 591, 'Qmin': 0, 'Vg': 1.02, 'RU': 31, 'RD': 31},
    15: {'pmax': 215, 'pmin': 66.25, 'b': 21, 'Qmax': 215, 'Qmin': -100, 'Vg': 1.014, 'RU': 31, 'RD': 31},
    16: {'pmax': 155, 'pmin': 54.25, 'b': 10.52, 'Qmax': 155, 'Qmin': -50, 'Vg': 1.017, 'RU': 31, 'RD': 31},
    18: {'pmax': 400, 'pmin': 100, 'b': 5.47, 'Qmax': 400, 'Qmin': -50, 'Vg': 1.05, 'RU': 70, 'RD': 70},
    21: {'pmax': 400, 'pmin': 100, 'b': 5.47, 'Qmax': 400, 'Qmin': -50, 'Vg': 1.05, 'RU': 70, 'RD': 70},
    22: {'pmax': 300, 'pmin': 0, 'b': 0, 'Qmax': 300, 'Qmin': -60, 'Vg': 1.05, 'RU': 53, 'RD': 53},
    23: {'pmax': 360, 'pmin': 248.5, 'b': 10.52, 'Qmax': 310, 'Qmin': -125, 'Vg': 1.05, 'RU': 31, 'RD': 31}
}



model.GenD = Param(model.GB, ['pmax', 'pmin', 'b', 'Qmax', 'Qmin', 'Vg', 'RU', 'RD'], initialize=lambda model, i, j: GenD[i][j])

# Demandas de cada bus em MW e MVAr
BD = {
    1: {'Pd': 108, 'Qd': 22},
    2: {'Pd': 97, 'Qd': 20},
    3: {'Pd': 180, 'Qd': 37},
    4: {'Pd': 74, 'Qd': 15},
    5: {'Pd': 71, 'Qd': 14},
    6: {'Pd': 136, 'Qd': 28},
    7: {'Pd': 125, 'Qd': 25},
    8: {'Pd': 171, 'Qd': 35},
    9: {'Pd': 175, 'Qd': 36},
    10: {'Pd': 195, 'Qd': 40},
    11: {'Pd': 0, 'Qd': 0},
    12: {'Pd': 0, 'Qd': 0},
    13: {'Pd': 265, 'Qd': 54},
    14: {'Pd': 194, 'Qd': 39},
    15: {'Pd': 317, 'Qd': 64},
    16: {'Pd': 100, 'Qd': 20},
    17: {'Pd': 0, 'Qd': 0},
    18: {'Pd': 333, 'Qd': 68},
    19: {'Pd': 181, 'Qd': 37},
    20: {'Pd': 128, 'Qd': 26},
    21: {'Pd': 0, 'Qd': 0},
    22: {'Pd': 0, 'Qd': 0},
    23: {'Pd': 0, 'Qd': 0},
    24: {'Pd': 0, 'Qd': 0},
}

model.BD = Param(model.i, ['Pd', 'Qd'], initialize=lambda model, i, j: BD[i][j])


LN = { # Características da rede
    (1, 2): {'r': 0.0026, 'x': 0.0139, 'b': 0.4611, 'limit': 175},
    (1, 3): {'r': 0.0546, 'x': 0.2112, 'b': 0.0572, 'limit': 175},
    (1, 5): {'r': 0.0218, 'x': 0.0845, 'b': 0.0229, 'limit': 175},
    (2, 4): {'r': 0.0328, 'x': 0.1267, 'b': 0.0343, 'limit': 175},
    (2, 6): {'r': 0.0497, 'x': 0.1920, 'b': 0.0520, 'limit': 175},
    (3, 9): {'r': 0.0308, 'x': 0.1190, 'b': 0.0322, 'limit': 175},
    (3, 24): {'r': 0.0023, 'x': 0.0839, 'b': 0.0, 'limit': 400},
    (4, 9): {'r': 0.0268, 'x': 0.1037, 'b': 0.0281, 'limit': 175},
    (5, 10): {'r': 0.0228, 'x': 0.0883, 'b': 0.0239, 'limit': 175},
    (6, 10): {'r': 0.0139, 'x': 0.0605, 'b': 2.459, 'limit': 175},
    (7, 8): {'r': 0.0159, 'x': 0.0614, 'b': 0.0166, 'limit': 175},
    (8, 9): {'r': 0.0427, 'x': 0.1651, 'b': 0.0447, 'limit': 175},
    (8, 10): {'r': 0.0427, 'x': 0.1651, 'b': 0.0447, 'limit': 175},
    (9, 11): {'r': 0.0023, 'x': 0.0839, 'b': 0.0, 'limit': 400},
    (9, 12): {'r': 0.0023, 'x': 0.0839, 'b': 0.0, 'limit': 400},
    (10, 11): {'r': 0.0023, 'x': 0.0839, 'b': 0.0, 'limit': 400},
    (10, 12): {'r': 0.0023, 'x': 0.0839, 'b': 0.0, 'limit': 400},
    (11, 13): {'r': 0.0061, 'x': 0.0476, 'b': 0.0999, 'limit': 500},
    (11, 14): {'r': 0.0054, 'x': 0.0418, 'b': 0.0879, 'limit': 500},
    (12, 13): {'r': 0.0061, 'x': 0.0476, 'b': 0.0999, 'limit': 500},
    (12, 23): {'r': 0.0124, 'x': 0.0966, 'b': 0.2030, 'limit': 500},
    (13, 23): {'r': 0.0111, 'x': 0.0865, 'b': 0.1818, 'limit': 500},
    (14, 16): {'r': 0.0050, 'x': 0.0389, 'b': 0.0818, 'limit': 500},
    (15, 16): {'r': 0.0022, 'x': 0.0173, 'b': 0.0364, 'limit': 500},
    (15, 21): {'r': 0.00315, 'x': 0.0245, 'b': 0.2060, 'limit': 1000},
    (15, 24): {'r': 0.0067, 'x': 0.0519, 'b': 0.1091, 'limit': 500},
    (16, 17): {'r': 0.0033, 'x': 0.0259, 'b': 0.0545, 'limit': 500},
    (16, 19): {'r': 0.0030, 'x': 0.0231, 'b': 0.0485, 'limit': 500},
    (17, 18): {'r': 0.0018, 'x': 0.0144, 'b': 0.0303, 'limit': 500},
    (17, 22): {'r': 0.0135, 'x': 0.1053, 'b': 0.2212, 'limit': 500},
    (18, 21): {'r': 0.00165, 'x': 0.01295, 'b': 0.1090, 'limit': 1000},
    (19, 20): {'r': 0.00255, 'x': 0.0198, 'b': 0.1666, 'limit': 1000},
    (20, 23): {'r': 0.0014, 'x': 0.0108, 'b': 0.0910, 'limit': 1000},
    (21, 22): {'r': 0.0087, 'x': 0.0678, 'b': 0.1424, 'limit': 500},
}
#calculate 'th' for each line
for (i, j), params in LN.items():
    r, x, b, limit = params['r'], params['x'], params['b'], params['limit']
    LN[(i, j)]['th'] = np.arctan2(x, r)  # angle in radians # np.arctan(x/r) 
    LN[(i, j)]['z'] = np.sqrt(r**2 + x**2)  # impedance magnitude

model.LN = Param(model.i, model.i, ['r', 'x', 'b', 'limit','z','th'], initialize=lambda model, i, j, k: LN.get((i, j), {}).get(k, 0))

# Dados de geração eólica e demanda
WD = {
    't1': {'w': 0.0786666666666667, 'd': 0.684511335492475},
    't2': {'w': 0.0866666666666667, 'd': 0.644122690036197},
    't3': {'w': 0.117333333333333, 'd': 0.6130691560297},
    't4': {'w': 0.258666666666667, 'd': 0.599733282530006},
    't5': {'w': 0.361333333333333, 'd': 0.588874071251667},
    't6': {'w': 0.566666666666667, 'd': 0.5980186702229},
    't7': {'w': 0.650666666666667, 'd': 0.626786054486569},
    't8': {'w': 0.566666666666667, 'd': 0.651743189178891},
    't9': {'w': 0.484, 'd': 0.706039245570585},
    't10': {'w': 0.548, 'd': 0.787007048961707},
    't11': {'w': 0.757333333333333, 'd': 0.839016955610593},
    't12': {'w': 0.710666666666667, 'd': 0.852733854067441},
    't13': {'w': 0.870666666666667, 'd': 0.870642027052772},
    't14': {'w': 0.932, 'd': 0.834254143646409},
    't15': {'w': 0.966666666666667, 'd': 0.816536483139646},
    't16': {'w': 1, 'd': 0.819394170318156},
    't17': {'w': 0.869333333333333, 'd': 0.874071251666984},
    't18': {'w': 0.665333333333333, 'd': 1},
    't19': {'w': 0.656, 'd': 0.983615926843208},
    't20': {'w': 0.561333333333333, 'd': 0.936368832158506},
    't21': {'w': 0.565333333333333, 'd': 0.887597637645266},
    't22': {'w': 0.556, 'd': 0.809297008954087},
    't23': {'w': 0.724, 'd': 0.74585635359116},
    't24': {'w': 0.84, 'd': 0.733473042484283},
}
model.WD = Param(model.t, ['w', 'd'], initialize=lambda model, t, j: WD[f't{t}'][j])



# Capacidade da energia eólica
Wcap = {8: 200, 19: 150, 21: 100}
model.Wcap = Param(model.i, initialize=lambda model, i: Wcap.get(i, 0))

# Variáveis

#model.Pij = Var(model.i, model.i, model.t)  # Power flow variables
#model.Qij = Var(model.i, model.i, model.t)  # Reactive power flow variables

def Pij_bounds(model, i, j, t):
    if (i, j) in LN:
        return (-LN[i, j]['limit'] / model.Sbase, LN[i, j]['limit'] / model.Sbase)
    return (0, 0)  # No flow for non-existent lines
model.Pij = Var(model.i, model.i, model.t, within=Reals, bounds=Pij_bounds)  # Active power flow
# Aplicando limites às variáveis com lambda

def Qij_bounds(model, i, j, t):
    if (i, j) in LN:
        return (-LN[i, j]['limit'] / model.Sbase, LN[i, j]['limit'] / model.Sbase)
    return (0, 0)  # No flow for non-existent lines
model.Qij = Var(model.i, model.i, model.t, within=Reals, bounds=Qij_bounds)  # Reactive power flow

def Pg_bounds(model, i, t):
    return (model.GenD[i,'pmin'] / model.Sbase, model.GenD[i,'pmax'] / model.Sbase)
model.Pg = Var(model.GB, model.t, within=NonNegativeReals, bounds=Pg_bounds)  # Active generation

def Qg_bounds(model, i, t):
    return (model.GenD[i,'Qmin'] / model.Sbase, model.GenD[i,'Qmax'] / model.Sbase)
model.Qg = Var(model.GB, model.t, within=Reals, bounds=Qg_bounds)
#model.Qg = Var(model.GB, model.t, within=NonNegativeReals)  # Reactive generation

def Va_bounds(model, i, t):
    if i == model.slack:
        return (0, 0)  # Slack bus angle is fixed at 0
    return (-np.pi/2, np.pi/2)  # Voltage angle bounds for other buses
model.Va = Var(model.i, model.t, bounds=Va_bounds, initialize=lambda model, i, t: 0 if i == model.slack else 0)  # Voltage angle
#fixando o ângulo do slack bus

#model.Va[model.slack, :] = 0  # Define o ângulo para o slack bus

def V_bounds(model, i, t):
    if i == model.slack:
        return (1, 1)  # Slack bus voltage is fixed at 1.035 p.u.
    return (0.9, 1.1)  # Voltage magnitude bounds for other buses
model.V = Var(model.i, model.t, bounds=V_bounds, initialize=lambda model, i, t: 1 if i == model.slack else 1.0)  # Voltage magnitude

def Pw_bounds(model, i, t):
    if i in model.Wcap:
        return (0, model.WD[t,'w'] * model.Wcap[i] / model.Sbase)  # Wind power bounds
    return (0, 0)  # No wind power for other buses
model.Pw = Var(model.i, model.t, within=NonNegativeReals, bounds=Pw_bounds)  # Wind power

# Equações
def eq1(model, i, j, t): # Eq. 6.8f - Active power flow of lines
    if (i, j) in LN:
        return model.Pij[i, j, t] == (model.V[i, t]**2 * cos(model.LN[i, j, 'th']) 
                                        - model.V[i, t] * model.V[j, t] * cos(model.Va[i, t] 
                                        - model.Va[j, t] + model.LN[i, j, 'th'])) / model.LN[i, j, 'z']
    return Constraint.Skip

def eq2(model, i, j, t): # Eq. 6.8g - Reactive power flow of lines
    if (i, j) in LN: 
        return model.Qij[i, j, t] == (model.V[i, t]**2 * sin(model.LN[i, j, 'th'])
                    - model.V[i, t] * model.V[j, t] * sin(model.Va[i, t] - model.Va[j, t] 
                    + model.LN[i, j, 'th'])) / model.LN[i, j, 'z'] - model.LN[i, j, 'b'] * model.V[i, t]**2 / 2
    return Constraint.Skip

def eq3(model, i, t):  # Eq.6.8b - Active power balance
    lhs = model.Pw[i, t] + (model.Pg[i, t] if i in model.GB else 0) - model.WD[t, 'd'] * model.BD[i, 'Pd'] / model.Sbase
    rhs = sum(model.Pij[j, i, t] for j in model.i if (i, j) in LN)
    
    expression = (lhs == rhs)
    if expression is True:  # Check if the expression is literally the boolean True
        return Constraint.Feasible
    return expression

def eq4(model, i, t): # Eq.6.8c - Reactive power balance
    lhs = (model.Qg[i, t] if i in model.GB else 0) - model.WD[t, 'd'] * model.BD[i, 'Qd'] / model.Sbase
    rhs = sum(model.Qij[j, i, t] for j in model.i if (i, j) in LN)

    expression = (lhs == rhs)
    if expression is True:  # Check if the expression is literally the boolean True
        return Constraint.Feasible
    return expression

'''def eq5(model):
    model.OF == sum(model.Pg[i, t] * model.GenD[i, 'b'] * 100 for i in model.GB for t in model.t)
    return model.OF'''

model.OF = Objective(
    expr=sum(model.Pg[i, t] * model.GenD[i, 'b'] * model.Sbase for i in model.GB for t in model.t),
    sense=minimize)


def eq6(model, i, t):
    if model.GenD[i,'pmax'] and t > 1:
        return model.Pg[i, t] - model.Pg[i, t-1] <= model.GenD[i,'RU'] / model.Sbase
    return Constraint.Skip

def eq7(model, i, t):
    if model.GenD[i,'pmax'] and t < (len(model.t)):
        return model.Pg[i, t] - model.Pg[i, t+1] <= model.GenD[i,'RD'] / model.Sbase
    return Constraint.Skip

model.Eq1 = Constraint(model.i, model.i, model.t, rule=eq1) # active power flow of lines
model.Eq2 = Constraint(model.i, model.i, model.t, rule=eq2) # reactive power flow of lines
model.Eq3 = Constraint(model.i, model.t, rule=eq3) # active power balance; Optimal found, = 0
model.Eq4 = Constraint(model.i, model.t, rule=eq4) # reactive power balance
#model.Eq5 = Objective(rule=eq5, sense=minimize, name='ObjectiveFunction') # Objective function11
model.Eq6 = Constraint(model.GB, model.t, rule=eq6) # Eq.6.8d - Ramp up constraint
model.Eq7 = Constraint(model.GB, model.t, rule=eq7) # Eq.6.8e - Ramp down constraint



# Resolvendo o modelo
solver = SolverFactory('ipopt')
#set for more iterations
#solver.options['max_iter'] = 100000
#solver.options['tol'] = 1e-6
#solver.options['print_level'] = 5  # Increase verbosity for debugging


from amplpy import modules


#solver "c:\ampl\conopt.exe"

#solver = SolverFactory('conopt',
#                       executable=r'C:\ampl\conopt.exe')
results = solver.solve(model, tee=True, load_solutions=True)


#solver_name = "conopt"  # "highs", "cbc",  "couenne", "bonmin", "ipopt", "scip", or "gcg".
#solver = SolverFactory(solver_name+"nl", executable=modules.find(solver_name), solve_io="nl")




#convert model.Pij to a DataFrame
#Pij_data = {(i, j, t): model.Pij[i, j, t].value for i in model.i for j in model.i for t in model.t if (i, j) in LN}
#Pij_df = pd.DataFrame.from_dict(Pij_data, orient='index', columns=['Pij']).reset_index()

Pg_data = {(i, t): model.Pg[i, t].value for i in model.GB for t in model.t}
Pg_df = pd.DataFrame.from_dict(Pg_data, orient='index', columns=['Pg']).reset_index()
Pg_df['Time'] = Pg_df['index'].apply(lambda x: x[1])  # Extract time from the index
Pg_df['Generator'] = Pg_df['index'].apply(lambda x: x[0])  # Extract generator from the index
Pg_df = Pg_df.pivot(index='Time', columns='Generator', values='Pg')
Pg_df = Pg_df*100  # Convert to MW

#construct a df for  Pij
Pij_data = {(i, j, t): model.Pij[i, j, t].value for i in model.i for j in model.i for t in model.t if (i, j) in LN}
Pij_df = pd.DataFrame.from_dict(Pij_data, orient='index', columns=['Pij']).reset_index()
Pij_df.Pij = Pij_df.Pij*100 # Convert to MW

# Iterar sobre todas as variáveis do modelo
for var in model.component_objects(Var, active=True):
    # Criar um nome de arquivo baseado no nome da variável
    var_name = var.name
    filename = f"{var_name}.txt"
    
    # Salvar a variável em um arquivo .txt
    with open(filename, 'w') as f:
        f.write(f"Resultados para a variável: {var_name}\n")
        var.pprint(ostream=f)

# Iterar sobre todas as restrições do modelo
for constraint in model.component_objects(Constraint, active=True):
    # Criar um nome de arquivo baseado no nome da restrição
    constraint_name = constraint.name
    filename = f"{constraint_name}.txt"
    
    # Salvar a restrição em um arquivo .txt
    with open(filename, 'w') as f:
        f.write(f"Resultados para a restrição: {constraint_name}\n")
        constraint.pprint(ostream=f)
        
# Relatórios
report = pd.DataFrame(index=model.t, columns=model.i)
for i in model.i:
    for t in model.t:
        report.loc[t, i] = model.V[i, t].value

#model.OF.pprint()

print(model.OF.expr())
#report.to_excel('exemplos/results.xlsx', sheet_name='Results')
#print(report)

#plot Pg_df.plot() with grid on
plt.figure(figsize=(12, 6))
Pg_df.plot(grid=True)
plt.title('Active Power Generation (Pg) by Generator Over Time')
plt.xlabel('Time Period')
plt.ylabel('Active Power (MW)')
plt.legend(title='Generator')
plt.tight_layout()
plt.savefig('exemplos/active_power_generation.png')

a=1
