from pyomo.environ import *
import pandas as pd
import numpy as np

# Create a model
model = ConcreteModel()

# Sets
model.i = RangeSet(1, 24)  # network buses
model.slack = Set(initialize=[13])  # slack bus
model.t = RangeSet(1, 24)  # time periods
model.GB = Set(initialize=[1, 2, 7, 15, 16, 18, 21, 22, 23])  # generating buses

# Parameters
Sbase = 100
VOLL = 10000
VOLW = 50

# Generating units characteristics
GenD = {
    1: (152, 30.4, 13.32, 192, -50, 1.035, 21, 21),
    2: (152, 30.4, 13.32, 192, -50, 1.035, 21, 21),
    7: (350, 75, 20.7, 300, 0, 1.025, 43, 43),
    13: (591, 206.85, 20.93, 591, 0, 1.02, 31, 31),
    15: (215, 66.25, 21, 215, -100, 1.014, 31, 31),
    16: (155, 54.25, 10.52, 155, -50, 1.017, 31, 31),
    18: (400, 100, 5.47, 400, -50, 1.05, 70, 70),
    21: (400, 100, 5.47, 400, -50, 1.05, 70, 70),
    22: (300, 0, 0, 300, -60, 1.05, 53, 53),
    23: (360, 248.5, 10.52, 310, -125, 1.05, 31, 31),
}

# Demands of each bus
BD = {
    1: (108, 22),
    2: (97, 20),
    3: (180, 37),
    4: (74, 15),
    5: (71, 14),
    6: (136, 28),
    7: (125, 25),
    8: (171, 35),
    9: (175, 36),
    10: (195, 40),
    11: (0, 0),
    12: (0, 0),
    13: (265, 54),
    14: (194, 39),
    15: (317, 64),
    16: (100, 20),
    17: (0, 0),
    18: (333, 68),
    19: (181, 37),
    20: (128, 26),
    21: (0, 0),
    22: (0, 0),
    23: (0, 0),
    24: (0, 0),
}

# Network technical characteristics
LN = {
    (1, 2): (0.2, 0.0026, 0.0139, 0.4611, 175),
    (1, 3): (0.0546, 0.2112, 0.0572, 175),
    (1, 5): (0.0218, 0.0845, 0.0229, 175),
    (2, 4): (0.0328, 0.1267, 0.0343, 175),
    (2, 6): (0.0497, 0.192, 0.052, 175),
    (3, 9): (0.0308, 0.119, 0.0322, 175),
    (3, 24): (0.0023, 0.0839, 0, 400),
    (4, 9): (0.0268, 0.1037, 0.0281, 175),
    (5, 10): (0.0228, 0.0883, 0.0239, 175),
    (6, 10): (0.0139, 0.0605, 2.459, 175),
    (7, 8): (0.0159, 0.0614, 0.0166, 175),
    (8, 9): (0.0427, 0.1651, 0.0447, 175),
    (8, 10): (0.0427, 0.1651, 0.0447, 175),
    (9, 11): (0.0023, 0.0839, 0, 400),
    (9, 12): (0.0023, 0.0839, 0, 400),
    (10, 11): (0.0023, 0.0839, 0, 400),
    (10, 12): (0.0023, 0.0839, 0, 400),
    (11, 13): (0.0061, 0.0476, 0.0999, 500),
    (11, 14): (0.0054, 0.0418, 0.0879, 500),
    (12, 13): (0.0061, 0.0476, 0.0999, 500),
    (12, 23): (0.0124, 0.0966, 0.203, 500),
    (13, 23): (0.0111, 0.0865, 0.1818, 500),
    (14, 16): (0.005, 0.0389, 0.0818, 500),
    (15, 16): (0.0022, 0.0173, 0.0364, 500),
    (15, 21): (0.00315, 0.0245, 0.206, 1000),
    (15, 24): (0.0067, 0.0519, 0.1091, 500),
    (16, 17): (0.0033, 0.0259, 0.0545, 500),
    (16, 19): (0.003, 0.0231, 0.0485, 500),
    (17, 18): (0.0018, 0.0144, 0.0303, 500),
    (17, 22): (0.0135, 0.1053, 0.2212, 500),
    (18, 21): (0.00165, 0.01295, 0.109, 1000),
    (19, 20): (0.00255, 0.0198, 0.1666, 1000),
    (20, 23): (0.0014, 0.0108, 0.091, 1000),
    (21, 22): (0.0087, 0.0678, 0.1424, 500),
}

# Wind and demand characteristics
WD = {
    1: (0, 1),
    2: (0.0786666666666667, 0.684511335492475),
    3: (0.0866666666666667, 0.644122690036197),
    4: (0.117333333333333, 0.6130691560297),
    5: (0.258666666666667, 0.599733282530006),
    6: (0.361333333333333, 0.588874071251667),
    7: (0.566666666666667, 0.5980186702229),
    8: (0.650666666666667, 0.626786054486569),
    9: (0.566666666666667, 0.651743189178891),
    10: (0.484, 0.706039245570585),
    11: (0.548, 0.787007048961707),
    12: (0.757333333333333, 0.839016955610593),
    13: (0.710666666666667, 0.852733854067441),
    14: (0.870666666666667, 0.870642027052772),
    15: (0.932, 0.834254143646409),
    16: (0.966666666666667, 0.816536483139646),
    17: (1, 0.819394170318156),
    18: (0.869333333333333, 0.874071251666984),
    19: (0.665333333333333, 1),
    20: (0.656, 0.983615926843208),
    21: (0.561333333333333, 0.936368832158506),
    22: (0.565333333333333, 0.887597637645266),
    23: (0.556, 0.809297008954087),
    24: (0.724, 0.74585635359116),
}

# Wind capacity
Wcap = {8: 200, 19: 150, 21: 100}

# Variable declarations
model.OF = Var(domain=NonNegativeReals)  # Objective function
model.Pij = Var(model.i, model.i, model.t, domain=Reals)  # Line flow
model.Qij = Var(model.i, model.i, model.t, domain=Reals)  # Reactive power flow
model.Pg = Var(model.GB, model.t, bounds=(0, None))  # Active power generation
model.Qg = Var(model.GB, model.t, bounds=(0, None))  # Reactive power generation
model.Va = Var(model.i, model.t, bounds=(-np.pi/2, np.pi/2))  # Voltage angle
model.V = Var(model.i, model.t, bounds=(0.9, 1.1))  # Voltage magnitude
model.Pw = Var(model.GB, model.t, bounds=(0, None))  # Wind power


# Equations
def flow_eq1(model, i, j, t):
    if (i, j) in LN:
        return model.Pij[i, j, t] == (
            model.V[i, t]**2 * cos(LN[i, j][2]) - 
            model.V[i, t] * model.V[j, t] * cos(model.Va[i, t] - model.Va[j, t] + LN[i, j][2])
        ) / (sqrt(LN[i, j][1]**2 + LN[i, j][0]**2))
    
model.flow_eq1 = Constraint(model.i, model.i, model.t, rule=flow_eq1)

def flow_eq2(model, i, j, t):
    if (i, j) in LN:
        return model.Qij[i, j, t] == (
            model.V[i, t]**2 * sin(LN[i, j][2]) - 
            model.V[i, t] * model.V[j, t] * sin(model.Va[i, t] - model.Va[j, t] + LN[i, j][2])
        ) / (sqrt(LN[i, j][1]**2 + LN[i, j][0]**2)) - (LN[i, j][0] * model.V[i, t]**2) / 2

model.flow_eq2 = Constraint(model.i, model.i, model.t, rule=flow_eq2)

def power_balance(model, i, t):
    return sum(model.Pij[i, j, t] for j in model.i if (j, i) in LN) + model.Pg[i, t] - WD[t][1] * BD[i][0] / Sbase == 0

model.power_balance = Constraint(model.i, model.t, rule=power_balance)

def reactive_power_balance(model, i, t):
    return sum(model.Qij[i, j, t] for j in model.i if (j, i) in LN) + model.Qg[i, t] - WD[t][1] * BD[i][1] / Sbase == 0

model.reactive_power_balance = Constraint(model.i, model.t, rule=reactive_power_balance)

def objective_function(model):
    return sum(model.Pg[i, t] * GenD[i][1] * Sbase for i in model.GB for t in model.t)

model.objective = Objective(rule=objective_function, sense=minimize)

# Generation limits
for i in model.GB:
    for t in model.t:
        model.Pg[i, t].setlb(GenD[i][1] / Sbase)  # Pmin
        model.Pg[i, t].setub(GenD[i][0] / Sbase)  # Pmax
        model.Qg[i, t].setlb(GenD[i][4] / Sbase)  # Qmin
        model.Qg[i, t].setub(GenD[i][3] / Sbase)  # Qmax

# Voltage angle limits
for i in model.i:
    for t in model.t:
        model.Va[i, t].setlb(-np.pi/2)
        model.Va[i, t].setub(np.pi/2)
        model.V[i, t].setlb(0.9)
        model.V[i, t].setub(1.1)

# Solve the model
solver = SolverFactory('ipopt')
results = solver.solve(model, tee=True)

# Prepare output report
report = pd.DataFrame(index=model.t, columns=model.i)
for i in model.i:
    for t in model.t:
        report.loc[t, i] = model.V[i, t].value

report.to_excel('results.xlsx', index=True, sheet_name='voltages')