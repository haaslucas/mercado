
import pyomo.environ as pyo
import pandas as pd
import os
import glob

####################### LENDO OS DADOS #######################

df = pd.read_excel('inputdata.xlsx')

#read all excel files inside folder data and subfolders
path = os.path.join(os.getcwd(), 'data')
all_files = glob.glob(os.path.join(path, '*.xlsx'))
#name of the dataframe should be the name of the file without the extension
files = []
for file in all_files:
    file_name = os.path.splitext(os.path.basename(file))[0]
    file_name = file_name.lower()
    file_name = file_name.replace('-', '_')
    globals()[file_name] = pd.read_excel(file)
    files.append(file_name)

SBASE = 100 #MVA
VBASE = 24.9 #kV

############## Cenários ##############
df = scenario_distribution_connected_loads
df.columns = ['#Scenario', 'Hour', 'Loading Factor'] * 5 
df = pd.concat([df.iloc[:,:3], df.iloc[:,3:6], df.iloc[:,6:9], df.iloc[:,9:12], df.iloc[:,12:15]], axis=0)
scenario_distribution_connected_loads = df.copy()

df = scenario_distribution_connected_wind_farms
df.columns = ['#Scenario', 'Hour', 'Power Availability (pu)'] * 5
df = pd.concat([df.iloc[:,:3], df.iloc[:,3:6], df.iloc[:,6:9], df.iloc[:,9:12], df.iloc[:,12:15]], axis=0)
scenario_distribution_connected_wind_farms = df.copy()

df = scenario_solar_park
df.columns = ['#Scenario', 'Hour', 'Power Availability (pu)'] * 5
df = pd.concat([df.iloc[:,:3], df.iloc[:,3:6], df.iloc[:,6:9], df.iloc[:,9:12], df.iloc[:,12:15]], axis=0)
scenario_solar_park = df.copy()

df = scenario_transmission_connected_loads
df.columns = ['#Scenario', 'Hour', 'Loading Factor'] * 5
df = pd.concat([df.iloc[:,:3], df.iloc[:,3:6], df.iloc[:,6:9], df.iloc[:,9:12], df.iloc[:,12:15]], axis=0)
scenario_transmission_connected_loads = df.copy()

df = scenario_transmission_connected_wind_farms
df.columns = ['#Scenario', 'Hour', 'Power Availability (pu)'] * 5
df = pd.concat([df.iloc[:,:3], df.iloc[:,3:6], df.iloc[:,6:9], df.iloc[:,9:12], df.iloc[:,12:15]], axis=0)
scenario_transmission_connected_wind_farms = df.copy()


files = ['scenario_distribution_connected_loads', 'scenario_distribution_connected_wind_farms', 'scenario_solar_park', 
         'scenario_transmission_connected_loads', 'scenario_transmission_connected_wind_farms',
'distribution_line_data', 'distribution_node_data', 'dso_connected_dispatchable_generators', 'dso_connected_non_dispatchable_generators',
'load_shape_data', 'transmission_connected_dispatchable_generators', 'transmission_connected_non_dispatchable_generators',
'transmission_line_data', 'transmission_node_data']

# agora salve todas as variáveis em csvs.
for file in files:
    df = globals()[file]
    df.to_csv(file + '.csv', index=False)
    
    
    
def create_superior_model():
    model = pyo.ConcreteModel()

    # Conjuntos
    model.S = pyo.Set(initialize=['s1'])  # Cenários
    model.T = pyo.Set(initialize=['t1', 't2'])  # Períodos de tempo
    model.G_D = pyo.Set(initialize=['G1', 'G2'])  # Geradores no sistema de distribuição
    model.N_D = pyo.Set(initialize=['i1', 'i2'])  # Nodos do sistema de distribuição
    model.L_D = pyo.Set(initialize=[('i1', 'i2'), ('i2', 'i1')])  # Linhas de distribuição

    # Parâmetros
    model.C_g_b = pyo.Param(model.G_D, within=pyo.NonNegativeReals, initialize={'G1': 10, 'G2': 15})  # Custo de geração
    model.R_ij = pyo.Param(model.L_D, within=pyo.NonNegativeReals, initialize={('i1', 'i2'): 0.1, ('i2', 'i1'): 0.1})  # Resistência de linhas
    
    # Variáveis de decisão
    model.P_g = pyo.Var(model.G_D, model.T, model.S, within=pyo.NonNegativeReals)  # Potência gerada
    model.P_i = pyo.Var(model.N_D, model.T, model.S, within=pyo.NonNegativeReals)  # Potência injetada nos nodos
    model.L_i = pyo.Var(model.N_D, model.T, model.S, within=pyo.NonNegativeReals)  # Demanda nos nodos
    model.psi = pyo.Var(model.T, model.S, within=pyo.NonNegativeReals)  # Penalidade
    model.lambda_i = pyo.Var(model.N_D, model.T, model.S, within=pyo.NonNegativeReals)  # Preço sombra

    # Função objetivo
    def objective_rule(model):
        return sum(model.P_g[g, t, s] * model.C_g_b[g] for g in model.G_D for t in model.T for s in model.S) + \
                sum(model.psi[t, s] for t in model.T for s in model.S) + \
                sum(model.lambda_i[i, t, s] * model.P_i[i, t, s] for i in model.N_D for t in model.T for s in model.S)
    model.obj = pyo.Objective(rule=objective_rule, sense=pyo.minimize)

    # Restrições
    def power_balance_rule(model, i, t, s):
        return (sum(model.P_g[g, t, s] for g in model.G_D) +
                sum(model.P_i[k, t, s] for k, _ in model.L_D) -
                sum(model.P_i[j, t, s] for _, j in model.L_D) ==
                model.P_i[i, t, s] + model.L_i[i, t, s] for i in model.N_D)
    model.power_balance_constraint = pyo.Constraint(model.N_D, model.T, model.S, rule=power_balance_rule)

    def power_balance_infinite_limit_rule(model, i, t, s):
        return (sum(model.P_g[g, t, s] for g in model.G_D) +
                sum(model.P_i[k, t, s] for k, _ in model.L_D) -
                sum(model.P_i[j, t, s] for _, j in model.L_D) ==
                model.P_i[i, t, s] + model.L_i[i, t, s] for i in model.N_D)
    model.power_balance_infinite_limit_constraint = pyo.Constraint(model.N_D, model.T, model.S, rule=power_balance_infinite_limit_rule)

    def reactive_power_balance_rule(model, i, t, s):
        return (sum(model.Q_g[g, t, s] for g in model.G_D) +
                sum(model.Q_i[k, t, s] for k, _ in model.L_D) -
                sum(model.Q_i[j, t, s] for _, j in model.L_D) ==
                model.K_i * (model.P_i[i, t, s] + model.L_i[i, t, s]) for i in model.N_D)
    model.reactive_power_balance_constraint = pyo.Constraint(model.N_D, model.T, model.S, rule=reactive_power_balance_rule)

    return model

if __name__ == "__main__":
    model = create_superior_model()
    # Aqui você pode adicionar a resolução do modelo


'''
função pra criar o dist_load_df


# Supondo que Dist_Load seja um dicionário para armazenar os valores
Dist_Load = {}

for n in dnodes['Node'].unique():  # Acesse valores únicos de 'Node'
    for t in dload['Hour'].unique():
        for s in dload_scn['#Scenario'].unique():
            # Obtenha o tipo de nó e a demanda máxima
            node_type = dnodes['Load Shape'][dnodes['Node'] == n].values[0]
            max_demand = dnodes['Maximum Demand (MW)'][dnodes['Node'] == n].values[0]

            if node_type == 'IND1':
                Dist_Load[(n, t, s)] = max_demand * dload['ind1'][dload['Hour'] == t].values[0] * dload_scn['Loading Factor'][dload_scn['#Scenario'] == s].values[0] / SBASE
            elif node_type == 'IND2':
                Dist_Load[(n, t, s)] = max_demand * dload['ind2'][dload['Hour'] == t].values[0] * dload_scn['Loading Factor'][dload_scn['#Scenario'] == s].values[0] / SBASE
            elif node_type == 'IND3':
                Dist_Load[(n, t, s)] = max_demand * dload['ind3'][dload['Hour'] == t].values[0] * dload_scn['Loading Factor'][dload_scn['#Scenario'] == s].values[0] / SBASE
            elif node_type == 'RES1':
                Dist_Load[(n, t, s)] = max_demand * dload['res1'][dload['Hour'] == t].values[0] * dload_scn['Loading Factor'][dload_scn['#Scenario'] == s].values[0] / SBASE
            elif node_type == 'RES4':
                Dist_Load[(n, t, s)] = max_demand * dload['res4'][dload['Hour'] == t].values[0] * dload_scn['Loading Factor'][dload_scn['#Scenario'] == s].values[0] / SBASE
            elif node_type == 'RES5':
                Dist_Load[(n, t, s)] = max_demand * dload['res5'][dload['Hour'] == t].values[0] * dload_scn['Loading Factor'][dload_scn['#Scenario'] == s].values[0] / SBASE
            elif node_type == 'OFF1':
                Dist_Load[(n, t, s)] = max_demand * dload['off1'][dload['Hour'] == t].values[0] * dload_scn['Loading Factor'][dload_scn['#Scenario'] == s].values[0] / SBASE
            elif node_type == 'OFF4':
                Dist_Load[(n, t, s)] = max_demand * dload['off4'][dload['Hour'] == t].values[0] * dload_scn['Loading Factor'][dload_scn['#Scenario'] == s].values[0] / SBASE
                
# Converter o dicionário Dist_Load para um DataFrame
dist_load_df = pd.DataFrame.from_dict(
    Dist_Load,
    orient='index',  # Os índices do dicionário se tornam as linhas do DataFrame
    columns=['Load']  # Coluna para os valores do dicionário
)

# Adiciona colunas adicionais de Node, Hour e Scenario
dist_load_df.reset_index(inplace=True)
dist_load_df[['Node', 'Hour', 'Scenario']] = pd.DataFrame(dist_load_df['index'].tolist(), index=dist_load_df.index)
dist_load_df.drop(columns='index', inplace=True)

# Reorganiza as colunas
dist_load_df = dist_load_df[['Node', 'Hour', 'Scenario', 'Load']]

# sort dist_load_df by Scenario, Node and Hour
dist_load_df.sort_values(by=['Scenario', 'Hour', 'Node'], inplace=True)
dist_load_df.to_csv(pasta + 'dist_load.csv', index=False)

df = pd.DataFrame.from_dict(Dist_Load, orient='index', columns=['Dist_Load'])
df.reset_index(inplace=True)




'''