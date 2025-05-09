
import numpy as np
import pandas as pd

import pyomo.environ as pyo

####################### LENDO OS DADOS #######################
files = ['scenario_distribution_connected_loads', 'scenario_distribution_connected_wind_farms', 'scenario_solar_park', 
         'scenario_transmission_connected_loads', 'scenario_transmission_connected_wind_farms',
'distribution_line_data', 'distribution_node_data', 'dso_connected_dispatchable_generators', 'dso_connected_non_dispatchable_generators',
'load_shape_data', 'transmission_connected_dispatchable_generators', 'transmission_connected_non_dispatchable_generators',
'transmission_line_data', 'transmission_node_data','matriz_incidencia_distribuicao', 'matriz_incidencia_transmissao'] 

#leia todos files
for file in files:
    globals()[file] = pd.read_csv('data/' + file + '.csv')


# converte o dataframe em uma lista de listas
matriz_incidencia_transmissao = matriz_incidencia_transmissao.values.tolist()

#matriz_incidencia_transmissao = matriz_incidencia_transmissao
#matriz_incidencia_transmissao.index = matriz_incidencia_transmissao.index + 1
#matriz_incidencia_dict = { (int(idx[0]), int(idx[1])): value for idx, value in matriz_incidencia_transmissao.stack().items()}




################################################################

####################### mO SUPERIOR ########################


m = pyo.ConcreteModel()

# Conjuntos



m.G = pyo.Set(initialize=[i for i in range(1, 11)])  # lances fornecidos por geradores

m.D = pyo.Set(initialize=[i for i in range(1, 35)])  # D é um conjunto de geradores conectados à distribuição

m.I = pyo.Set(initialize=[i for i in range(1, 35)])  # I é um conjunto de geradores conectados ao nó i

m.T = pyo.Set(initialize=[i for i in range(1, 11)])  # T é um conjunto de geradores conectados à transmissão

m.LD = pyo.Set(initialize=[i for i in range(1, 21)])  # LD é um conjunto de linhas de distribuição

m.LT = pyo.Set(initialize=[i for i in range(1, 21)])  # LT é um conjunto de linhas de transmissão

m.ND = pyo.Set(initialize=[i for i in range(1, 35)])  # ND é um conjunto de nós de distribuição

m.NT = pyo.Set(initialize=[i for i in range(1, 15)])  # NT é um conjunto de nós de transmissão

m.N_inf = pyo.Set(initialize=[i for i in range(1, 10)])  # N_inf é um conjunto de nós de interface

m.s = pyo.RangeSet(1, 5)  # S é um conjunto de 5 cenários
m.t = pyo.RangeSet(1, 24)  # T é um conjunto de períodos de tempo


m.j = pyo.RangeSet(1, 20)  # Conjunto de linhas de distribuição

bid_values = np.ones((10, 20))

m.i, m.j = range(0,14), range(0,14) # Conjuntos de linhas de distribuição

m.l = range(0,20)
m.g = range(0,8)
m.b = range(0,5)

m.A = pyo.Param(m.i, m.l, initialize= lambda m, col, lin: matriz_incidencia_transmissao[lin][col]) # Matriz de incidência de linhas de transmissão
m.C_D = pyo.Param(m.g, m.b, initialize = lambda m, g, b: 1)  # Oferta do gerador de distribuição g para o bloco b
m.C_T = pyo.Param(m.g, m.b, initialize = lambda m, g, b: 1)   # Oferta do gerador de transmissão g para o bloco b



# dL_i é a demanda de potência no nó de distribuição i
m.dL = pyo.Param(m.D, m.t, m.s, initialize=demand_values)  # Exemplo de valores

# dT_i é a demanda de potência no nó de transmissão i
m.dT = pyo.Param(m.T, m.t, m.s, initialize=demand_values)  # Exemplo de valores

# p_s,t peso do cenário s no período t
m.p = pyo.Param(m.S, m.T, initialize=weights)  # Exemplo de valores

# P_D_g,b é o tamanho do bloco de oferta de eletricidade do gerador g
m.P_D = pyo.Param(m.G, m.b, initialize=block_size_values)  # Exemplo de valores

# P_T_g,b é o tamanho do bloco de oferta de eletricidade do gerador g na transmissão
m.P_T = pyo.Param(m.G, m.b, initialize=block_size_values)  # Exemplo de valores

# P_ESS_i é a potência máxima de carga/descarga do sistema de armazenamento
m.P_ESS = pyo.Param(m.i, initialize=storage_capacity_values)  # Exemplo de valores

# SOC_i é o estado de carga inicial dos sistemas de armazenamento
m.SOC_initial = pyo.Param(m.i, initialize=initial_soc_values)  # Exemplo de valores

# SOC_max_i é o estado de carga máximo dos sistemas de armazenamento
m.SOC_max = pyo.Param(m.i, initialize=max_soc_values)  # Exemplo de valores

# SOC_min_i é o estado de carga mínimo dos sistemas de armazenamento
m.SOC_min = pyo.Param(m.i, initialize=min_soc_values)  # Exemplo de valores

# U representa a rigidez da regulamentação de emissões de carbono
m.U = pyo.Param(initialize=regulation_flexibility)  # Exemplo

# Função objetivo
def objective_rule(model):
    return sum(m.P_g[g, t, s] * m.C_g_b[g] for g in m.G_D for t in m.T for s in m.S) + \
            sum(m.psi[t, s] for t in m.T for s in m.S) + \
            sum(m.lambda_i[i, t, s] * m.P_i[i, t, s] for i in m.N_D for t in m.T for s in m.S)
m.obj = pyo.Objective(rule=objective_rule, sense=pyo.minimize)



