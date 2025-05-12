import pandas as pd

# Passo 1: Carregar os dados do CSV
csv_file_path = 'data/transmission_line_data.csv'  # Substitua pelo caminho do seu arquivo CSV
data = pd.read_csv(csv_file_path)

# Passo 2: Identificar os nós únicos
nodes = sorted(set(data['From']).union(set(data['To'])))
node_count = len(nodes)

# Criar um dicionário para mapear nós a índices
node_index = {node: idx for idx, node in enumerate(nodes)}

# Passo 3: Construir a matriz de incidência
incidence_matrix = pd.DataFrame(0, index=range(len(data)), columns=nodes)

for i, row in data.iterrows():
    from_node = row['From']
    to_node = row['To']
    
    # Adiciona 1 para o nó de partida e -1 para o nó de chegada
    incidence_matrix.at[i, from_node] = 1
    incidence_matrix.at[i, to_node] = -1

incidence_matrix.to_csv('data/matriz_incidencia_transmissao.csv', index=False)