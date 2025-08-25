
import pandas as pd
from amplpy import AMPL

def create_variables_dataframe(ampl: AMPL) -> pd.DataFrame:
    """
    Extrai todas as variáveis de um modelo AMPL resolvido e as consolida em 
    um único DataFrame do Pandas.

    A função assume que o último índice de cada variável é o período de tempo,
    que será transformado nas colunas do DataFrame resultante.

    Exemplo de transformação:
    - Uma variável AMPL `P[i, t]`
    - Se torna linhas no DataFrame com um MultiIndex `('P', i)`
    - E colunas correspondentes aos valores de `t`.

    Args:
        ampl: Uma instância do objeto amplpy.AMPL com um modelo resolvido.

    Returns:
        Um DataFrame do Pandas contendo os valores de todas as variáveis,
        formatado com o tempo como colunas.
    """
    list_of_reshaped_dfs = []

    print("Iniciando extração e formatação das variáveis...")

    # Itera sobre todas as variáveis definidas no modelo
    for name, var in ampl.get_variables():
        try:
            df = var.get_values().to_pandas()

            # Pula variáveis sem valores
            if df.empty:
                print(f"  - Variável '{name}' está vazia. Pulando.")
                continue

            # Guarda os nomes originais dos índices (ex: ['n', 't'])
            original_index_names = list(df.index.names)
            
            # Caso 1: Variável unidimensional, como Var[t]
            if df.index.nlevels == 1:
                # Transpõe o DataFrame para que o tempo vire colunas
                reshaped = df.T
                # Define o índice da linha como o nome da variável
                reshaped.index = pd.MultiIndex.from_tuples(
                    [(name, 'value')], names=['variable', 'index']
                )

            # Caso 2: Variável multidimensional, como Var[i, j, t]
            else:
                # Usa 'unstack' para pivotar o último nível do índice (tempo) para as colunas
                reshaped = df.unstack(level=-1)
                
                # As colunas agora são um MultiIndex (ex: ('P_DSO', 'T1'), ('P_DSO', 'T2'))
                # Simplificamos para ter apenas o valor do tempo ('T1', 'T2')
                reshaped.columns = reshaped.columns.get_level_values(1)
                
                # Adiciona o nome da variável ao índice da linha
                reshaped = reshaped.reset_index()
                reshaped.insert(0, 'variable', name)
                
                # Define o novo índice hierárquico
                # (ex: ['variable', 'n', 'm'] para uma variável Var[n, m, t])
                final_index_cols = ['variable'] + original_index_names[:-1]
                reshaped = reshaped.set_index(final_index_cols)
            
            list_of_reshaped_dfs.append(reshaped)
            print(f"  - Variável '{name}' processada com sucesso.")

        except Exception as e:
            print(f"  - Erro ao processar a variável '{name}': {e}")

    if not list_of_reshaped_dfs:
        print("Nenhuma variável com dados foi encontrada para criar o DataFrame.")
        return pd.DataFrame()

    # Concatena todos os DataFrames individuais em um só
    final_df = pd.concat(list_of_reshaped_dfs)
    
    # Renomeia os eixos para clareza
    final_df.index.names = [name if name is not None else 'index' for name in final_df.index.names]
    final_df.columns.name = "Time Period"

    print("DataFrame final criado com sucesso!")
    return final_df