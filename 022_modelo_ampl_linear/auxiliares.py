
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


def _create_entity_dataframe(ampl: AMPL, entity_type: str) -> pd.DataFrame:
    """
    Função interna genérica para extrair variáveis ou parâmetros de um modelo AMPL
    e consolidá-los em um único DataFrame do Pandas.
    """
    if entity_type == 'variable':
        get_entities_func = ampl.get_variables
        entity_name_str = 'variable'
        print("Iniciando extração e formatação das Variáveis...")
    elif entity_type == 'parameter':
        get_entities_func = ampl.get_parameters
        entity_name_str = 'parameter'
        print("Iniciando extração e formatação dos Parâmetros...")
    else:
        raise ValueError("Tipo de entidade deve ser 'variable' ou 'parameter'")

    list_of_reshaped_dfs = []

    for entity in get_entities_func():
        name = entity[0]
        try:
            # Para parâmetros, usamos get_values(), assim como para variáveis
            df = entity[1].get_values().to_pandas()

            original_index_names = list(df.index.names)
            

            
            # Caso 1: Entidade unidimensional, como Var[t] ou Param[t]
            if df.index.nlevels == 1:
                reshaped = df.T
                reshaped.index = pd.MultiIndex.from_tuples(
                    [(name, 'value')], names=[entity_name_str, 'index']
                )

            # Caso 2: Entidade multidimensional, como Var[i, j, t]
            else:
                reshaped = df.unstack(level=-1)
                
                if isinstance(reshaped.columns, pd.MultiIndex):
                    reshaped.columns = reshaped.columns.get_level_values(1)
                
                reshaped = reshaped.reset_index()
                reshaped.insert(0, entity_name_str, name)
                
                final_index_cols = [entity_name_str] + original_index_names[:-1]
                reshaped = reshaped.set_index(final_index_cols)
            
            list_of_reshaped_dfs.append(reshaped)
            print(f"  - {entity_name_str.capitalize()} '{name}' processada com sucesso.")

        except Exception as e:
            print(f"  - Erro ao processar {entity_name_str} '{name}': {e}")

    if not list_of_reshaped_dfs:
        print(f"Nenhum {entity_name_str} com dados indexados por tempo foi encontrado.")
        return pd.DataFrame()

    final_df = pd.concat(list_of_reshaped_dfs)
    
    final_df.index.names = [n if n is not None else 'index' for n in final_df.index.names]
    final_df.columns.name = "Time Period"

    print(f"DataFrame final de {entity_name_str}s criado com sucesso!")
    return final_df

def get_ampl_dataframes(ampl: AMPL) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extrai todas as variáveis e parâmetros de um modelo AMPL resolvido e os 
    retorna como dois DataFrames separados.

    Args:
        ampl: Uma instância do objeto amplpy.AMPL com um modelo resolvido.

    Returns:
        Uma tupla contendo dois DataFrames do Pandas: (variables_df, parameters_df).
    """
    variables_df = _create_entity_dataframe(ampl, 'variable')
    print("-" * 50)
    parameters_df = _create_entity_dataframe(ampl, 'parameter')
    return variables_df, parameters_df

