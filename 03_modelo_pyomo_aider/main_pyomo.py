import pyomo.environ as pyo
from modelo_pyomo import model as abstract_model # Import the abstract model structure
import data_parser # Import the parser utility
import os

# --- Control Parameters from execute.run ---
# These would typically be command-line arguments or read from a config file
U_param = 1.00
Remove_Energy_Mkt_param = 1
Remove_Carbon_Mkt_param = 1
Remove_Ramp_param = 0

# --- Path to data files ---
# Assuming the script is run from the root of the repository
BASE_DIR = os.path.dirname(os.path.abspath(__file__)) # Or adjust if main_pyomo.py is elsewhere
INPUT_DAT_PATH = os.path.join(BASE_DIR, "input.dat")
SCENARIOS_DAT_PATH = os.path.join(BASE_DIR, "Reducao_Cenarios","Scenarios.dat") # Corrected "Reducao_Cenarios" to "Reducao Cenarios"


def load_data():
    """
    Loads all data from .dat files into Python dictionaries.
    This is a simplified version and will need to be expanded.
    """
    data = {}
    
    # Load from Scenarios.dat
    # Assuming S contains integer-like scenario numbers that might be parsed as strings.
    s_values_from_parser = data_parser.parse_ampl_set(SCENARIOS_DAT_PATH, "S")
    try:
        data['S'] = [int(s) for s in s_values_from_parser]
        # print(f"DEBUG: Set S loaded as integers: {data['S'][:5]}...") # Optional: for debugging
    except ValueError:
        print(f"Warning: Could not convert all elements of set S to integers. Using original values: {s_values_from_parser[:5]}...")
        data['S'] = s_values_from_parser # Fallback to original if conversion fails

    data['wind1'] = data_parser.parse_ampl_simple_param_2d(SCENARIOS_DAT_PATH, "wind1", ['s', 't'], 'val')
    data['wind2'] = data_parser.parse_ampl_simple_param_2d(SCENARIOS_DAT_PATH, "wind2", ['s', 't'], 'val')
    data['load1'] = data_parser.parse_ampl_simple_param_2d(SCENARIOS_DAT_PATH, "load1", ['s', 't'], 'val')
    data['load2'] = data_parser.parse_ampl_simple_param_2d(SCENARIOS_DAT_PATH, "load2", ['s', 't'], 'val')
    data['solar'] = data_parser.parse_ampl_simple_param_2d(SCENARIOS_DAT_PATH, "solar", ['s', 't'], 'val')
    data['probability'] = data_parser.parse_ampl_simple_param_2d(SCENARIOS_DAT_PATH, "probability", ['s', 't'], 'val')

    # Load from input.dat
    # param: T: table
    t_table_cols = ['res1', 'res4', 'res5', 'off1', 'off4', 'ind1', 'ind2', 'ind3', 'L2', 'L3', 'L4', 'L6', 'L8']
    t_data_rows = data_parser.parse_ampl_indexed_table(INPUT_DAT_PATH, "param:	T:", "T_idx", t_table_cols)
    # Assuming T contains integer-like period numbers that might be parsed as strings from table keys.
    t_keys_from_parser = list(t_data_rows.keys())
    try:
        data['T'] = [int(t) for t in t_keys_from_parser]
        # print(f"DEBUG: Set T loaded as integers: {data['T'][:5]}...") # Optional: for debugging
    except ValueError:
        print(f"Warning: Could not convert all elements of set T to integers. Using original values: {t_keys_from_parser[:5]}...")
        data['T'] = t_keys_from_parser # Fallback to original if conversion fails
    # Store individual parameters like model.ind1, model.res1, etc.
    for col_name in t_table_cols:
        data[col_name] = {idx: row_data[col_name] for idx, row_data in t_data_rows.items()}

    # param: Trans_Lines: table
    trans_lines_cols = ['idx','From', 'To', 'Trans_Reactance', 'Trans_Capacity', 'Trans_Status']
    # The parser needs to be adapted or a more general one used.
    # For now, let's assume parse_general_table can fetch rows.
    # The first column in the file is the index of Trans_Lines, not explicitly named 'idx' in AMPL syntax.
    # We'll need to handle this.
    # Let's assume the parser can extract the index and then the named columns.
    
    # Simplified parsing for Trans_Lines (example, needs robust parser)
    # For 'param: Trans_Lines: From To ... :=' the actual index is implicit (1, 2, ...)
    # The file has '1  1  2  0.05917  158.0  1' where '1' is the line index.
    
    # A more robust way is to read the table into a pandas DataFrame first if possible,
    # or write a more specific parser for this structure.
    # For now, this part is a placeholder for actual parsing of Trans_Lines table.
    data['Trans_Lines_data'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: Trans_Lines:", trans_lines_cols)
    
    # Extract sets and specific params from Trans_Lines_data
    data['Trans_Lines'] = [int(row['idx']) for row in data['Trans_Lines_data']]
    data['From'] = {int(row['idx']): int(row['From']) for row in data['Trans_Lines_data']}
    data['To'] = {int(row['idx']): int(row['To']) for row in data['Trans_Lines_data']}
    data['Trans_Reactance'] = {int(row['idx']): float(row['Trans_Reactance']) for row in data['Trans_Lines_data']}
    data['Trans_Capacity_initial'] = {int(row['idx']): float(row['Trans_Capacity']) for row in data['Trans_Lines_data']} # Store initial before scaling
    data['Trans_Status'] = {int(row['idx']): int(row['Trans_Status']) for row in data['Trans_Lines_data']}

    # param: Trans_Nodes: table
    trans_nodes_cols = ['idx', 'Trans_Shift_Min', 'Trans_Shift_Max', 'Trans_P_ESS_min', 'Trans_P_ESS_max', 'Trans_SOC_min', 'Trans_SOC_max', 'Trans_SOC_initial', 'Trans_Type', 'Trans_Inst_Cap']
    data['Trans_Nodes_data'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: Trans_Nodes:", trans_nodes_cols)
    data['Trans_Nodes'] = [int(row['idx']) for row in data['Trans_Nodes_data']]
    # ... populate other Trans_Nodes related params ...
    data['Trans_Type'] = {int(row['idx']): row['Trans_Type'] for row in data['Trans_Nodes_data']}
    data['Trans_Inst_Cap'] = {int(row['idx']): float(row['Trans_Inst_Cap']) for row in data['Trans_Nodes_data']}


    # param: N: table
    n_cols = ['idx', 'Vmin', 'Vmax', 'Vnom', 'Dist_Type', 'Dist_Inst_Cap']
    data['N_data'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: N:", n_cols)
    data['N'] = [int(row['idx']) for row in data['N_data']]
    data['Vmax'] = {int(row['idx']): float(row['Vmax']) for row in data['N_data']}
    data['Vmin'] = {int(row['idx']): float(row['Vmin']) for row in data['N_data']}
    # ... populate other N related params ...
    data['Dist_Type'] = {int(row['idx']): row['Dist_Type'] for row in data['N_data']}
    data['Dist_Inst_Cap'] = {int(row['idx']): float(row['Dist_Inst_Cap']) for row in data['N_data']}


    # param: L: table (resistance, reactance, etc.)
    # This table has (idx1, idx2) as a compound key.
    # param: L: R X Dist_Status Imax :=
    # 1 2 0.0086 0.008 1 850
    l_cols = ['idx1', 'idx2', 'R', 'X', 'Dist_Status', 'Imax']
    data['L_data_rows'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: L:", l_cols)
    
    # Initialize L set and related parameters
    data['L'] = []
    data['R'] = {}
    data['X'] = {}
    data['Dist_Status_val'] = {} # temp name to avoid conflict with model.Dist_Status
    data['Imax_val'] = {}      # temp name
    
    for row in data['L_data_rows']:
        key = (int(row['idx1']), int(row['idx2']))
        data['L'].append(key)
        data['R'][key] = float(row['R'])
        data['X'][key] = float(row['X'])
        data['Dist_Status_val'][key] = int(row['Dist_Status'])
        data['Imax_val'][key] = float(row['Imax'])
        
    # ... (Load other sets and params: G_T, G_D, LS, ESS, RES_T, RES_D, LDA, etc.)
    # This is a substantial amount of parsing logic.
    # For brevity, I'll stop here for the data loading demonstration.
    # A complete implementation would parse all tables from input.dat.

    # param: G_T: table
    gt_cols = ['g_t_idx', 'G_Node', 'a_t', 'b_t', 'c_t', 'Pmax_t', 'Pmin_t', 'Carbon_Cost_t', 'G_Owner_t']
    data['G_T_data_rows'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: G_T:", gt_cols)
    data['G_T'] = [int(row['g_t_idx']) for row in data['G_T_data_rows']]
    data['G_Node'] = {int(row['g_t_idx']): int(row['G_Node']) for row in data['G_T_data_rows']}
    data['a_t_initial'] = {int(row['g_t_idx']): float(row['a_t']) for row in data['G_T_data_rows']}
    data['b_t_initial'] = {int(row['g_t_idx']): float(row['b_t']) for row in data['G_T_data_rows']}
    data['c_t'] = {int(row['g_t_idx']): float(row['c_t']) for row in data['G_T_data_rows']}
    data['Pmax_t_initial'] = {int(row['g_t_idx']): float(row['Pmax_t']) for row in data['G_T_data_rows']}
    data['Pmin_t_initial'] = {int(row['g_t_idx']): float(row['Pmin_t']) for row in data['G_T_data_rows']}
    data['Carbon_Cost_t_initial'] = {int(row['g_t_idx']): float(row['Carbon_Cost_t']) for row in data['G_T_data_rows']}
    data['G_Owner_t'] = {int(row['g_t_idx']): row['G_Owner_t'] for row in data['G_T_data_rows']} # Assuming G_Owner_t can be non-integer

    # param: G_D: table
    gd_cols = ['g_d_idx', 'a_d', 'b_d', 'c_d', 'Pmax_d', 'Pmin_d', 'G_LDA', 'Carbon_Cost_d', 'G_LDA_Node', 'G_Owner_d']
    data['G_D_data_rows'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: G_D:", gd_cols)
    data['G_D'] = [int(row['g_d_idx']) for row in data['G_D_data_rows']]
    data['a_d_initial'] = {int(row['g_d_idx']): float(row['a_d']) for row in data['G_D_data_rows']}
    data['b_d_initial'] = {int(row['g_d_idx']): float(row['b_d']) for row in data['G_D_data_rows']}
    data['c_d'] = {int(row['g_d_idx']): float(row['c_d']) for row in data['G_D_data_rows']}
    data['Pmax_d_initial'] = {int(row['g_d_idx']): float(row['Pmax_d']) for row in data['G_D_data_rows']}
    data['Pmin_d_initial'] = {int(row['g_d_idx']): float(row['Pmin_d']) for row in data['G_D_data_rows']}
    data['G_LDA_val'] = {int(row['g_d_idx']): int(row['G_LDA']) for row in data['G_D_data_rows']} # temp name
    data['Carbon_Cost_d_initial'] = {int(row['g_d_idx']): float(row['Carbon_Cost_d']) for row in data['G_D_data_rows']}
    data['G_LDA_Node_val'] = {int(row['g_d_idx']): int(row['G_LDA_Node']) for row in data['G_D_data_rows']} # temp name
    data['G_Owner_d'] = {int(row['g_d_idx']): row['G_Owner_d'] for row in data['G_D_data_rows']}

    # param: LS: table
    ls_cols = ['ls_idx', 'LS_Node', 'Dist_Shift_Min', 'Dist_Shift_Max']
    data['LS_data_rows'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: LS:", ls_cols)
    data['LS'] = [int(row['ls_idx']) for row in data['LS_data_rows']]
    data['LS_Node_val'] = {int(row['ls_idx']): int(row['LS_Node']) for row in data['LS_data_rows']} # temp name
    data['Dist_Shift_Min_val'] = {int(row['ls_idx']): float(row['Dist_Shift_Min']) for row in data['LS_data_rows']} # temp name
    data['Dist_Shift_Max_val'] = {int(row['ls_idx']): float(row['Dist_Shift_Max']) for row in data['LS_data_rows']} # temp name

    # param: ESS: table
    ess_cols = ['ess_idx', 'ESS_Node', 'Dist_P_ESS_min', 'Dist_P_ESS_max', 'Dist_SOC_min', 'Dist_SOC_max', 'Dist_SOC_initial']
    data['ESS_data_rows'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: ESS:", ess_cols)
    data['ESS'] = [int(row['ess_idx']) for row in data['ESS_data_rows']]
    data['ESS_Node_val'] = {int(row['ess_idx']): int(row['ESS_Node']) for row in data['ESS_data_rows']} # temp name
    data['Dist_P_ESS_min_initial'] = {int(row['ess_idx']): float(row['Dist_P_ESS_min']) for row in data['ESS_data_rows']}
    data['Dist_P_ESS_max_initial'] = {int(row['ess_idx']): float(row['Dist_P_ESS_max']) for row in data['ESS_data_rows']}
    data['Dist_SOC_min_initial'] = {int(row['ess_idx']): float(row['Dist_SOC_min']) for row in data['ESS_data_rows']}
    data['Dist_SOC_max_initial'] = {int(row['ess_idx']): float(row['Dist_SOC_max']) for row in data['ESS_data_rows']}
    data['Dist_SOC_initial_val'] = {int(row['ess_idx']): float(row['Dist_SOC_initial']) for row in data['ESS_data_rows']} # temp name

    # param: RES_T: table
    res_t_cols = ['res_t_idx', 'RES_type_t', 'RES_Node_t', 'RES_Cap_t', 'Wind_min_t', 'Wind_max_t', 'Wind_nom_t']
    data['RES_T_data_rows'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: RES_T:", res_t_cols)
    data['RES_T'] = [int(row['res_t_idx']) for row in data['RES_T_data_rows']]
    data['RES_type_t_val'] = {int(row['res_t_idx']): row['RES_type_t'] for row in data['RES_T_data_rows']}
    data['RES_Node_t_val'] = {int(row['res_t_idx']): int(row['RES_Node_t']) for row in data['RES_T_data_rows']}
    data['RES_Cap_t_val'] = {int(row['res_t_idx']): float(row['RES_Cap_t']) for row in data['RES_T_data_rows']}
    # Wind_min/max/nom_t can have '.' which means None or needs specific handling.
    # Assuming parse_general_table returns them as strings; convert to float or None.
    def safe_float_or_none(val_str):
        if val_str == '.': return None
        try: return float(val_str)
        except ValueError: return val_str # Or raise error

    data['Wind_min_t_val'] = {int(row['res_t_idx']): safe_float_or_none(row['Wind_min_t']) for row in data['RES_T_data_rows']}
    data['Wind_max_t_val'] = {int(row['res_t_idx']): safe_float_or_none(row['Wind_max_t']) for row in data['RES_T_data_rows']}
    data['Wind_nom_t_val'] = {int(row['res_t_idx']): safe_float_or_none(row['Wind_nom_t']) for row in data['RES_T_data_rows']}

    # param: RES_D: table
    res_d_cols = ['res_d_idx', 'RES_type_d', 'RES_Node_d', 'RES_Cap_d', 'Wind_min_d', 'Wind_max_d', 'Wind_nom_d']
    data['RES_D_data_rows'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: RES_D:", res_d_cols)
    data['RES_D'] = [int(row['res_d_idx']) for row in data['RES_D_data_rows']]
    data['RES_type_d_val'] = {int(row['res_d_idx']): row['RES_type_d'] for row in data['RES_D_data_rows']}
    data['RES_Node_d_val'] = {int(row['res_d_idx']): int(row['RES_Node_d']) for row in data['RES_D_data_rows']}
    data['RES_Cap_d_val'] = {int(row['res_d_idx']): float(row['RES_Cap_d']) for row in data['RES_D_data_rows']}
    data['Wind_min_d_val'] = {int(row['res_d_idx']): safe_float_or_none(row['Wind_min_d']) for row in data['RES_D_data_rows']}
    data['Wind_max_d_val'] = {int(row['res_d_idx']): safe_float_or_none(row['Wind_max_d']) for row in data['RES_D_data_rows']}
    data['Wind_nom_d_val'] = {int(row['res_d_idx']): safe_float_or_none(row['Wind_nom_d']) for row in data['RES_D_data_rows']}

    # param: LDA: PCC:=
    lda_cols = ['lda_idx', 'PCC']
    data['LDA_data_rows'] = data_parser.parse_general_table(INPUT_DAT_PATH, "param: LDA:", lda_cols)
    data['LDA'] = [int(row['lda_idx']) for row in data['LDA_data_rows'] if row['lda_idx'] != '0'] # Exclude '0 0' if it's a comment or placeholder
    data['PCC_val'] = {int(row['lda_idx']): int(row['PCC']) for row in data['LDA_data_rows'] if row['lda_idx'] != '0'}

    # Set WS (Wholesale market participants) - complex set, deferring full parsing if not immediately needed for params
    # For now, an empty set or placeholder if its structure is too complex for current parsers
    data['WS'] = [] # Placeholder

    return data

def create_model_instance(data_dict):
    """
    Creates a concrete instance of the Pyomo model and populates it with data.
    """
    # Create a new ConcreteModel instance from the abstract_model structure
    # This isn't strictly necessary if abstract_model is already Concrete,
    # but good practice if abstract_model were truly abstract.
    # For now, we are directly modifying the imported 'abstract_model' which is already Concrete.
    instance = abstract_model 

    # Populate Sets by clearing and adding new elements
    # Using a list of tuples (set_component, data_key) to avoid hashing issues with Set objects as dict keys.
    set_initialization_list = [
        (instance.S, 'S'),
        (instance.T, 'T'), # ordered=True is already set in model definition
        (instance.Trans_Lines, 'Trans_Lines'),
        (instance.Trans_Nodes, 'Trans_Nodes'),
        (instance.N, 'N'),
        (instance.L, 'L'), # dimen=2 is part of model definition
        (instance.G_T, 'G_T'),
        (instance.G_D, 'G_D'),
        (instance.LS, 'LS'),
        (instance.ESS, 'ESS'),
        (instance.RES_T, 'RES_T'),
        (instance.RES_D, 'RES_D'),
        (instance.LDA, 'LDA'),
        # (instance.WS, 'WS'), # dimen=3 is part of model definition
    ]

    for set_component, data_key in set_initialization_list:
        set_data = data_dict.get(data_key, [])
        set_component.clear()
        if set_data is not None: # data_dict.get with default [] makes this check somewhat redundant for None
            for item in set_data:
                set_component.add(item)
        # For ordered sets, ensure the order is preserved if necessary,
        # though adding in order usually suffices.
        # If T was not ordered, we might need: instance.T.ordered = pyo.Set.SortedOrder

    # Populate Params
    # Scalar params like VBASE, SBASE, IBASE, SE_Capacity are already initialized in modelo_pyomo.py
    
    # Indexed Params that are directly assigned
    param_map_simple_2d = {
        'load1': instance.load1, 'load2': instance.load2,
        'wind1': instance.wind1, 'wind2': instance.wind2,
        'solar': instance.solar, 'probability': instance.probability,
    }
    for name, param_obj in param_map_simple_2d.items():
        for k, v in data_dict.get(name, {}).items():
            param_obj[k] = v
            
    # Params from T_table (ind1, res1, etc.)
    t_table_cols_in_model = { # map file column to model param
        'ind1': instance.ind1, 'ind2': instance.ind2, 'ind3': instance.ind3,
        'res1': instance.res1, 'res4': instance.res4, 'res5': instance.res5,
        'off1': instance.off1, 'off4': instance.off4,
        'L2': instance.L2, 'L3': instance.L3, 'L4': instance.L4, 'L6': instance.L6, 'L8': instance.L8,
    }
    for col_name, param_obj in t_table_cols_in_model.items():
        for k,v in data_dict.get(col_name, {}).items():
            param_obj[k] = v

    # Trans_Lines related params
    if 'From' in data_dict: instance.From.store_values(data_dict['From'])
    if 'To' in data_dict: instance.To.store_values(data_dict['To'])
    if 'Trans_Reactance' in data_dict: instance.Trans_Reactance.store_values(data_dict['Trans_Reactance'])
    if 'Trans_Capacity_initial' in data_dict: instance.Trans_Capacity.store_values(data_dict['Trans_Capacity_initial'])
    if 'Trans_Status' in data_dict: instance.Trans_Status.store_values(data_dict['Trans_Status'])

    # Trans_Nodes related params
    if 'Trans_Type' in data_dict: instance.Trans_Type.store_values(data_dict['Trans_Type'])
    if 'Trans_Inst_Cap' in data_dict: instance.Trans_Inst_Cap.store_values(data_dict['Trans_Inst_Cap'])
    
    # N related params
    if 'Vmax' in data_dict: instance.Vmax.store_values(data_dict['Vmax'])
    if 'Vmin' in data_dict: instance.Vmin.store_values(data_dict['Vmin'])
    if 'Dist_Type' in data_dict: instance.Dist_Type.store_values(data_dict['Dist_Type'])
    if 'Dist_Inst_Cap' in data_dict: instance.Dist_Inst_Cap.store_values(data_dict['Dist_Inst_Cap'])

    # L related params
    if 'R' in data_dict: instance.R.store_values(data_dict['R'])
    if 'X' in data_dict: instance.X.store_values(data_dict['X'])
    if 'Dist_Status_val' in data_dict: instance.Dist_Status.store_values(data_dict['Dist_Status_val'])
    if 'Imax_val' in data_dict: instance.Imax.store_values(data_dict['Imax_val'])
    
    # G_T related params
    if 'G_Node' in data_dict: instance.G_Node.store_values(data_dict['G_Node'])
    if 'a_t_initial' in data_dict: instance.a_t.store_values(data_dict['a_t_initial'])
    if 'b_t_initial' in data_dict: instance.b_t.store_values(data_dict['b_t_initial'])
    if 'c_t' in data_dict: instance.c_t.store_values(data_dict['c_t'])
    if 'Pmax_t_initial' in data_dict: instance.Pmax_t.store_values(data_dict['Pmax_t_initial'])
    if 'Pmin_t_initial' in data_dict: instance.Pmin_t.store_values(data_dict['Pmin_t_initial'])
    if 'Carbon_Cost_t_initial' in data_dict: instance.Carbon_Cost_t.store_values(data_dict['Carbon_Cost_t_initial'])
    if 'G_Owner_t' in data_dict: instance.G_Owner_t.store_values(data_dict['G_Owner_t'])

    # G_D related params
    if 'a_d_initial' in data_dict: instance.a_d.store_values(data_dict['a_d_initial'])
    if 'b_d_initial' in data_dict: instance.b_d.store_values(data_dict['b_d_initial'])
    if 'c_d' in data_dict: instance.c_d.store_values(data_dict['c_d'])
    if 'Pmax_d_initial' in data_dict: instance.Pmax_d.store_values(data_dict['Pmax_d_initial'])
    if 'Pmin_d_initial' in data_dict: instance.Pmin_d.store_values(data_dict['Pmin_d_initial'])
    if 'G_LDA_val' in data_dict: instance.G_LDA.store_values(data_dict['G_LDA_val'])
    if 'Carbon_Cost_d_initial' in data_dict: instance.Carbon_Cost_d.store_values(data_dict['Carbon_Cost_d_initial'])
    if 'G_LDA_Node_val' in data_dict: instance.G_LDA_Node.store_values(data_dict['G_LDA_Node_val'])
    if 'G_Owner_d' in data_dict: instance.G_Owner_d.store_values(data_dict['G_Owner_d'])

    # LS related params
    if 'LS_Node_val' in data_dict: instance.LS_Node.store_values(data_dict['LS_Node_val'])
    if 'Dist_Shift_Min_val' in data_dict: instance.Dist_Shift_Min.store_values(data_dict['Dist_Shift_Min_val'])
    if 'Dist_Shift_Max_val' in data_dict: instance.Dist_Shift_Max.store_values(data_dict['Dist_Shift_Max_val'])

    # ESS related params
    if 'ESS_Node_val' in data_dict: instance.ESS_Node.store_values(data_dict['ESS_Node_val'])
    if 'Dist_P_ESS_min_initial' in data_dict: instance.Dist_P_ESS_min.store_values(data_dict['Dist_P_ESS_min_initial'])
    if 'Dist_P_ESS_max_initial' in data_dict: instance.Dist_P_ESS_max.store_values(data_dict['Dist_P_ESS_max_initial'])
    if 'Dist_SOC_min_initial' in data_dict: instance.Dist_SOC_min.store_values(data_dict['Dist_SOC_min_initial'])
    if 'Dist_SOC_max_initial' in data_dict: instance.Dist_SOC_max.store_values(data_dict['Dist_SOC_max_initial'])
    if 'Dist_SOC_initial_val' in data_dict: instance.Dist_SOC_initial.store_values(data_dict['Dist_SOC_initial_val'])

    # RES_T related params
    if 'RES_type_t_val' in data_dict: instance.RES_type_t.store_values(data_dict['RES_type_t_val'])
    if 'RES_Node_t_val' in data_dict: instance.RES_Node_t.store_values(data_dict['RES_Node_t_val'])
    if 'RES_Cap_t_val' in data_dict: instance.RES_Cap_t.store_values(data_dict['RES_Cap_t_val'])
    if 'Wind_min_t_val' in data_dict: instance.Wind_min_t.store_values({k:v for k,v in data_dict['Wind_min_t_val'].items() if v is not None})
    if 'Wind_max_t_val' in data_dict: instance.Wind_max_t.store_values({k:v for k,v in data_dict['Wind_max_t_val'].items() if v is not None})
    if 'Wind_nom_t_val' in data_dict: instance.Wind_nom_t.store_values({k:v for k,v in data_dict['Wind_nom_t_val'].items() if v is not None})
    
    # RES_D related params
    if 'RES_type_d_val' in data_dict: instance.RES_type_d.store_values(data_dict['RES_type_d_val'])
    if 'RES_Node_d_val' in data_dict: instance.RES_Node_d.store_values(data_dict['RES_Node_d_val'])
    if 'RES_Cap_d_val' in data_dict: instance.RES_Cap_d.store_values(data_dict['RES_Cap_d_val'])
    if 'Wind_min_d_val' in data_dict: instance.Wind_min_d.store_values({k:v for k,v in data_dict['Wind_min_d_val'].items() if v is not None})
    if 'Wind_max_d_val' in data_dict: instance.Wind_max_d.store_values({k:v for k,v in data_dict['Wind_max_d_val'].items() if v is not None})
    if 'Wind_nom_d_val' in data_dict: instance.Wind_nom_d.store_values({k:v for k,v in data_dict['Wind_nom_d_val'].items() if v is not None})

    # LDA related params
    if 'PCC_val' in data_dict: instance.PCC.store_values(data_dict['PCC_val'])

    return instance

def apply_execute_run_logic(model_instance, data_dict):
    """
    Applies the logic from execute.run to modify the model instance.
    """
    # Scale Trans_Capacity
    for l in model_instance.Trans_Lines:
        if model_instance.Trans_Capacity[l] is not None and model_instance.SBASE.value != 0:
            model_instance.Trans_Capacity[l] = model_instance.Trans_Capacity[l].value / model_instance.SBASE.value
        else:
            # Handle cases where initial capacity might be missing or SBASE is zero
            pass 
            
    # Initialize Trans_Incidencia
    for n_node in model_instance.Trans_Nodes:
        for l_line in model_instance.Trans_Lines:
            model_instance.Trans_Incidencia[n_node, l_line] = 0
            
    for l_line in model_instance.Trans_Lines:
        # From and To are 1-indexed in AMPL data, ensure consistency
        from_node = model_instance.From[l_line] # This should be the node ID
        to_node = model_instance.To[l_line]     # This should be the node ID
        
        if from_node in model_instance.Trans_Nodes:
             model_instance.Trans_Incidencia[from_node, l_line] = 1
        if to_node in model_instance.Trans_Nodes:
             model_instance.Trans_Incidencia[to_node, l_line] = -1
             
    # Calculate Trans_Load
    # for{n in Trans_Nodes, t in T, s in S} {	
	# 	if (Trans_Type[n] == "L2") then
	# 		let Trans_Load[n,t,s] := Trans_Inst_Cap[n] * L2[t] * load2[s,t] / SBASE; ...
    # }
    # This requires L2, L3 etc. params and load2 param to be loaded.
    # And Trans_Load param to be mutable.
    if hasattr(model_instance, 'Trans_Load') and model_instance.Trans_Load.is_indexed():
        for n_node in model_instance.Trans_Nodes:
            for t_period in model_instance.T:
                for s_scenario in model_instance.S:
                    load_val = 0
                    trans_type = model_instance.Trans_Type[n_node]
                    inst_cap = model_instance.Trans_Inst_Cap[n_node]
                    sbase_val = model_instance.SBASE.value
                    
                    profile_factor = 0
                    if trans_type == "L2" and t_period in model_instance.L2: profile_factor = model_instance.L2[t_period]
                    elif trans_type == "L3" and t_period in model_instance.L3: profile_factor = model_instance.L3[t_period]
                    elif trans_type == "L4" and t_period in model_instance.L4: profile_factor = model_instance.L4[t_period]
                    elif trans_type == "L6" and t_period in model_instance.L6: profile_factor = model_instance.L6[t_period]
                    elif trans_type == "L8" and t_period in model_instance.L8: profile_factor = model_instance.L8[t_period]
                    # Add other types if necessary
                    
                    load2_factor = model_instance.load2.get((s_scenario, t_period), 1.0) # Default to 1 if missing

                    if sbase_val != 0 and inst_cap is not None:
                        load_val = inst_cap * profile_factor * load2_factor / sbase_val
                    model_instance.Trans_Load[n_node, t_period, s_scenario] = load_val


    # Calculate Dist_Load
    # for{n in N, t in T, s in S} {  	
	# 	if (Dist_Type[n] == "IND1") then
	# 		let Dist_Load[n,t,s] := Dist_Inst_Cap[n] * ind1[t] * load1[s,t] / SBASE; ...
    # }
    if hasattr(model_instance, 'Dist_Load') and model_instance.Dist_Load.is_indexed():
        for n_node in model_instance.N:
            for t_period in model_instance.T:
                for s_scenario in model_instance.S:
                    load_val = 0
                    dist_type = model_instance.Dist_Type[n_node]
                    inst_cap = model_instance.Dist_Inst_Cap[n_node]
                    sbase_val = model_instance.SBASE.value

                    profile_factor = 0
                    if dist_type == "IND1" and t_period in model_instance.ind1: profile_factor = model_instance.ind1[t_period]
                    elif dist_type == "IND2" and t_period in model_instance.ind2: profile_factor = model_instance.ind2[t_period]
                    elif dist_type == "IND3" and t_period in model_instance.ind3: profile_factor = model_instance.ind3[t_period]
                    elif dist_type == "RES1" and t_period in model_instance.res1: profile_factor = model_instance.res1[t_period]
                    elif dist_type == "RES4" and t_period in model_instance.res4: profile_factor = model_instance.res4[t_period]
                    elif dist_type == "RES5" and t_period in model_instance.res5: profile_factor = model_instance.res5[t_period]
                    elif dist_type == "OFF1" and t_period in model_instance.off1: profile_factor = model_instance.off1[t_period]
                    elif dist_type == "OFF4" and t_period in model_instance.off4: profile_factor = model_instance.off4[t_period]
                    
                    load1_factor = model_instance.load1.get((s_scenario, t_period), 1.0)

                    if sbase_val != 0 and inst_cap is not None:
                         load_val = inst_cap * profile_factor * load1_factor / sbase_val
                    model_instance.Dist_Load[n_node, t_period, s_scenario] = load_val
                    
    # Calculate Z and scale Imax
    # let {(n,m) in L} Z[n,m]  := sqrt(R[n,m]^2+X[n,m]^2);
    # let {(n,m) in L} Imax[n,m] := (Imax[n,m]/(SBASE*1e3/VBASE)); (IBASE already defined)
    for line_key in model_instance.L:
        if model_instance.R[line_key] is not None and model_instance.X[line_key] is not None:
            model_instance.Z[line_key] = (model_instance.R[line_key].value**2 + model_instance.X[line_key].value**2)**0.5
        if model_instance.Imax[line_key] is not None and model_instance.Imax[line_key].value is not None and model_instance.IBASE.value != 0:
             model_instance.Imax[line_key] = model_instance.Imax[line_key].value / model_instance.IBASE.value
             
    # Scale ESS parameters
    sbase_val = model_instance.SBASE.value
    if sbase_val != 0:
        for n_ess in model_instance.ESS:
            if model_instance.Dist_P_ESS_min[n_ess].value is not None:
                model_instance.Dist_P_ESS_min[n_ess] = model_instance.Dist_P_ESS_min[n_ess].value / sbase_val
            if model_instance.Dist_P_ESS_max[n_ess].value is not None:
                model_instance.Dist_P_ESS_max[n_ess] = model_instance.Dist_P_ESS_max[n_ess].value / sbase_val
            if model_instance.Dist_SOC_min[n_ess].value is not None:
                model_instance.Dist_SOC_min[n_ess] = model_instance.Dist_SOC_min[n_ess].value / sbase_val
            if model_instance.Dist_SOC_max[n_ess].value is not None:
                model_instance.Dist_SOC_max[n_ess] = model_instance.Dist_SOC_max[n_ess].value / sbase_val
            if model_instance.Dist_SOC_initial[n_ess].value is not None:
                model_instance.Dist_SOC_initial[n_ess] = model_instance.Dist_SOC_initial[n_ess].value / sbase_val

    # Scale Generator parameters (G_T)
    if sbase_val != 0:
        for g in model_instance.G_T:
            if model_instance.Pmax_t[g].value is not None:
                model_instance.Pmax_t[g] = model_instance.Pmax_t[g].value / sbase_val
            if model_instance.Pmin_t[g].value is not None:
                model_instance.Pmin_t[g] = model_instance.Pmin_t[g].value / sbase_val
            if model_instance.a_t[g].value is not None:
                model_instance.a_t[g] = model_instance.a_t[g].value * (sbase_val**2)
            if model_instance.b_t[g].value is not None:
                model_instance.b_t[g] = model_instance.b_t[g].value * sbase_val
            if model_instance.Carbon_Cost_t[g].value is not None:
                model_instance.Carbon_Cost_t[g] = model_instance.Carbon_Cost_t[g].value * sbase_val
                
    # Scale Generator parameters (G_D)
    if sbase_val != 0:
        for g in model_instance.G_D:
            if model_instance.Pmax_d[g].value is not None:
                model_instance.Pmax_d[g] = model_instance.Pmax_d[g].value / sbase_val
            if model_instance.Pmin_d[g].value is not None:
                model_instance.Pmin_d[g] = model_instance.Pmin_d[g].value / sbase_val
            if model_instance.a_d[g].value is not None:
                model_instance.a_d[g] = model_instance.a_d[g].value * (sbase_val**2)
            if model_instance.b_d[g].value is not None:
                model_instance.b_d[g] = model_instance.b_d[g].value * sbase_val
            if model_instance.Carbon_Cost_d[g].value is not None:
                model_instance.Carbon_Cost_d[g] = model_instance.Carbon_Cost_d[g].value * sbase_val

    # Calculate Emission_Weighted_Average
    # Pmax_t and Carbon_Cost_t are now scaled.
    sum_pmax_carbon_cost_t = sum(model_instance.Pmax_t[g].value * model_instance.Carbon_Cost_t[g].value 
                                 for g in model_instance.G_T 
                                 if model_instance.Pmax_t[g].value is not None and model_instance.Carbon_Cost_t[g].value is not None)
    sum_pmax_t = sum(model_instance.Pmax_t[j].value for j in model_instance.G_T if model_instance.Pmax_t[j].value is not None)
    
    if sum_pmax_t != 0:
        ewa = (sum_pmax_carbon_cost_t / sum_pmax_t) * U_param
    else:
        ewa = 0 # Avoid division by zero
    model_instance.Emission_Weighted_Average = ewa

    # Calculate Carbon_Limit_Trans and Carbon_Limit_Dist
    for t_period in model_instance.T:
        # Carbon_Limit_Trans
        sum_max_trans_load_s = 0
        for n_node in model_instance.Trans_Nodes:
            max_load_s = 0
            if any(model_instance.Trans_Load[n_node, t_period, s_scenario].value is not None for s_scenario in model_instance.S):
                 max_load_s = max(model_instance.Trans_Load[n_node, t_period, s_scenario].value 
                                 for s_scenario in model_instance.S 
                                 if model_instance.Trans_Load[n_node, t_period, s_scenario].value is not None)
            sum_max_trans_load_s += max_load_s
        model_instance.Carbon_Limit_Trans[t_period] = model_instance.Emission_Weighted_Average.value * sum_max_trans_load_s

        # Carbon_Limit_Dist
        sum_max_dist_load_s = 0
        for n_node in model_instance.N:
            max_load_s = 0
            if any(model_instance.Dist_Load[n_node, t_period, s_scenario].value is not None for s_scenario in model_instance.S):
                max_load_s = max(model_instance.Dist_Load[n_node, t_period, s_scenario].value 
                                for s_scenario in model_instance.S
                                if model_instance.Dist_Load[n_node, t_period, s_scenario].value is not None)
            sum_max_dist_load_s += max_load_s
        model_instance.Carbon_Limit_Dist[t_period] = model_instance.Emission_Weighted_Average.value * sum_max_dist_load_s
        
    # Conditional logic
    if Remove_Carbon_Mkt_param == 1:
        for t_period in model_instance.T:
            for s_scenario in model_instance.S:
                model_instance.Carbon_SE[t_period, s_scenario].fix(0)

    if Remove_Energy_Mkt_param == 1:
        for t_period in model_instance.T:
            for s_scenario in model_instance.S:
                model_instance.P_DSO[t_period, s_scenario].fix(0)

    if Remove_Ramp_param == 1:
        model_instance.Trans_Gen_Upper_Ramp.deactivate()
        model_instance.Trans_Gen_Lower_Ramp.deactivate()
        model_instance.Deriv_Potencia_with_Ramp.deactivate()
        # model_instance.Ramp_Upper.deactivate() # These are dual variables, not constraints to deactivate directly
        # model_instance.Ramp_Lower.deactivate()
        model_instance.slack_NL11.deactivate()
        model_instance.slack_NL12.deactivate()
    else:
        model_instance.Deriv_Potencia_without_Ramp.deactivate()

    print("Finished applying execute.run logic (partially).")


if __name__ == "__main__":
    print("Loading data...")
    data_dict_loaded = load_data()
    print("Creating model instance...")
    model_inst = create_model_instance(data_dict_loaded)
    print("Applying execute.run logic...")
    apply_execute_run_logic(model_inst, data_dict_loaded)
    
    print("Model setup complete (partial). Next steps: finish data loading, execute.run logic, and solve.")

    # Example: Displaying some loaded data from the model instance
    # print("Set S:", list(model_inst.S))
    # print("Set T:", list(model_inst.T))
    # if model_inst.Trans_Lines:
    #    first_line = next(iter(model_inst.Trans_Lines))
    #    print(f"Trans_Capacity for line {first_line}: {model_inst.Trans_Capacity[first_line].value if model_inst.Trans_Capacity[first_line] is not None else 'N/A'}")
    #    print(f"Trans_Incidencia for node 1, line {first_line}: {model_inst.Trans_Incidencia[1, first_line].value if (1,first_line) in model_inst.Trans_Incidencia else 'N/A'}")

    # To solve the model (example, solver needs to be available):
    # solver = pyo.SolverFactory('glpk') # or 'ipopt', 'conopt' if available and model is NLP/MPEC
    # results = solver.solve(model_inst, tee=True)
    # print(results)
