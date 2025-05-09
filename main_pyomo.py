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
INPUT_DAT_PATH = os.path.join(BASE_DIR, "ampl", "input.dat")
SCENARIOS_DAT_PATH = os.path.join(BASE_DIR, "ampl", "Reducao Cenarios", "Reducao Cenarios", "Scenarios.dat")


def load_data():
    """
    Loads all data from .dat files into Python dictionaries.
    This is a simplified version and will need to be expanded.
    """
    data = {}
    
    # Load from Scenarios.dat
    data['S'] = data_parser.parse_ampl_set(SCENARIOS_DAT_PATH, "S")
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
    data['T'] = list(t_data_rows.keys()) # Set T
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

    # Populate Sets
    instance.S = pyo.Set(initialize=data_dict.get('S', []))
    instance.T = pyo.Set(initialize=data_dict.get('T', []), ordered=True) # Ensure T is ordered
    instance.Trans_Lines = pyo.Set(initialize=data_dict.get('Trans_Lines', []))
    instance.Trans_Nodes = pyo.Set(initialize=data_dict.get('Trans_Nodes', []))
    instance.N = pyo.Set(initialize=data_dict.get('N', []))
    instance.L = pyo.Set(initialize=data_dict.get('L', []), dimen=2) # dimen=2 for pairs

    # ... (Populate other sets: G_T, G_D, LS, ESS, RES_T, RES_D, LDA)

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
    
    # ... (Populate other params)

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
        if model_instance.Imax[line_key] is not None and model_instance.IBASE.value != 0:
             model_instance.Imax[line_key] = model_instance.Imax[line_key].value / model_instance.IBASE.value
             
    # ... (Implement scaling for ESS params, Generator params, Carbon Costs)
    # ... (Implement Emission_Weighted_Average, Carbon_Limit_Trans, Carbon_Limit_Dist calculations)
    
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
