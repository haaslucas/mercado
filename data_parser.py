import pandas as pd
import re

def parse_ampl_table_param(file_path, param_name, index_cols, value_cols):
    """
    Parses a specific AMPL-style table parameter from a .dat file.
    Example: param: G_T: G_Node a_t ... := ... ;
    Assumes data ends with a semicolon.
    """
    data_dict = {}
    capture = False
    header = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(f"param: {param_name}:") or \
               (line.startswith(f"param {param_name}:=") and not index_cols): # For simple params like 'param wind1 :='
                # Try to capture header from the same line if present after ':='
                if ":=" in line:
                    header_part = line.split(":=")[0].strip()
                    # param: NAME: IDX1 IDX2 ... VAL1 VAL2 :=
                    # or param NAME := (for non-indexed or 2D array style)
                    if f"param: {param_name}:" in header_part: # Table with explicit index and value columns
                         # param: G_T: G_Node	a_t ... Pmin_t	Carbon_Cost_t	G_Owner_t:=
                        potential_header = header_part.replace(f"param: {param_name}:", "").strip().split()
                        if len(potential_header) == len(index_cols) + len(value_cols):
                            header = potential_header
                    elif f"param {param_name} :=" in line and not index_cols and len(value_cols) > 0 : # 2D array style like wind1
                        # This case is for params like 'param wind1 :=' where indices are implicit in data rows
                        # The value_cols would represent the actual data columns after implicit indices.
                        # For 'param wind1 :=' the implicit indices are S and T, value is the last column.
                        # This simplified parser might need specific handling for such cases.
                        # For now, assume index_cols handles the named index columns.
                        pass

                capture = True
                continue

            if capture:
                if line.endswith(";"):
                    line = line[:-1].strip() # Remove semicolon
                    if line: # Process line if not empty after removing semicolon
                        process_line(line, header, index_cols, value_cols, data_dict, param_name)
                    break # End of param block
                if not line or line.startswith("#"):
                    continue # Skip empty lines or comments
                
                process_line(line, header, index_cols, value_cols, data_dict, param_name)
    return data_dict

def process_line(line, header, index_cols, value_cols, data_dict, param_name):
    """ Helper to process a single data line for a table. """
    parts = re.split(r'\s+', line.strip())
    if not header: # If header wasn't on the 'param:' line, assume this is it (simplification)
        # This is a very basic heuristic, real AMPL parsing is more complex
        if len(parts) == len(index_cols) + len(value_cols):
            header = parts # This line is likely a header if not already parsed
            return 
        # If it's not a header and we don't have one, this is tricky.
        # For now, we'll assume headers are found or are implicit.

    # For 'param wind1 :=' style, parts are [idx1, idx2, val1]
    # For table style, parts correspond to header columns.
    
    current_values = {}
    # Assuming parts directly map to index_cols followed by value_cols
    # This is a simplification. AMPL .dat files can have more complex structures.
    
    num_indices = len(index_cols)
    key_parts = []
    for i in range(num_indices):
        # Convert to int if possible, else keep as string
        try:
            key_parts.append(int(parts[i]))
        except ValueError:
            key_parts.append(parts[i])
    
    # Create tuple key if multiple indices, else single value
    key = tuple(key_parts) if num_indices > 1 else key_parts[0] if num_indices == 1 else None


    if len(value_cols) == 1:
        # Convert to float if possible
        try:
            data_dict[key] = float(parts[num_indices])
        except ValueError:
            data_dict[key] = parts[num_indices] # Keep as string if not float
    else:
        # Multiple value columns
        val_dict = {}
        for i, col_name in enumerate(value_cols):
            try:
                val_dict[col_name] = float(parts[num_indices + i])
            except ValueError:
                val_dict[col_name] = parts[num_indices + i]
        data_dict[key] = val_dict


def parse_ampl_set(file_path, set_name):
    """
    Parses a specific AMPL-style set from a .dat file.
    Example: set S := 1 2 3 4 5 ;
    Returns a list of set elements.
    """
    set_elements = []
    capture = False
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(f"set {set_name} :="):
                capture = True
                line_content = line.split(":=")[1].strip()
            elif capture and line.startswith("set") and ":=" in line: # Start of a new set
                break # Stop if we encounter another set definition
            elif capture:
                line_content = line
            else:
                continue

            if capture:
                # Remove comments from the line content
                line_content = line_content.split("#")[0].strip()
                
                elements_on_line = line_content.split()
                for elem in elements_on_line:
                    if elem == ";":
                        capture = False
                        break
                    if elem: # Ensure not empty string
                        try:
                            set_elements.append(int(elem))
                        except ValueError:
                            set_elements.append(elem) # Keep as string if not int
                if not capture: # If semicolon was found
                    break
    return set_elements

def parse_ampl_simple_param_2d(file_path, param_name, index_names, value_name):
    """
    Parses AMPL parameters like: param wind1 := 1 1 0.28 ... ;
    Returns a dictionary where keys are (idx1, idx2) and value is the parameter value.
    """
    data = {}
    capture = False
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(f"param {param_name} :="):
                capture = True
                # Process rest of the line if data starts on the same line
                line_content = line.split(":=", 1)[1].strip()
                if line_content:
                    parts = re.split(r'\s+', line_content)
                    # Consume parts in groups of 3 (idx1, idx2, val) until ';'
                    while parts:
                        if parts[0] == ';':
                            capture = False
                            break
                        if len(parts) < 3: break # Not enough parts for a full entry

                        idx1, idx2, val_str = parts[0], parts[1], parts[2]
                        try:
                            key = (int(idx1), int(idx2))
                            val = float(val_str)
                            data[key] = val
                        except ValueError:
                            # Handle cases where conversion might fail or structure is different
                            # For simplicity, we'll skip malformed entries here
                            pass 
                        parts = parts[3:]
                    if not capture: break # Semicolon found
                continue

            if capture:
                if line.endswith(";"):
                    line = line[:-1].strip() # Remove semicolon
                    parts = re.split(r'\s+', line)
                    while parts:
                        if not parts[0]: # Handle potential empty strings from multiple spaces
                            parts.pop(0)
                            continue
                        if len(parts) < 3: break
                        idx1, idx2, val_str = parts[0], parts[1], parts[2]
                        try:
                            key = (int(idx1), int(idx2))
                            val = float(val_str)
                            data[key] = val
                        except ValueError:
                            pass
                        parts = parts[3:]
                    break # End of param block
                
                if not line or line.startswith("#"):
                    continue

                parts = re.split(r'\s+', line)
                while parts:
                    if not parts[0]: 
                        parts.pop(0)
                        continue
                    if len(parts) < 3: break
                    idx1, idx2, val_str = parts[0], parts[1], parts[2]
                    try:
                        key = (int(idx1), int(idx2))
                        val = float(val_str)
                        data[key] = val
                    except ValueError:
                        pass
                    parts = parts[3:]
    return data

# More parsers can be added for different AMPL data structures as needed.
# For example, for tables like 'param: T:'
# param:	T: res1 res4 ... :=
# 1	0.90	0.36 ...
# This requires knowing the column headers.

def parse_ampl_indexed_table(file_path, table_marker, index_col_name, value_col_names):
    """
    Parses tables like:
    param:	T: res1 res4 ... :=
    1	0.90	0.36 ...
    index_col_name is the name of the first column (e.g., 'T_idx')
    value_col_names is a list of names for the subsequent value columns.
    """
    data = {}
    header = []
    capture_data = False
    
    with open(file_path, 'r') as f:
        for line_num, line_text in enumerate(f):
            line = line_text.strip()

            if line.startswith(table_marker):
                # Header is usually on the same line after ':=', or on the next line
                if ":=" in line:
                    header_line_content = line.split(":=")[0].strip()
                    # param: T: res1 res4 res5 ... :=
                    # We need to extract 'res1', 'res4', etc.
                    # This part is tricky as the header is part of the marker itself.
                    # For 'param: T:', the actual value columns are listed.
                    # Let's assume value_col_names are provided correctly.
                    header = [index_col_name] + value_col_names # Construct expected header
                capture_data = True
                continue

            if capture_data:
                if not line or line.startswith("#"): # Skip comments or empty lines
                    continue
                if line.endswith(";"):
                    line = line[:-1].strip() # Remove semicolon
                    if line: # Process line if not empty
                        parts = re.split(r'\s+', line)
                        if len(parts) == len(header):
                            idx_val = int(parts[0])
                            data[idx_val] = {}
                            for i, col_name in enumerate(value_col_names):
                                try:
                                    data[idx_val][col_name] = float(parts[i+1])
                                except ValueError:
                                    data[idx_val][col_name] = parts[i+1]
                    break # End of table

                parts = re.split(r'\s+', line)
                if len(parts) == len(header): # Ensure line matches header structure
                    idx_val = int(parts[0])
                    data[idx_val] = {}
                    for i, col_name in enumerate(value_col_names):
                        try:
                            data[idx_val][col_name] = float(parts[i+1])
                        except ValueError:
                            data[idx_val][col_name] = parts[i+1]
    return data

def parse_general_table(file_path, param_marker, column_names):
    """
    A more general table parser.
    param_marker: e.g., "param: G_T:"
    column_names: list of column names in order. The first N are indices, rest are values.
                  Caller needs to specify how many are index columns.
    Returns a list of dictionaries, where each dictionary is a row.
    """
    rows = []
    capture = False
    header_parsed = False

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(param_marker):
                capture = True
                # Check if column names are on this line after ':='
                if ":=" in line:
                    # Example: param: G_T: G_Node a_t ... :=
                    # The provided column_names should match these.
                    # This parser assumes column_names is the ground truth.
                    header_parsed = True # Assume header is implicitly known by column_names
                continue

            if capture:
                if not line or line.startswith("#"):
                    continue
                if line.endswith(";"):
                    line = line[:-1].strip()
                    if line:
                        parts = re.split(r'\s+', line)
                        if len(parts) == len(column_names):
                            rows.append(dict(zip(column_names, parts)))
                    break
                
                parts = re.split(r'\s+', line)
                if not header_parsed: # If column names were not on the param_marker line
                                      # and this line looks like a header
                    if len(parts) == len(column_names) and not any(c.isdigit() for c in parts[0]): # Heuristic for header
                        # This is a very weak check for a header row
                        header_parsed = True 
                        continue # Skip this line as it's treated as header
                
                if len(parts) == len(column_names):
                    rows.append(dict(zip(column_names, parts)))
    return rows

def convert_df_to_param_dict(df, index_cols, value_col):
    """Converts a DataFrame to a dictionary suitable for Pyomo Param initialization."""
    param_dict = {}
    if isinstance(index_cols, str): # Single index column
        for _, row in df.iterrows():
            try:
                idx = int(row[index_cols])
            except ValueError:
                idx = row[index_cols]
            try:
                param_dict[idx] = float(row[value_col])
            except ValueError:
                 param_dict[idx] = row[value_col]

    else: # Multiple index columns
        for _, row in df.iterrows():
            idx_parts = []
            for col in index_cols:
                try:
                    idx_parts.append(int(row[col]))
                except ValueError:
                    idx_parts.append(row[col])
            idx = tuple(idx_parts)
            try:
                param_dict[idx] = float(row[value_col])
            except ValueError:
                param_dict[idx] = row[value_col]
    return param_dict

def convert_df_to_multi_value_param_dict(df, index_cols, value_cols_map):
    """
    Converts a DataFrame to a dictionary for Params with multiple values per index.
    value_cols_map: a dictionary like {'pyomo_param_name': 'df_column_name'}
    Result: { (idx_tuple): {'pyomo_param_name1': val1, 'pyomo_param_name2': val2} }
    """
    param_dict = {}
    if isinstance(index_cols, str): # Single index column
        for _, row in df.iterrows():
            try:
                idx = int(row[index_cols])
            except ValueError:
                idx = row[index_cols]
            
            values = {}
            for pyomo_name, df_col in value_cols_map.items():
                try:
                    values[pyomo_name] = float(row[df_col])
                except ValueError:
                    values[pyomo_name] = row[df_col]
            param_dict[idx] = values
    else: # Multiple index columns
        for _, row in df.iterrows():
            idx_parts = []
            for col in index_cols:
                try:
                    idx_parts.append(int(row[col]))
                except ValueError:
                    idx_parts.append(row[col])
            idx = tuple(idx_parts)

            values = {}
            for pyomo_name, df_col in value_cols_map.items():
                try:
                    values[pyomo_name] = float(row[df_col])
                except ValueError:
                    values[pyomo_name] = row[df_col]
            param_dict[idx] = values
    return param_dict
