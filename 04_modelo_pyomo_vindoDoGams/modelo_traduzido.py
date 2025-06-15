"""
## Pyomo Translation of GAMSSOURCE: https://www.gams.com/latest/psoptlib_ml/libhtml/psoptlib_MultiperiodACOPF24bus.html

Multi-period AC-OPF for IEEE 24-bus network considering wind and load shedding

Original GAMS model details:
Soroudi, Alireza. Power System Optimization Modeling in GAMS. Springer, 2017.
--------------------------------------------------------------------------------
Model type: NLP
--------------------------------------------------------------------------------
Contributed by
Dr. Alireza Soroudi
IEEE Senior Member
email: alireza.soroudi@gmail.com
"""

from __future__ import annotations

import math
import pandas as pd
import pyomo.environ as pyo

# Helper function from original script
def reformat_df(dataframe):
    return dataframe.reset_index().melt(
        id_vars="index", var_name="Category", value_name="Value"
    )

# Data records function from original script
def data_records():
    # GenD records table
    cols = ["pmax", "pmin", "b", "Qmax", "Qmin", "Vg", "RU", "RD"]
    inds = ["1", "2", "7", "13", "15", "16", "18", "21", "22", "23"]
    data = [
        [152, 30.4, 13.32, 192, -50, 1.035, 21, 21],
        [152, 30.4, 13.32, 192, -50, 1.035, 21, 21],
        [350, 75.0, 20.7, 300, 0, 1.025, 43, 43],
        [591, 206.85, 20.93, 591, 0, 1.02, 31, 31],
        [215, 66.25, 21.0, 215, -100, 1.014, 31, 31],
        [155, 54.25, 10.52, 155, -50, 1.017, 31, 31],
        [400, 100.0, 5.47, 400, -50, 1.05, 70, 70],
        [400, 100.0, 5.47, 400, -50, 1.05, 70, 70],
        [300, 0.0, 0.0, 300, -60, 1.05, 53, 53],
        [360, 248.5, 10.52, 310, -125, 1.05, 31, 31],
    ]
    GenD_recs = reformat_df(pd.DataFrame(data, columns=cols, index=inds))

    # BD records table
    cols = ["Pd", "Qd"]
    inds = [str(ii) for ii in range(1, 25)]
    data = [
        [108, 22], [97, 20], [180, 37], [74, 15], [71, 14],
        [136, 28], [125, 25], [171, 35], [175, 36], [195, 40],
        [0, 0], [0, 0], [265, 54], [194, 39], [317, 64],
        [100, 20], [0, 0], [333, 68], [181, 37], [128, 26],
        [0, 0], [0, 0], [0, 0], [0, 0],
    ]
    BD_recs = reformat_df(pd.DataFrame(data, columns=cols, index=inds))

    # LN records table
    cols_ln = ["r", "x", "b", "limit"]
    inds_ln = [
        ("1", "2"), ("1", "3"), ("1", "5"), ("2", "4"), ("2", "6"),
        ("3", "9"), ("3", "24"), ("4", "9"), ("5", "10"), ("6", "10"),
        ("7", "8"), ("8", "9"), ("8", "10"), ("9", "11"), ("9", "12"),
        ("10", "11"), ("10", "12"), ("11", "13"), ("11", "14"), ("12", "13"),
        ("12", "23"), ("13", "23"), ("14", "16"), ("15", "16"), ("15", "21"),
        ("15", "24"), ("16", "17"), ("16", "19"), ("17", "18"), ("17", "22"),
        ("18", "21"), ("19", "20"), ("20", "23"), ("21", "22"),
    ]
    data_ln = [
        [0.0026, 0.0139, 0.4611, 175], [0.0546, 0.2112, 0.0572, 175],
        [0.0218, 0.0845, 0.0229, 175], [0.0328, 0.1267, 0.0343, 175],
        [0.0497, 0.192, 0.052, 175], [0.0308, 0.119, 0.0322, 175],
        [0.0023, 0.0839, 0.0, 400], [0.0268, 0.1037, 0.0281, 175],
        [0.0228, 0.0883, 0.0239, 175], [0.0139, 0.0605, 2.459, 175],
        [0.0159, 0.0614, 0.0166, 175], [0.0427, 0.1651, 0.0447, 175],
        [0.0427, 0.1651, 0.0447, 175], [0.0023, 0.0839, 0.0, 400],
        [0.0023, 0.0839, 0.0, 400], [0.0023, 0.0839, 0.0, 400],
        [0.0023, 0.0839, 0.0, 400], [0.0061, 0.0476, 0.0999, 500],
        [0.0054, 0.0418, 0.0879, 500], [0.0061, 0.0476, 0.0999, 500],
        [0.0124, 0.0966, 0.203, 500], [0.0111, 0.0865, 0.1818, 500],
        [0.005, 0.0389, 0.0818, 500], [0.0022, 0.0173, 0.0364, 500],
        [0.00315, 0.0245, 0.206, 1000], [0.0067, 0.0519, 0.1091, 500],
        [0.0033, 0.0259, 0.0545, 500], [0.003, 0.0231, 0.0485, 500],
        [0.0018, 0.0144, 0.0303, 500], [0.0135, 0.1053, 0.2212, 500],
        [0.00165, 0.01295, 0.109, 1000], [0.00255, 0.0198, 0.1666, 1000],
        [0.0014, 0.0108, 0.091, 1000], [0.0087, 0.0678, 0.1424, 500],
    ]
    # Corrected LN_recs creation to match GAMSpy version
    _inds_ln_multi = pd.MultiIndex.from_tuples(inds_ln, names=["Index1", "Index2"])
    _df_ln = pd.DataFrame(data_ln, columns=cols_ln, index=_inds_ln_multi)
    _df_ln.reset_index(inplace=True)
    LN_recs = _df_ln.melt(
        id_vars=["Index1", "Index2"], value_vars=["r", "x", "b", "limit"],
        var_name="Category", value_name="Value" # Match column names from reformat_df
    )


    # WD records table
    cols_wd = ["w", "d"]
    inds_wd = [f"t{tt}" for tt in range(1, 25)]
    data_wd = [
        [0.0786666666666667, 0.684511335492475], [0.0866666666666667, 0.644122690036197],
        [0.117333333333333, 0.6130691560297], [0.258666666666667, 0.599733282530006],
        [0.361333333333333, 0.588874071251667], [0.566666666666667, 0.5980186702229],
        [0.650666666666667, 0.626786054486569], [0.566666666666667, 0.651743189178891],
        [0.484, 0.706039245570585], [0.548, 0.787007048961707],
        [0.757333333333333, 0.839016955610593], [0.710666666666667, 0.852733854067441],
        [0.870666666666667, 0.870642027052772], [0.932, 0.834254143646409],
        [0.966666666666667, 0.816536483139646], [1.0, 0.819394170318156],
        [0.869333333333333, 0.874071251666984], [0.665333333333333, 1.0],
        [0.656, 0.983615926843208], [0.561333333333333, 0.936368832158506],
        [0.565333333333333, 0.887597637645266], [0.556, 0.809297008954087],
        [0.724, 0.74585635359116], [0.84, 0.733473042484283],
    ]
    WD_recs = reformat_df(pd.DataFrame(data_wd, columns=cols_wd, index=inds_wd))

    return GenD_recs, BD_recs, LN_recs, WD_recs

def main():
    model = pyo.ConcreteModel(name="MultiPeriodACOPF")

    # SETS #
    model.i_nodes = pyo.Set(initialize=[str(ii) for ii in range(1, 25)], doc="network buses")
    model.slack_nodes = pyo.Set(initialize=['13'], within=model.i_nodes, doc="slack buses")
    model.t_periods = pyo.Set(initialize=[f"t{tt}" for tt in range(1, 25)], ordered=True, doc="time periods")

    # ALIAS for j (often just use model.i_nodes again) #
    model.j_nodes = pyo.Set(initialize=model.i_nodes.data(), ordered=model.i_nodes.ordered)


    # SCALARS #
    model.Sbase = pyo.Param(initialize=100, doc="System base MVA")

    # Data Loading and Preprocessing for Parameters
    _GenD_recs, _BD_recs, _LN_recs, _WD_recs = data_records()

    # PARAMETERS
    _GenD_dict = {(_rec['index'], _rec['Category']): _rec['Value'] for _, _rec in _GenD_recs.iterrows()}
    model.GenD_attrs = pyo.Set(initialize=_GenD_recs['Category'].unique().tolist())
    model.GenD_gen_nodes = pyo.Set(initialize=_GenD_recs['index'].unique().tolist(), within=model.i_nodes)
    model.GenD = pyo.Param(model.GenD_gen_nodes, model.GenD_attrs, initialize=_GenD_dict, doc="generating units characteristics")

    _BD_dict = {(_rec['index'], _rec['Category']): _rec['Value'] for _, _rec in _BD_recs.iterrows()}
    model.BD_attrs = pyo.Set(initialize=_BD_recs['Category'].unique().tolist())
    model.BD_load_nodes = pyo.Set(initialize=_BD_recs['index'].unique().tolist(), within=model.i_nodes) # All nodes have demand data
    model.BD = pyo.Param(model.BD_load_nodes, model.BD_attrs, initialize=_BD_dict, doc="demands of each bus")

    _WD_dict = {(_rec['index'], _rec['Category']): _rec['Value'] for _, _rec in _WD_recs.iterrows()}
    model.WD_attrs = pyo.Set(initialize=_WD_recs['Category'].unique().tolist())
    model.WD = pyo.Param(model.t_periods, model.WD_attrs, initialize=_WD_dict, doc="wind and demand factors per period")

    _Wcap_data_pyomo = {"8": 200, "19": 150, "21": 100}
    model.Wcap_wind_nodes = pyo.Set(initialize=list(_Wcap_data_pyomo.keys()), within=model.i_nodes)
    model.Wcap = pyo.Param(model.Wcap_wind_nodes, initialize=_Wcap_data_pyomo, default=0, doc="wind capacity at buses")

    # LN Parameter Preprocessing
    _processed_LN_data = {(_rec['Index1'], _rec['Index2'], _rec['Category']): _rec['Value']
                           for _, _rec in _LN_recs.iterrows()}
    
    _temp_LN_storage = _processed_LN_data.copy() # Work on a copy for iterative updates

    for i_n in model.i_nodes:
        for j_n in model.j_nodes:
            if i_n == j_n: continue
            for attr in ['x', 'r', 'b', 'limit']:
                val_ij = _temp_LN_storage.get((i_n, j_n, attr), 0)
                if val_ij == 0:
                    val_ji = _temp_LN_storage.get((j_n, i_n, attr), 0)
                    if val_ji != 0:
                         _temp_LN_storage[i_n, j_n, attr] = val_ji
    
    for i_n in model.i_nodes:
        for j_n in model.j_nodes:
            if i_n == j_n: continue
            limit_ij = _temp_LN_storage.get((i_n, j_n, 'limit'), 0)
            if limit_ij != 0:
                x_ij = _temp_LN_storage.get((i_n, j_n, 'x'), 0)
                r_ij = _temp_LN_storage.get((i_n, j_n, 'r'), 0)

                if x_ij != 0:
                    _temp_LN_storage[i_n, j_n, 'bij'] = 1 / x_ij
                
                _temp_LN_storage[i_n, j_n, 'z'] = math.sqrt(x_ij**2 + r_ij**2)

                if x_ij != 0 and r_ij != 0:
                    _temp_LN_storage[i_n, j_n, 'th'] = math.atan(x_ij / r_ij)
                elif x_ij != 0 and r_ij == 0: # x is non-zero, r is zero
                    _temp_LN_storage[i_n, j_n, 'th'] = math.pi / 2
                elif r_ij != 0 and x_ij == 0: # r is non-zero, x is zero
                    _temp_LN_storage[i_n, j_n, 'th'] = 0
                # else th remains unset (or 0 if x_ij and r_ij are both 0 and atan(0/0) is an issue)
                # GAMS logic implies th=0 if x=0 and r!=0. If both x,r=0, th is also 0.

    # Symmetrize z and th
    for i_n in model.i_nodes:
        for j_n in model.j_nodes:
            if i_n == j_n: continue
            # Symmetrize z
            z_ij = _temp_LN_storage.get((i_n, j_n, 'z'))
            if z_ij is not None and _temp_LN_storage.get((j_n, i_n, 'z'), 0) == 0:
                _temp_LN_storage[j_n, i_n, 'z'] = z_ij
            
            # Symmetrize th if line i-j exists (has limit_ij)
            limit_ij = _temp_LN_storage.get((i_n, j_n, 'limit'), 0)
            if limit_ij != 0:
                th_ij = _temp_LN_storage.get((i_n, j_n, 'th'))
                if th_ij is not None: # Ensure th_ij was computed
                     _temp_LN_storage[j_n, i_n, 'th'] = th_ij
    
    _ln_attrs_list = set(k[2] for k in _temp_LN_storage)
    model.LN_attrs = pyo.Set(initialize=list(_ln_attrs_list))
    model.LN = pyo.Param(model.i_nodes, model.j_nodes, model.LN_attrs, initialize=_temp_LN_storage, default=0.0, doc="network technical characteristics")

    # Line Existence Set (cx parameter in GAMS)
    model.lines = pyo.Set(within=model.i_nodes * model.j_nodes, initialize=[])
    _cx_dict = {}
    for i_n in model.i_nodes:
        for j_n in model.j_nodes:
            if i_n == j_n: continue
            # Use .get on _temp_LN_storage as model.LN might not be fully queryable yet
            limit_ij_cx = _temp_LN_storage.get((i_n, j_n, 'limit'), 0)
            limit_ji_cx = _temp_LN_storage.get((j_n, i_n, 'limit'), 0)
            if limit_ij_cx != 0 and limit_ji_cx != 0:
                _cx_dict[i_n, j_n] = 1
    
    # Symmetrize cx based on GAMS logic: cx[i,j].where[cx[j,i]] = 1
    for i_n in model.i_nodes:
        for j_n in model.j_nodes:
            if i_n == j_n: continue
            if _cx_dict.get((j_n, i_n)) == 1 and _cx_dict.get((i_n, j_n)) is None:
                 _cx_dict[i_n, j_n] = 1
    
    for k_cx, v_cx in _cx_dict.items():
        if v_cx == 1:
            model.lines.add(k_cx)

    # VARIABLES #
    model.OF = pyo.Var(doc="Objective function value")
    model.Pij = pyo.Var(model.i_nodes, model.j_nodes, model.t_periods, doc="Active power flow")
    model.Qij = pyo.Var(model.i_nodes, model.j_nodes, model.t_periods, doc="Reactive power flow")
    model.Pg = pyo.Var(model.i_nodes, model.t_periods, doc="Active power generation")
    model.Qg = pyo.Var(model.i_nodes, model.t_periods, doc="Reactive power generation")
    model.Va = pyo.Var(model.i_nodes, model.t_periods, doc="Voltage angle")
    model.V = pyo.Var(model.i_nodes, model.t_periods, doc="Voltage magnitude")
    model.Pw = pyo.Var(model.i_nodes, model.t_periods, doc="Wind power generation")

    # EQUATIONS (CONSTRAINTS) #
    def eq1_rule(model, i_n, j_n, t_p): # Pij
        if (i_n, j_n) not in model.lines: return pyo.Constraint.Skip
        z_ji = model.LN[j_n, i_n, 'z']
        if abs(z_ji) < 1e-9: return model.Pij[i_n,j_n,t_p] == 0 # Or skip, depends on desired behavior
        
        th_ji = model.LN[j_n, i_n, 'th'] # Parameters, use math.cos/sin
        return model.Pij[i_n, j_n, t_p] == \
               (model.V[i_n, t_p]**2 * math.cos(th_ji) -
                model.V[i_n, t_p] * model.V[j_n, t_p] *
                pyo.cos(model.Va[i_n, t_p] - model.Va[j_n, t_p] + th_ji)
               ) / z_ji
    model.eq1 = pyo.Constraint(model.i_nodes, model.j_nodes, model.t_periods, rule=eq1_rule)

    def eq2_rule(model, i_n, j_n, t_p): # Qij
        if (i_n, j_n) not in model.lines: return pyo.Constraint.Skip
        z_ji = model.LN[j_n, i_n, 'z']
        if abs(z_ji) < 1e-9: return model.Qij[i_n,j_n,t_p] == 0

        th_ji = model.LN[j_n, i_n, 'th']
        b_ji_half = model.LN[j_n, i_n, 'b'] / 2 # Line charging susceptance (half)
        return model.Qij[i_n, j_n, t_p] == \
               (model.V[i_n, t_p]**2 * math.sin(th_ji) -
                model.V[i_n, t_p] * model.V[j_n, t_p] *
                pyo.sin(model.Va[i_n, t_p] - model.Va[j_n, t_p] + th_ji)
               ) / z_ji - b_ji_half * model.V[i_n, t_p]**2
    model.eq2 = pyo.Constraint(model.i_nodes, model.j_nodes, model.t_periods, rule=eq2_rule)

    def eq3_rule(model, i_n, t_p): # Active Power Balance
        lhs = 0
        if i_n in model.Wcap_wind_nodes: lhs += model.Pw[i_n, t_p]
        if i_n in model.GenD_gen_nodes: lhs += model.Pg[i_n, t_p]
        
        # Demand at bus i_n
        # BD[i,"pd"] is defined for model.BD_load_nodes. Assume all i_nodes are in BD_load_nodes or BD returns 0.
        pd_val = model.BD.get((i_n, 'Pd'), 0) # Use .get for safety if BD not dense for all i_nodes
        lhs -= model.WD[t_p, 'd'] * pd_val / model.Sbase
        
        sum_Pij = sum(model.Pij[i_n, j_n, t_p] for j_n in model.j_nodes if (i_n, j_n) in model.lines)
        return lhs == sum_Pij
    model.eq3 = pyo.Constraint(model.i_nodes, model.t_periods, rule=eq3_rule)

    def eq4_rule(model, i_n, t_p): # Reactive Power Balance
        lhs = 0
        if i_n in model.GenD_gen_nodes: lhs += model.Qg[i_n, t_p]
        
        qd_val = model.BD.get((i_n, 'Qd'), 0)
        lhs -= model.WD[t_p, 'd'] * qd_val / model.Sbase
        
        sum_Qij = sum(model.Qij[i_n, j_n, t_p] for j_n in model.j_nodes if (i_n, j_n) in model.lines)
        return lhs == sum_Qij
    model.eq4 = pyo.Constraint(model.i_nodes, model.t_periods, rule=eq4_rule)

    def eq5_rule(model): # Objective Function Constraint
        cost = sum(model.Pg[i_n, t_p] * model.GenD[i_n, 'b'] * model.Sbase
                   for i_n in model.GenD_gen_nodes for t_p in model.t_periods)
        return cost <= model.OF
    model.eq5 = pyo.Constraint(rule=eq5_rule)

    def eq6_rule(model, i_n, t_p): # Ramp-Up
        if i_n not in model.GenD_gen_nodes: return pyo.Constraint.Skip
        if t_p == model.t_periods.first(): return pyo.Constraint.Skip
        t_prev = model.t_periods.prev(t_p)
        return model.Pg[i_n, t_p] - model.Pg[i_n, t_prev] <= model.GenD[i_n, 'RU'] / model.Sbase
    model.eq6 = pyo.Constraint(model.i_nodes, model.t_periods, rule=eq6_rule)

    def eq7_rule(model, i_n, t_p): # Ramp-Down
        if i_n not in model.GenD_gen_nodes: return pyo.Constraint.Skip
        if t_p == model.t_periods.last(): return pyo.Constraint.Skip
        t_next = model.t_periods.next(t_p)
        return model.Pg[i_n, t_p] - model.Pg[i_n, t_next] <= model.GenD[i_n, 'RD'] / model.Sbase
    model.eq7 = pyo.Constraint(model.i_nodes, model.t_periods, rule=eq7_rule)

    # OBJECTIVE #
    model.objective = pyo.Objective(expr=model.OF, sense=pyo.minimize)

    # VARIABLE BOUNDS AND INITIAL VALUES #
    for i_n in model.i_nodes:
        for t_p in model.t_periods:
            if i_n in model.GenD_gen_nodes:
                model.Pg[i_n, t_p].setlb(model.GenD[i_n, 'Pmin'] / model.Sbase)
                model.Pg[i_n, t_p].setub(model.GenD[i_n, 'Pmax'] / model.Sbase)
                model.Pg[i_n, t_p].value = 43000 / model.Sbase # Initial guess from GAMS
                model.Qg[i_n, t_p].setlb(model.GenD[i_n, 'Qmin'] / model.Sbase)
                model.Qg[i_n, t_p].setub(model.GenD[i_n, 'Qmax'] / model.Sbase)
            else: # Node is not a generator, fix Pg, Qg to 0
                model.Pg[i_n, t_p].fix(0)
                model.Qg[i_n, t_p].fix(0)

            model.Va[i_n, t_p].setlb(-math.pi / 2)
            model.Va[i_n, t_p].setub(math.pi / 2)
            model.Va[i_n, t_p].value = 0
            if i_n in model.slack_nodes:
                model.Va[i_n, t_p].fix(0)

            model.V[i_n, t_p].setlb(0.9)
            model.V[i_n, t_p].setub(1.1)
            model.V[i_n, t_p].value = 1.0
            
            if i_n in model.Wcap_wind_nodes:
                model.Pw[i_n, t_p].setlb(0)
                model.Pw[i_n, t_p].setub(model.WD[t_p, 'w'] * model.Wcap[i_n] / model.Sbase)
            else: # Node has no wind capacity, fix Pw to 0
                model.Pw[i_n, t_p].fix(0)

    for i_n in model.i_nodes:
        for j_n in model.j_nodes:
            for t_p in model.t_periods:
                if (i_n, j_n) in model.lines:
                    limit_val = model.LN[i_n, j_n, 'limit'] / model.Sbase
                    model.Pij[i_n, j_n, t_p].setlb(-limit_val)
                    model.Pij[i_n, j_n, t_p].setub(limit_val)
                    model.Qij[i_n, j_n, t_p].setlb(-limit_val)
                    model.Qij[i_n, j_n, t_p].setub(limit_val)
                else: # No line, fix flow to 0
                    model.Pij[i_n, j_n, t_p].fix(0)
                    model.Qij[i_n, j_n, t_p].fix(0)
    
    # SOLVE #
    # Specify CONOPT path if not in system PATH, e.g. executable=r'C:\GAMS\win64\24.7\conopt.exe'
    # Ensure you have a CONOPT license accessible by Pyomo.
    # If CONOPT is not available, try IPOPT: solver = pyo.SolverFactory('ipopt')
    solver_path_comment = r"c:\ampl\conopt.exe" # From original GAMSpy file comment
    try:
        solver = pyo.SolverFactory('conopt', executable=solver_path_comment)
    except Exception:
        print(f"Failed to load CONOPT from {solver_path_comment}, trying 'conopt' from PATH or 'ipopt'.")
        try:
            solver = pyo.SolverFactory('conopt') # Try from PATH
        except Exception:
            print("CONOPT not found in PATH, trying IPOPT.")
            solver = pyo.SolverFactory('ipopt') # Fallback

    results = solver.solve(model, tee=True)
    
    # For duals (LMPs)
    model.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)

    # REPORTING #
    print("\nSolver Status:", results.solver.status)
    print("Termination Condition:", results.solver.termination_condition)
    print("Objective Function Value: ", pyo.value(model.OF))

    report_rows = []
    for t_p in model.t_periods:
        for i_n in model.i_nodes:
            row = {'t': t_p, 'i': i_n}
            row['V'] = pyo.value(model.V[i_n, t_p])
            row['Angle'] = pyo.value(model.Va[i_n, t_p]) * 180 / math.pi # Convert angle to degrees for report
            row['Pd'] = model.BD.get((i_n, 'Pd'),0) 
            row['Qd'] = model.BD.get((i_n, 'Qd'),0)
            
            lmp_p_val = float('nan')
            lmp_q_val = float('nan')
            try: # Duals might not be available for all constraints or if solve failed
                if model.eq3[i_n, t_p].active: # Check if constraint is active
                    lmp_p_val = model.dual[model.eq3[i_n, t_p]] / model.Sbase
                if model.eq4[i_n, t_p].active:
                    lmp_q_val = model.dual[model.eq4[i_n, t_p]] / model.Sbase
            except KeyError: # Constraint might not exist (e.g. skipped)
                pass 
            except Exception as e: # Other issues like dual not found
                print(f"Warning: Could not retrieve dual for ({i_n},{t_p}): {e}")

            row['LMP_P'] = lmp_p_val
            row['LMP_Q'] = lmp_q_val
            report_rows.append(row)
    
    df_report_final = pd.DataFrame(report_rows)
    if not df_report_final.empty:
        df_report_final = df_report_final.set_index(['t', 'i']).round(4)
    print("\nReport (V, Angle(deg), Pd, Qd, LMP_P, LMP_Q):\n", df_report_final)

    # Excel Export
    output_excel_file = "results_pyomo_translation.xlsx"
    with pd.ExcelWriter(output_excel_file, engine="openpyxl") as writer_pyomo:
        if not df_report_final.empty:
            df_report_final.to_excel(writer_pyomo, sheet_name="MainReport")

        # Pg levels report (similar to GAMSpy df_Pg)
        pg_records = []
        for i_n in model.GenD_gen_nodes:
            for t_p in model.t_periods:
                pg_records.append({'i_node': i_n, 't_period': t_p, 'Pg_level': pyo.value(model.Pg[i_n, t_p]) * model.Sbase}) # Scale back to MW
        df_Pg_pyomo = pd.DataFrame(pg_records)
        
        if not df_Pg_pyomo.empty:
            df_Pg_pyomo['generator_id'] = 'G' + df_Pg_pyomo['i_node'].astype(str)
            df_Pg_pivot_pyomo = (
                df_Pg_pyomo.pivot(index='t_period', columns='generator_id', values='Pg_level')
                .sort_index(axis=1)
                .reset_index()
                .round(4)
            )
            df_Pg_pivot_pyomo.to_excel(writer_pyomo, sheet_name="Pg_MW_Levels", index=False)
        print(f"\nResults exported to {output_excel_file}")

if __name__ == "__main__":
    main()
    
a=1
