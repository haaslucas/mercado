from pyomo.environ import *
import pandas as pd
import numpy as np
import sys
sys.path.append('.')
from auxiliares2 import *

def ieee34bus(  prelog=False, poslog=False, analyze_model=True, Y = 40):  # 

    m = ConcreteModel(name='34 bus',
            doc='34 bus')

    m.dual = Suffix(direction=Suffix.IMPORT)

    SBASE = 100  # Base power in MVA
    VBASE = 24.9  # Base voltage in kV
    IBASE = SBASE * 1e3 / VBASE  # Base current in A
    

    ############ PROCESSAMENTO DE DADOS ############
    PASTA = './modelos/data/'
    DLINES  = pd.read_csv(PASTA + 'distribution_line_data.csv')
    DLINES = DLINES.set_index(['From', 'To'])

    DLINES['Z (pu)'] = np.sqrt(np.sqrt((DLINES['R (pu)']**2 + DLINES['X (pu)']**2).values))
    DLINES['Line Ampacity (pu)'] = DLINES['Line Ampacity (A)'] / IBASE
    #dnodes = pd.read_csv(PASTA + 'distribution_node_data.csv', index_col=0)
    #dnodes = pd.read_csv(PASTA + 'distribution_node_data.csv')
    ess = pd.read_csv(PASTA + 'ess.csv',index_col=0)

    ess[['Max ESS Charge (pu)', 'Max ESS Discharge (pu)', 'Min ESS SOC (pu)',
       'Max ESS SOC (pu)', 'Initial SOC (pu)',]] =     ess[['Max ESS Charge (MW)', 'Max ESS Discharge (MW)', 'Min ESS SOC (MW)',
       'Max ESS SOC (MW)', 'Initial SOC (MW)',]]/ SBASE

    load = pd.read_csv(PASTA + 'dist_load.csv')
    load.Load  = load.Load*3.5
    load = load.set_index(['Node','Hour'])
    load = load[load.Scenario == 1]  
    load.drop(columns=['Scenario'], inplace=True)
    
    gen_nd = pd.read_csv(PASTA + 'dso_connected_non_dispatchable_generators.csv')

    gen_d  = pd.read_csv(PASTA + 'dso_connected_dispatchable_generators.csv')

    
    gen_d['P_max (pu)'] = gen_d['P_max (MW)'] / SBASE
    gen_d['P_min (pu)'] = gen_d['P_min (MW)'] / SBASE
    gen_nd['P_max (pu)'] = gen_nd['P_max (MW)'] / SBASE
    
    #concat  gdd and gdnd. They have column Node in common
    gen_d.index = gen_d['Node']
    gen_d.drop(columns=['Node','Unnamed: 0'], inplace=True)
    gen_d.sort_index(inplace=True)
    ls = pd.read_csv(PASTA + 'ls.csv')

    gen_d.a = gen_d.a * SBASE**2
    gen_d.b = gen_d.b * SBASE
    
    gen_nd.index = gen_nd['Node']
    gen_nd.drop(columns=['Node','Unnamed: 0', 'Wind min Velocity (m/s)', 'Wind max Velocity (m/s)',
                        'Wind Rated Velocity (m/s)'], inplace=True)
    ####################################################
    
    #################### Conjuntos #####################
        
    m.i = RangeSet(1, 34, name='i', doc='Nós da rede de distribuição')
    m.t = RangeSet(1, 24, name='t', doc='Horas do dia')
    m.GB = Set(initialize=[5, 15, 19, 24, 34], name='GB', doc='Geradores despacháveis')
    m.NGB = Set(initialize=[13, 23, 33], name='NGB', doc='Geradores não despacháveis')
    m.l = Set(initialize=[(int(row.name[0]), int(row.name[1])) for _, row in DLINES.iterrows()], name='l', doc='Linhas de distribuição')
    #LMPs não nulos: 5,9,15,19,22,34
    
    m.SBASE = Param(initialize=100, name='SBASE', doc='Base de potência do sistema') # ,units='MVA'
    m.VsqBASE = Param(initialize=24.90, name='VBASE', doc='Tensão base do sistema') # ,units='kV'
    m.IsqBASE = Param(initialize=IBASE, name='IBASE', doc='Base de corrente do sistema') # ,units='A'

    m.slack = Param(initialize=1, name='slack', doc='Nó de referência da tensão (slack bus)')
    #m.ESS_CAPACITY = Param(initialize=50/SBASE, name='SE_Capacity', doc='Capacidade do ESS') # ,units='pu'


    #################### Parametros ####################

    m.SOC0 = Param(initialize=0.01, name='SOC0', doc='Estado inicial de carga do ESS') # ,units='pu'
    
    m.GenD = Param(m.GB, ['a', 'b', 'P_max (pu)', 'P_min (pu)'], 
                        initialize=lambda model, i, j: gen_d.loc[i][j],
                        name='GenD', doc='Parâmetros dos geradores despacháveis')
    
    m.GenND = Param(m.NGB, ['P_max (pu)'], 
                        initialize=lambda model, i,j: gen_nd.loc[i][j],
                        name='GenND', doc='Parâmetros dos geradores não despacháveis')
    
    m.Load = Param(m.i, m.t, ['Load'] , 
                        initialize=lambda model, i, t, l: load.loc[(i, t), l],
                        name='Load', doc='Carga nos nós da rede de distribuição',
                        )

    m.Pij = Var(m.l, m.t, within=Reals, 
                    doc='Potência ativa nas linhas de distribuição',
                    name='Pij', 
                    )
    
    m.Qij = Var(m.l, m.t, within=Reals,
                    doc='Potência reativa nas linhas de distribuição',
                    name='Qij', 
                    )
    
    m.Isq = Var(m.l, m.t, within=NonNegativeReals,
                doc='Corrente nas linhas de distribuição',
                name='I') # , units='pu'
    
    m.Vsq = Var(m.i, m.t, 
                doc='Tensão nos nós da rede de distribuição',
                name='V') # , units='pu'

    m.Pg = Var(m.GB, m.t, 
                doc='Potência ativa dos geradores despacháveis',
                name='Pg') #, units='pu'

    m.Qg = Var(m.GB, m.t,
                doc='Potência reativa dos geradores despacháveis',
                name='Qg') # , units='pu'
    
    m.Ppv = Var(m.NGB, m.t,
                doc='Potência ativa dos geradores fotovoltaicos',
                name='Ppv') # , units='pu'
    
    m.Pbess = Var(ess.index, m.t,
                doc='Potência ativa do ESS',
                name='Pbess') # , units='pu'
    
    m.SOC = Var(ess.index, m.t,
                doc='Estado de carga do ESS',
                name='SOC') # , units='pu'
    
    #m.GenCost = Var(m.GB, m.t,
    #            doc='Custo de geração dos geradores despacháveis',
    #            name='GenCost', units='R$')
    
    #################### RESTRIÇÕES ####################
    
    m.Isq_BOUND_UPPER = Constraint(m.l, m.t,
                rule=lambda model, i, j, t: m.Isq[i, j, t] <= DLINES.loc[(i, j), 'Line Ampacity (pu)']**2,
                doc='Upper bound for current flow on distribution lines',
                name='I_BOUND_UPPER')
    
    m.Isq_BOUND_LOWER = Constraint(m.l, m.t,
                rule=lambda model, i, j, t: m.Isq[i, j, t] >= 0,
                doc='Lower bound for current flow on distribution lines',
                name='I_BOUND_LOWER')
    
    m.Vsq_BOUND_UPPER = Constraint(m.i, m.t,
                rule=lambda model, i, t: m.Vsq[i, t] <= 1.05**2,
                doc='Upper bound for voltage at each bus',
                name='V_BOUND_UPPER')
    
    m.Vsq_BOUND_LOWER = Constraint(m.i, m.t,
                rule=lambda model, i, t: m.Vsq[i, t] >= 0.95**2,
                doc='Lower bound for voltage at each bus',
                name='V_BOUND_LOWER')

    m.Vsq_SLACK = Constraint(m.t,
                rule=lambda model, t: m.Vsq[m.slack, t] == 1,
                doc='Slack bus voltage constraint',
                name='V_SLACK')
    
    m.PG_BOUNDS_UPPER = Constraint(m.GB, m.t,
                rule=lambda model, i, t: m.Pg[i, t] <= m.GenD[i, 'P_max (pu)'],
                doc='Upper bound for generation at each bus',
                name='PG_BOUNDS_UPPER')
    
    m.PG_BOUNDS_LOWER = Constraint(m.GB, m.t,
                rule=lambda model, i, t: m.Pg[i, t] >= m.GenD[i, 'P_min (pu)'],
                doc='Lower bound for generation at each bus',
                name='PG_BOUNDS_LOWER')
    
    m.QG_BOUNDS_UPPER = Constraint(m.GB, m.t,
                rule=lambda model, i, t: m.Qg[i, t] <= m.GenD[i, 'P_max (pu)'],
                doc='Upper bound for reactive power generation at each bus',
                name='QG_BOUNDS_UPPER')
    
    m.QG_BOUNDS_LOWER = Constraint(m.GB, m.t,
                rule=lambda model, i, t: m.Qg[i, t] >= -m.GenD[i, 'P_max (pu)'],
                doc='Lower bound for reactive power generation at each bus',
                name='QG_BOUNDS_LOWER')

    m.PPV_BOUNDS_UPPER = Constraint(m.NGB, m.t,
                rule=lambda model, i, t: m.Ppv[i, t] <= m.GenND[i, 'P_max (pu)'],
                doc='Upper bound for non-dispatchable generation at each bus',
                name='PPV_BOUNDS_UPPER')
    
    m.PPV_BOUNDS_LOWER = Constraint(m.NGB, m.t,
                rule=lambda model, i, t: m.Ppv[i, t] >= 0,
                doc='Lower bound for non-dispatchable generation at each bus',
                name='PPV_BOUNDS_LOWER')
    m.PBESS_BOUNDS_UPPER = Constraint(ess.index, m.t,
                rule=lambda model, i, t: m.Pbess[i, t] <= ess.loc[i, 'Max ESS Charge (pu)'],
                doc='Upper bound for ESS charging at each bus',
                name='PBESS_BOUNDS_UPPER')
    
    m.PBESS_BOUNDS_LOWER = Constraint(ess.index, m.t,
                rule=lambda model, i, t: m.Pbess[i, t] >= -ess.loc[i, 'Max ESS Charge (pu)'],
                doc='Lower bound for ESS discharging at each bus',
                name='PBESS_BOUNDS_LOWER')

    m.SOC_BOUNDS_UPPER = Constraint(ess.index, m.t,
                rule=lambda model, i, t: m.SOC[i, t] <= ess.loc[i, 'Max ESS SOC (pu)'],
                doc='Upper bound for ESS state of charge at each bus',
                name='SOC_BOUNDS_UPPER')
    
    m.SOC_BOUNDS_LOWER = Constraint(ess.index, m.t,
                rule=lambda model, i, t: m.SOC[i, t] >= ess.loc[i, 'Min ESS SOC (pu)'],
                doc='Lower bound for ESS state of charge at each bus',
                name='SOC_BOUNDS_LOWER')


    
#m.ObjInferior é a soma dos custos dos geradores mais a soma da Var 
# lambda_BALANCO_POT_ATIVA (criada mais pra frente com a função create_dual_variables),
# vezes a potência de cada gerador fotovoltaico (m.NGB).

    def balanco_pot_ativa(model, i, t): 
        injecoes = (
            (m.Pg[i, t]   if i in m.GB        else 0)
        + (m.Ppv[i, t]  if i in m.NGB       else 0)
        + (m.Pbess[i, t] if i in ess.index  else 0)
        -  m.Load[i, t, 'Load']
        )

        saindo =   sum(m.Pij[i, j, t] + DLINES.loc[(i, j), 'R (pu)']
                    * m.Isq[i, j, t] for j in m.i if (i, j) in m.l)
        
        entrando = sum(m.Pij[j, i, t] for j in m.i if (j, i) in m.l)
        
        return injecoes == saindo - entrando

    def balanco_pot_reativa(model, i, t):
        injecoes =   ( m.Qg[i, t] if i in m.GB else 0 ) \
                     - m.Load[i, t,'Load']*np.tan(np.arccos(.9))
                   
        saindo   = sum(m.Qij[i, j, t] + DLINES.loc[(i, j), 'X (pu)'] 
                     * m.Isq[i, j, t] for j in m.i if (i, j) in m.l)
                        
        entrando = sum(m.Qij[j, i, t] for j in m.i if (j, i) in m.l)
        
        return injecoes == saindo - entrando

    m.BALANCO_POT_ATIVA = Constraint(m.i, m.t, rule=balanco_pot_ativa,
                doc='Balanço de potência ativa nos nós da rede de distribuição',
                name='balanco_pot_ativa')
    
    m.BALANCO_POT_REATIVA = Constraint(m.i, m.t, rule=balanco_pot_reativa,
                doc='Balanço de potência reativa nos nós da rede de distribuição',
                name='balanco_pot_reativa')
    
    def fluxo_Linhas(model, i, j, t):
        if (i, j) in DLINES.index:
            return m.Vsq[i, t] - m.Vsq[j, t] == 2 * (DLINES.loc[(i, j), 'R (pu)'] * m.Pij[i, j, t]
                + DLINES.loc[(i, j), 'X (pu)'] * m.Qij[i, j, t]) + DLINES.loc[(i,j),'Z (pu)']**2 * m.Isq[i, j, t]
        return Constraint.Skip

    m.FLUXO_LINHAS = Constraint(m.l, m.t, rule=fluxo_Linhas,
                doc='Equação de fluxo de potência ativa e reativa nas linhas de distribuição',
                name='eq6')

    # Applying the new piecewise linear approximation for signed variables
    add_piecewise_tangent_approximation_signed(
        model=m,
        P_to_linearize=m.Pij,
        P_max_abs_func=lambda i, j, t: DLINES.loc[(i, j), 'Line Ampacity (pu)'],
        Y=Y,
        base_name='Pij',
        output_var_name='Pij_sq_piecewised_signed'
    )

    add_piecewise_tangent_approximation_signed(
        model=m,
        P_to_linearize=m.Qij,
        P_max_abs_func=lambda i, j, t: DLINES.loc[(i, j), 'Line Ampacity (pu)'],
        Y=Y,
        base_name='Qij',
        output_var_name='Qij_sq_piecewised_signed'
    )
    
    def fluxo_Linhas2_modified(model, i, j, t):
        if (i, j) in DLINES.index:
            # Use the new output variables for the approximated squared terms
            return m.Pij_sq_piecewised_signed[i, j, t] + m.Qij_sq_piecewised_signed[i, j, t] <=m.Isq[i, j, t]
        return Constraint.Skip
        

    m.FLUXO_LINHAS2 = Constraint(m.l, m.t, rule=fluxo_Linhas2_modified,
                doc='Limite de potência aparente nas linhas de distribuição (linearized)',
                name='eq7_linearized')


    def state_of_charge(model, i, t):
        if t == 1:
            
            return m.SOC[i, t] == m.SOC0
        else:
            return m.SOC[i, t] == m.SOC[i, t-1] - (m.Pbess[i, t-1])

    m.STATE_OF_CHARGE = Constraint(ess.index, m.t, rule=state_of_charge,
                doc='Estado de carga do ESS ao longo do tempo',
                name='State of Charge')

    '''def fluxo_Linhas2_modified(model, i, j, t):
        if (i, j) in DLINES.index:
            # Use the new output variables for the approximated squared terms
            return m.Pij[i, j, t] + m.Qij[i, j, t] <= m.Isq[i, j, t] * 1
        return Constraint.Skip
        
    m.FLUXO_LINHAS2 = Constraint(m.l, m.t, rule=fluxo_Linhas2_modified,
                doc='Limite de potência aparente nas linhas de distribuição (linearized)',
                name='eq7_linearized')'''
    

    # The original code had these calls for Pg and Qg, which are commented out.
    # If Pg and Qg can be negative, these should also use the signed version.
    # For now, assuming Pg and Qg are positive or zero, the original function is fine.
    # If they can be negative, uncomment and use the signed version.

    
    add_piecewise_tangent_approximation(
        model=m,
        P_to_linearize=m.Pg,
        P_max_func=lambda i, t: m.GenD[i, 'P_max (pu)'],
        Y=Y,
        base_name='Pg',
        output_var_name='Pg_sq_piecewised'
    )
    add_piecewise_tangent_approximation(
        model=m,
        P_to_linearize=m.Qg,
        P_max_func=lambda i, t: m.GenD[i, 'P_max (pu)'],
        Y=Y,
        base_name='Qg',
        output_var_name='Qg_sq_piecewised'
    )

    def gd_upper_bound(model, g, t):
        # This constraint uses the output variables from the piecewise approximation.
        # Ensure these output variables are correctly named.
        return m.Pg_sq_piecewised[g, t] + m.Qg_sq_piecewised[g, t] <= m.GenD[g, 'P_max (pu)']**2

    m.GD_UPPER_BOUND = Constraint(gen_d.index, m.t, rule=gd_upper_bound,
                doc='Limite de potência aparente dos geradores despacháveis',
                name='eq10')
    
    '''def gnd_upper_bound(model, r, t):
        # NON-LINEAR: This is a quadratic constraint on the active power of non-dispatchable generators.
        return m.Ppv[r, t]**2 <= np.tan(np.arccos(.9))*m.GenND[r, 'P_max (pu)']**2
    
    m.GND_UPPER_BOUND = Constraint(gen_nd.index, m.t, rule=gnd_upper_bound,
                doc='Limite de potência aparente dos geradores não despacháveis',
                name='eq11')'''


    # --- KKT Conditions for the Lower-Level Problem (Manual Formulation) ---

    # The objective function for the lower-level problem remains the same.
    m.ObjInferior = Objective(
        expr=(
            sum(m.Pg_sq_piecewised[i, t] * m.GenD[i, 'a'] +
                m.Pg[i, t] * m.GenD[i, 'b']
                for i in m.GB for t in m.t)
        ),
        sense=minimize,
        name='OF'
    )

    # Big-M parameter for linearizing complementarity conditions
    m.M = Param(initialize=10000, doc="Big-M for KKT linearization")

    #################################################################
    # 1. DUAL VARIABLES DECLARATION
    #################################################################
    # For each constraint in the lower-level problem, we declare a corresponding dual variable.
    # lambda for equality constraints (domain: Reals)
    # mu for inequality constraints (domain: NonNegativeReals)

    # Power Flow and System Limits
    m.mu_Isq_UP = Var(m.l, m.t, within=NonNegativeReals)
    m.mu_Isq_LOW = Var(m.l, m.t, within=NonNegativeReals)
    m.mu_Vsq_UP = Var(m.i, m.t, within=NonNegativeReals)
    m.mu_Vsq_LOW = Var(m.i, m.t, within=NonNegativeReals)
    m.lambda_Vsq_SLACK = Var(m.t, within=Reals)
    m.lambda_BPA = Var(m.i, m.t, within=Reals)
    m.lambda_BPR = Var(m.i, m.t, within=Reals)
    m.lambda_FL = Var(m.l, m.t, within=Reals)
    m.mu_FL2 = Var(m.l, m.t, within=NonNegativeReals)

    # Generation and Storage Limits
    m.mu_Pg_UP = Var(m.GB, m.t, within=NonNegativeReals)
    m.mu_Pg_LOW = Var(m.GB, m.t, within=NonNegativeReals)
    m.mu_Qg_UP = Var(m.GB, m.t, within=NonNegativeReals)
    m.mu_Qg_LOW = Var(m.GB, m.t, within=NonNegativeReals)
    m.mu_Ppv_UP = Var(m.NGB, m.t, within=NonNegativeReals)
    m.mu_Ppv_LOW = Var(m.NGB, m.t, within=NonNegativeReals)
    m.mu_Pbess_UP = Var(ess.index, m.t, within=NonNegativeReals)
    m.mu_Pbess_LOW = Var(ess.index, m.t, within=NonNegativeReals)
    m.mu_SOC_UP = Var(ess.index, m.t, within=NonNegativeReals)
    m.mu_SOC_LOW = Var(ess.index, m.t, within=NonNegativeReals)
    m.lambda_SOC = Var(ess.index, m.t, within=Reals)
    m.mu_GD_UP = Var(gen_d.index, m.t, within=NonNegativeReals)

    # Duals for Piecewise Linearization Constraints
    # For Pij
    m.lambda_p_decomp_Pij = Var(m.l, m.t, within=Reals)
    m.lambda_abs_p_Pij = Var(m.l, m.t, within=Reals)
    m.mu_pw_tan_Pij = Var(m.l, m.t, m.Y_set_Pij_abs, within=NonNegativeReals)
    # For Qij
    m.lambda_p_decomp_Qij = Var(m.l, m.t, within=Reals)
    m.lambda_abs_p_Qij = Var(m.l, m.t, within=Reals)
    m.mu_pw_tan_Qij = Var(m.l, m.t, m.Y_set_Qij_abs, within=NonNegativeReals)
    # For Pg
    m.mu_pw_tan_Pg = Var(m.GB, m.t, m.Y_set_Pg, within=NonNegativeReals)
    # For Qg
    m.mu_pw_tan_Qg = Var(m.GB, m.t, m.Y_set_Qg, within=NonNegativeReals)

    # --- Build Lagrangian Expression ---
    # This expression is the foundation for the stationarity conditions.
    def lagrangian_rule(m):
        # Start with the objective function of the lower-level problem
        expr = m.ObjInferior.expr

        # Add terms for each constraint multiplied by its dual variable
        # The canonical_form helper converts constraints to 'g(x) <= 0' or 'h(x) == 0'
        
        # Power Flow and System Limits
        expr += sum(m.mu_Isq_UP[idx] * canonical_form(m.Isq_BOUND_UPPER[idx])[0] for idx in m.mu_Isq_UP)
        expr += sum(m.mu_Isq_LOW[idx] * canonical_form(m.Isq_BOUND_LOWER[idx])[0] for idx in m.mu_Isq_LOW)
        expr += sum(m.mu_Vsq_UP[idx] * canonical_form(m.Vsq_BOUND_UPPER[idx])[0] for idx in m.mu_Vsq_UP)
        expr += sum(m.mu_Vsq_LOW[idx] * canonical_form(m.Vsq_BOUND_LOWER[idx])[0] for idx in m.mu_Vsq_LOW)
        expr += sum(m.lambda_Vsq_SLACK[idx] * canonical_form(m.Vsq_SLACK[idx])[0] for idx in m.lambda_Vsq_SLACK)
        expr += sum(m.lambda_BPA[idx] * canonical_form(m.BALANCO_POT_ATIVA[idx])[0] for idx in m.lambda_BPA)
        expr += sum(m.lambda_BPR[idx] * canonical_form(m.BALANCO_POT_REATIVA[idx])[0] for idx in m.lambda_BPR)
        expr += sum(m.lambda_FL[idx] * canonical_form(m.FLUXO_LINHAS[idx])[0] for idx in m.lambda_FL)
        expr += sum(m.mu_FL2[idx] * canonical_form(m.FLUXO_LINHAS2[idx])[0] for idx in m.mu_FL2)

        # Generation and Storage Limits
        expr += sum(m.mu_Pg_UP[idx] * canonical_form(m.PG_BOUNDS_UPPER[idx])[0] for idx in m.mu_Pg_UP)
        expr += sum(m.mu_Pg_LOW[idx] * canonical_form(m.PG_BOUNDS_LOWER[idx])[0] for idx in m.mu_Pg_LOW)
        expr += sum(m.mu_Qg_UP[idx] * canonical_form(m.QG_BOUNDS_UPPER[idx])[0] for idx in m.mu_Qg_UP)
        expr += sum(m.mu_Qg_LOW[idx] * canonical_form(m.QG_BOUNDS_LOWER[idx])[0] for idx in m.mu_Qg_LOW)
        expr += sum(m.mu_Ppv_UP[idx] * canonical_form(m.PPV_BOUNDS_UPPER[idx])[0] for idx in m.mu_Ppv_UP)
        expr += sum(m.mu_Ppv_LOW[idx] * canonical_form(m.PPV_BOUNDS_LOWER[idx])[0] for idx in m.mu_Ppv_LOW)
        expr += sum(m.mu_Pbess_UP[idx] * canonical_form(m.PBESS_BOUNDS_UPPER[idx])[0] for idx in m.mu_Pbess_UP)
        expr += sum(m.mu_Pbess_LOW[idx] * canonical_form(m.PBESS_BOUNDS_LOWER[idx])[0] for idx in m.mu_Pbess_LOW)
        expr += sum(m.mu_SOC_UP[idx] * canonical_form(m.SOC_BOUNDS_UPPER[idx])[0] for idx in m.mu_SOC_UP)
        expr += sum(m.mu_SOC_LOW[idx] * canonical_form(m.SOC_BOUNDS_LOWER[idx])[0] for idx in m.mu_SOC_LOW)
        expr += sum(m.lambda_SOC[idx] * canonical_form(m.STATE_OF_CHARGE[idx])[0] for idx in m.lambda_SOC)
        expr += sum(m.mu_GD_UP[idx] * canonical_form(m.GD_UPPER_BOUND[idx])[0] for idx in m.mu_GD_UP)

        # Duals for Piecewise Linearization Constraints
        # For Pij
        expr += sum(m.lambda_p_decomp_Pij[idx] * canonical_form(m.p_decomposition_constraint_Pij[idx])[0] for idx in m.lambda_p_decomp_Pij)
        expr += sum(m.lambda_abs_p_Pij[idx] * canonical_form(m.abs_p_constraint_Pij[idx])[0] for idx in m.lambda_abs_p_Pij)
        expr += sum(m.mu_pw_tan_Pij[idx] * canonical_form(m.pw_tangent_constr_Pij_abs[idx])[0] for idx in m.mu_pw_tan_Pij)
        # For Qij
        expr += sum(m.lambda_p_decomp_Qij[idx] * canonical_form(m.p_decomposition_constraint_Qij[idx])[0] for idx in m.lambda_p_decomp_Qij)
        expr += sum(m.lambda_abs_p_Qij[idx] * canonical_form(m.abs_p_constraint_Qij[idx])[0] for idx in m.lambda_abs_p_Qij)
        expr += sum(m.mu_pw_tan_Qij[idx] * canonical_form(m.pw_tangent_constr_Qij_abs[idx])[0] for idx in m.mu_pw_tan_Qij)
        # For Pg
        expr += sum(m.mu_pw_tan_Pg[idx] * canonical_form(m.pw_tangent_constr_Pg[idx])[0] for idx in m.mu_pw_tan_Pg)
        # For Qg
        expr += sum(m.mu_pw_tan_Qg[idx] * canonical_form(m.pw_tangent_constr_Qg[idx])[0] for idx in m.mu_pw_tan_Qg)

        return expr

    m.Lagr = Expression(rule=lagrangian_rule)


    #################################################################
    # 2. STATIONARITY CONDITIONS
    #################################################################
    # For each primal variable, the derivative of the Lagrangian must be zero.
    # dL/dx = d(Obj)/dx + sum(dual_i * d(constraint_i)/dx) = 0

    # Note: Canonical form for constraints:
    # body <= upper  -->  body - upper <= 0
    # body >= lower  -->  lower - body <= 0
    # body == equal  -->  body - equal == 0

    # Derivative w.r.t. Pg
    def dLagr_dPg_rule(m, g, t):
        return (m.GenD[g, 'b']                                               # from ObjInferior
                + m.lambda_BPA[g, t]                                         # from BALANCO_POT_ATIVA
                + m.mu_Pg_UP[g, t]                                           # from PG_BOUNDS_UPPER
                - m.mu_Pg_LOW[g, t]                                          # from PG_BOUNDS_LOWER
                - sum(m.mu_pw_tan_Pg[g, t, k] * m.slope_Pg[g, t, k] for k in m.Y_set_Pg) # from pw_tangent_constr_Pg
               ) == 0
    m.dLagr_dPg = Constraint(m.GB, m.t, rule=dLagr_dPg_rule)

    # Derivative w.r.t. Pg_sq_piecewised
    def dLagr_dPg_sq_rule(m, g, t):
        return (m.GenD[g, 'a']                                               # from ObjInferior
                + m.mu_GD_UP[g, t]                                           # from GD_UPPER_BOUND
                - sum(m.mu_pw_tan_Pg[g, t, k] for k in m.Y_set_Pg)            # from pw_tangent_constr_Pg
               ) == 0
    m.dLagr_dPg_sq = Constraint(m.GB, m.t, rule=dLagr_dPg_sq_rule)

    # Derivative w.r.t. Pij
    def dLagr_dPij_rule(m, i, j, t):
        return (-m.lambda_BPA[i, t]                                          # from BALANCO_POT_ATIVA (saindo)
                + m.lambda_BPA[j, t]                                         # from BALANCO_POT_ATIVA (entrando)
                - m.lambda_FL[i, j, t] * 2 * DLINES.loc[(i, j), 'R (pu)']     # from FLUXO_LINHAS
                + m.lambda_p_decomp_Pij[i, j, t]                             # from p_decomposition_constraint_Pij
               ) == 0
    m.dLagr_dPij = Constraint(m.l, m.t, rule=dLagr_dPij_rule)

    # Derivative w.r.t. Qij
    def dLagr_dQij_rule(m, i, j, t):
        return (-m.lambda_BPR[i, t]                                          # from BALANCO_POT_REATIVA (saindo)
                + m.lambda_BPR[j, t]                                         # from BALANCO_POT_REATIVA (entrando)
                - m.lambda_FL[i, j, t] * 2 * DLINES.loc[(i, j), 'X (pu)']     # from FLUXO_LINHAS
                + m.lambda_p_decomp_Qij[i, j, t]                             # from p_decomposition_constraint_Qij
               ) == 0
    m.dLagr_dQij = Constraint(m.l, m.t, rule=dLagr_dQij_rule)

    # Derivative w.r.t. Isq
    def dLagr_dIsq_rule(m, i, j, t):
        return (m.mu_Isq_UP[i, j, t]                                         # from Isq_BOUND_UPPER
                - m.mu_Isq_LOW[i, j, t]                                      # from Isq_BOUND_LOWER
                - m.lambda_BPA[i, t] * DLINES.loc[(i, j), 'R (pu)']           # from BALANCO_POT_ATIVA
                - m.lambda_BPR[i, t] * DLINES.loc[(i, j), 'X (pu)']           # from BALANCO_POT_REATIVA
                - m.lambda_FL[i, j, t] * DLINES.loc[(i, j), 'Z (pu)']**2      # from FLUXO_LINHAS
                - m.mu_FL2[i, j, t]                                          # from FLUXO_LINHAS2
               ) == 0
    m.dLagr_dIsq = Constraint(m.l, m.t, rule=dLagr_dIsq_rule)

    # Derivative w.r.t. Vsq
    def dLagr_dVsq_rule(m, i, t):
        return (m.mu_Vsq_UP[i, t] - m.mu_Vsq_LOW[i, t]
                + (m.lambda_Vsq_SLACK[t] if i == m.slack else 0)
                + sum(m.lambda_FL[j, i, t] for j in m.i if (j, i) in m.l)
                - sum(m.lambda_FL[i, j, t] for j in m.i if (i, j) in m.l)
               ) == 0
    m.dLagr_dVsq = Constraint(m.i, m.t, rule=dLagr_dVsq_rule)

    # Derivative w.r.t. Qg
    def dLagr_dQg_rule(m, g, t):
        return (m.lambda_BPR[g, t]
                + m.mu_Qg_UP[g, t] - m.mu_Qg_LOW[g, t]
                - sum(m.mu_pw_tan_Qg[g, t, k] * m.slope_Qg[g, t, k] for k in m.Y_set_Qg)
               ) == 0
    m.dLagr_dQg = Constraint(m.GB, m.t, rule=dLagr_dQg_rule)

    # Derivative w.r.t. Qg_sq_piecewised
    def dLagr_dQg_sq_rule(m, g, t):
        return (m.mu_GD_UP[g, t]
                - sum(m.mu_pw_tan_Qg[g, t, k] for k in m.Y_set_Qg)
               ) == 0
    m.dLagr_dQg_sq = Constraint(m.GB, m.t, rule=dLagr_dQg_sq_rule)

    # Derivative w.r.t. Ppv
    def dLagr_dPpv_rule(m, r, t):
        return (m.lambda_BPA[r, t]
                + m.mu_Ppv_UP[r, t] - m.mu_Ppv_LOW[r, t]
               ) == 0
    m.dLagr_dPpv = Constraint(m.NGB, m.t, rule=dLagr_dPpv_rule)

    # Derivative w.r.t. Pbess
    def dLagr_dPbess_rule(m, e, t):
        return (m.lambda_BPA[e, t]
                + m.mu_Pbess_UP[e, t] - m.mu_Pbess_LOW[e, t]
                - (m.lambda_SOC[e, t + 1] if t < 24 else 0)
               ) == 0
    m.dLagr_dPbess = Constraint(ess.index, m.t, rule=dLagr_dPbess_rule)

    # Derivative w.r.t. SOC
    def dLagr_dSOC_rule(m, e, t):
        return (m.mu_SOC_UP[e, t] - m.mu_SOC_LOW[e, t]
                + m.lambda_SOC[e, t]
                - (m.lambda_SOC[e, t + 1] if t < 24 else 0)
               ) == 0
    m.dLagr_dSOC = Constraint(ess.index, m.t, rule=dLagr_dSOC_rule)

    # --- Derivatives for auxiliary variables from signed PWL approximation ---
    # For Pij
    def dLagr_dP_plus_Pij_rule(m, i, j, t):
        return m.lambda_p_decomp_Pij[i, j, t] + m.lambda_abs_p_Pij[i, j, t] == 0
    m.dLagr_dP_plus_Pij = Constraint(m.l, m.t, rule=dLagr_dP_plus_Pij_rule)

    def dLagr_dP_minus_Pij_rule(m, i, j, t):
        return -m.lambda_p_decomp_Pij[i, j, t] + m.lambda_abs_p_Pij[i, j, t] == 0
    m.dLagr_dP_minus_Pij = Constraint(m.l, m.t, rule=dLagr_dP_minus_Pij_rule)

    def dLagr_dabs_P_Pij_rule(m, i, j, t):
        return (-m.lambda_abs_p_Pij[i, j, t]
                - sum(m.mu_pw_tan_Pij[i, j, t, k] * m.slope_Pij_abs[i, j, t, k] for k in m.Y_set_Pij_abs)
               ) == 0
    m.dLagr_dabs_P_Pij = Constraint(m.l, m.t, rule=dLagr_dabs_P_Pij_rule)

    def dLagr_dPij_sq_signed_rule(m, i, j, t):
        return (m.mu_FL2[i, j, t]
                - sum(m.mu_pw_tan_Pij[i, j, t, k] for k in m.Y_set_Pij_abs)
               ) == 0
    m.dLagr_dPij_sq_signed = Constraint(m.l, m.t, rule=dLagr_dPij_sq_signed_rule)

    # For Qij
    def dLagr_dP_plus_Qij_rule(m, i, j, t):
        return m.lambda_p_decomp_Qij[i, j, t] + m.lambda_abs_p_Qij[i, j, t] == 0
    m.dLagr_dP_plus_Qij = Constraint(m.l, m.t, rule=dLagr_dP_plus_Qij_rule)

    def dLagr_dP_minus_Qij_rule(m, i, j, t):
        return -m.lambda_p_decomp_Qij[i, j, t] + m.lambda_abs_p_Qij[i, j, t] == 0
    m.dLagr_dP_minus_Qij = Constraint(m.l, m.t, rule=dLagr_dP_minus_Qij_rule)

    def dLagr_dabs_P_Qij_rule(m, i, j, t):
        return (-m.lambda_abs_p_Qij[i, j, t]
                - sum(m.mu_pw_tan_Qij[i, j, t, k] * m.slope_Qij_abs[i, j, t, k] for k in m.Y_set_Qij_abs)
               ) == 0
    m.dLagr_dabs_P_Qij = Constraint(m.l, m.t, rule=dLagr_dabs_P_Qij_rule)

    def dLagr_dQij_sq_signed_rule(m, i, j, t):
        return (m.mu_FL2[i, j, t]
                - sum(m.mu_pw_tan_Qij[i, j, t, k] for k in m.Y_set_Qij_abs)
               ) == 0
    m.dLagr_dQij_sq_signed = Constraint(m.l, m.t, rule=dLagr_dQij_sq_signed_rule)

    #################################################################
    # 3. COMPLEMENTARY SLACKNESS & DUAL FEASIBILITY
    #################################################################
    # For each inequality constraint g(x) <= 0 with dual mu:
    # 1. Dual feasibility: mu >= 0 (already handled by NonNegativeReals domain)
    # 2. Complementarity: -g(x) <= M*(1-z) AND mu <= M*z, where z is binary.

    # --- PG_BOUNDS_UPPER ---
    m.z_Pg_UP = Var(m.GB, m.t, within=Binary)
    def cs_Pg_UP1_rule(m, g, t):
        expr = m.GenD[g, 'P_max (pu)'] - m.Pg[g, t]
        return expr <= m.M * (1 - m.z_Pg_UP[g, t])
    m.cs_Pg_UP1 = Constraint(m.GB, m.t, rule=cs_Pg_UP1_rule)
    def cs_Pg_UP2_rule(m, g, t):
        return m.mu_Pg_UP[g, t] <= m.M * m.z_Pg_UP[g, t]
    m.cs_Pg_UP2 = Constraint(m.GB, m.t, rule=cs_Pg_UP2_rule)

    # --- PG_BOUNDS_LOWER ---
    m.z_Pg_LOW = Var(m.GB, m.t, within=Binary)
    def cs_Pg_LOW1_rule(m, g, t):
        expr = m.Pg[g, t] - m.GenD[g, 'P_min (pu)']
        return expr <= m.M * (1 - m.z_Pg_LOW[g, t])
    m.cs_Pg_LOW1 = Constraint(m.GB, m.t, rule=cs_Pg_LOW1_rule)
    def cs_Pg_LOW2_rule(m, g, t):
        return m.mu_Pg_LOW[g, t] <= m.M * m.z_Pg_LOW[g, t]
    m.cs_Pg_LOW2 = Constraint(m.GB, m.t, rule=cs_Pg_LOW2_rule)

    # --- QG_BOUNDS_UPPER ---
    m.z_Qg_UP = Var(m.GB, m.t, within=Binary)
    def cs_Qg_UP1_rule(m, g, t):
        expr = m.GenD[g, 'P_max (pu)'] - m.Qg[g, t]
        return expr <= m.M * (1 - m.z_Qg_UP[g, t])
    m.cs_Qg_UP1 = Constraint(m.GB, m.t, rule=cs_Qg_UP1_rule)
    def cs_Qg_UP2_rule(m, g, t):
        return m.mu_Qg_UP[g, t] <= m.M * m.z_Qg_UP[g, t]
    m.cs_Qg_UP2 = Constraint(m.GB, m.t, rule=cs_Qg_UP2_rule)

    # --- QG_BOUNDS_LOWER ---
    m.z_Qg_LOW = Var(m.GB, m.t, within=Binary)
    def cs_Qg_LOW1_rule(m, g, t):
        expr = m.Qg[g, t] + m.GenD[g, 'P_max (pu)']
        return expr <= m.M * (1 - m.z_Qg_LOW[g, t])
    m.cs_Qg_LOW1 = Constraint(m.GB, m.t, rule=cs_Qg_LOW1_rule)
    def cs_Qg_LOW2_rule(m, g, t):
        return m.mu_Qg_LOW[g, t] <= m.M * m.z_Qg_LOW[g, t]
    m.cs_Qg_LOW2 = Constraint(m.GB, m.t, rule=cs_Qg_LOW2_rule)

    # --- Isq_BOUND_UPPER ---
    m.z_Isq_UP = Var(m.l, m.t, within=Binary)
    def cs_Isq_UP1_rule(m, i, j, t):
        expr = DLINES.loc[(i, j), 'Line Ampacity (pu)']**2 - m.Isq[i, j, t]
        return expr <= m.M * (1 - m.z_Isq_UP[i, j, t])
    m.cs_Isq_UP1 = Constraint(m.l, m.t, rule=cs_Isq_UP1_rule)
    def cs_Isq_UP2_rule(m, i, j, t):
        return m.mu_Isq_UP[i, j, t] <= m.M * m.z_Isq_UP[i, j, t]
    m.cs_Isq_UP2 = Constraint(m.l, m.t, rule=cs_Isq_UP2_rule)

    # --- Isq_BOUND_LOWER ---
    m.z_Isq_LOW = Var(m.l, m.t, within=Binary)
    def cs_Isq_LOW1_rule(m, i, j, t):
        expr = m.Isq[i, j, t]
        return expr <= m.M * (1 - m.z_Isq_LOW[i, j, t])
    m.cs_Isq_LOW1 = Constraint(m.l, m.t, rule=cs_Isq_LOW1_rule)
    def cs_Isq_LOW2_rule(m, i, j, t):
        return m.mu_Isq_LOW[i, j, t] <= m.M * m.z_Isq_LOW[i, j, t]
    m.cs_Isq_LOW2 = Constraint(m.l, m.t, rule=cs_Isq_LOW2_rule)

    # --- Vsq_BOUND_UPPER ---
    m.z_Vsq_UP = Var(m.i, m.t, within=Binary)
    def cs_Vsq_UP1_rule(m, i, t):
        expr = 1.05**2 - m.Vsq[i, t]
        return expr <= m.M * (1 - m.z_Vsq_UP[i, t])
    m.cs_Vsq_UP1 = Constraint(m.i, m.t, rule=cs_Vsq_UP1_rule)
    def cs_Vsq_UP2_rule(m, i, t):
        return m.mu_Vsq_UP[i, t] <= m.M * m.z_Vsq_UP[i, t]
    m.cs_Vsq_UP2 = Constraint(m.i, m.t, rule=cs_Vsq_UP2_rule)

    # --- Vsq_BOUND_LOWER ---
    m.z_Vsq_LOW = Var(m.i, m.t, within=Binary)
    def cs_Vsq_LOW1_rule(m, i, t):
        expr = m.Vsq[i, t] - 0.95**2
        return expr <= m.M * (1 - m.z_Vsq_LOW[i, t])
    m.cs_Vsq_LOW1 = Constraint(m.i, m.t, rule=cs_Vsq_LOW1_rule)
    def cs_Vsq_LOW2_rule(m, i, t):
        return m.mu_Vsq_LOW[i, t] <= m.M * m.z_Vsq_LOW[i, t]
    m.cs_Vsq_LOW2 = Constraint(m.i, m.t, rule=cs_Vsq_LOW2_rule)

    # (This is a simplified example. A full implementation would require adding
    # binary variables and the two Big-M constraints for EVERY inequality constraint
    # in the model, which is a substantial number.)
    
    
    '''
    
    m.ObjInferior.deactivate()  # Desativa o objetivo do nível inferior
    
    m.OBJ_SUPERIOR = Objective(
    expr=sum(m.lambda_BALANCO_POT_ATIVA[i, t] * m.Ppv[i, t] for i in m.NGB for t in m.t),
    sense=maximize, 
    doc='Objective function for the upper level of the MPEC',
    name='OF_superior'
    )'''


    m.SYMBOL_MAP = {
    # Variáveis
    'Pg': 'P_g',
    'Pd': 'P_d',
    'Pij': 'P_{ij}',
    'Qij': 'Q_{ij}',
    'I': 'I_{ij}',
    'V': 'V_i',
    'Qg': 'Q_g',
    'Ppv': 'P_{pv}',
    'Pbess': 'P_{bess}',
    'SOC': 'SOC_i',
    'GenCost': 'C_g',
    # Parâmetros
    'SBASE': 'S_{base}',
    'VBASE': 'V_{base}',
    'IBASE': 'I_{base}',
    'SOC0': 'SOC_0',
    'GenD': 'G_d',
    'GenND': 'G_{nd}',
    'Load': 'L_i',
    # Conjuntos
    'i': 'N_i',
    't': 'T_t',
    # Restrições
    'I_BOUND_UPPER': 'I_{ij}^{max}',
    'I_BOUND_LOWER': 'I_{ij}^{min}',
    'V_BOUND_UPPER': 'V_i^{max}',
    'V_BOUND_LOWER': 'V_i^{min}',
    'V_SLACK': 'V_{slack}',
    'PG_BOUNDS_UPPER': 'P_g^{max}',
    'PG_BOUNDS_LOWER': 'P_g^{min}',
    'QG_BOUNDS_UPPER': 'Q_g^{max}',
    'QG_BOUNDS_LOWER': 'Q_g^{min}',
    'PPV_BOUNDS_UPPER': 'P_{pv}^{max}',
    'PPV_BOUNDS_LOWER': 'P_{pv}^{min}',
    'PBESS_BOUNDS_UPPER': 'P_{bess}^{max}',
    'PBESS_BOUNDS_LOWER': 'P_{bess}^{min}',
    'SOC_BOUNDS_UPPER': 'SOC_i^{max}',
    'SOC_BOUNDS_LOWER': 'SOC_i^{min}',
    'BALANCO_POT_ATIVA': 'Balanço de Potência Ativa',
    'BALANCO_POT_REATIVA': 'Balanço de Potência Reativa',
    'FLUXO_LINHAS': 'Fluxo de Linhas',
    'FLUXO_LINHAS2': 'Limite de Potência Aparente',
    'GD_UPPER_BOUND': 'Limite de Potência Aparente dos Geradores Despacháveis',
    'GND_BOUNDS_UPPER': 'Limite de Potência Aparente dos Geradores Não Despacháveis',
    'STATE_OF_CHARGE': 'Estado de Carga do ESS',
    'ObjInferior': 'Custo de Geração'
    }
    
    #if analyze_model:
    #    analyze_model_gurobi(m)
    
    return m
