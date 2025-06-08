#-----------------------------------------------------------------------
#  DECLARE SET OF BUS
#-----------------------------------------------------------------------
param VBASE  := 24.90;   # base voltage magnitude (kV) 
param SBASE  := 100 ;   # base apparent power (MVA)
param IBASE := (SBASE*1e3/VBASE);


param SE_Capacity;
let SE_Capacity := 50;

set N;                                  # set of buses
set G_T;                                # set of transmission-connected generators
set G_D;                                # set of distribution-connected generators
set LS;									# set of load shifting nodes
set ESS;								# set of energy storage system nodes
set RES_T;								# set of transmission-connected renewable energy sources
set RES_D;								# set of distribution-connected renewable energy sources
set L within N cross N; 				# set of branches
set T ordered;
set Trans_Lines;
set Trans_Nodes;
set LDA;
set WS within 1..3 cross Trans_Nodes cross (G_T union LDA);		# set of wholesale market participants
set S;


#-----------------------------------------------------------------------
# DECLARE ITEMS TO BE READ FROM THE .bus DATA FILE.
#-----------------------------------------------------------------------

param NBUS {N} symbolic; # bus name
param TBUS {N};          # bus type 
param Vmax {N};          # maximum voltage magnitude (kV) 
param Vmin {N};          # minimum voltage magnitude (kV)
param Vnom {N};
param Ve   {N};
param Pgmx {N};          # maximum real power (kW) 
param Pgmn {N};          # minimum real power (kW)
param R    {L};          # branch resistance (ohm)
param X    {L};          # branch reactance (ohm)
param Z    {L};          # branch impedance (ohm)
param Imax {L};          # branch maximum current (kA)

param load1{S,T};
param load2{S,T};
param wind1{S,T};
param wind2{S,T};
param solar{S,T};
param probability{S,T};

#load_shape_data.csv	
param ind1{T};
param ind2{T};
param ind3{T};

param res1{T};
param res4{T};
param res5{T};

param off1{T};
param off4{T};

param L2{T};
param L3{T};
param L4{T};
param L6{T};
param L8{T};

param PCC{LDA};

#transmission_node_data.csv
param Trans_Load{Trans_Nodes, T, S};
param Trans_Shift_Max{Trans_Nodes};
param Trans_Shift_Min{Trans_Nodes};
param Trans_P_ESS_min{Trans_Nodes};
param Trans_P_ESS_max{Trans_Nodes};
param Trans_SOC_max{Trans_Nodes};
param Trans_SOC_min{Trans_Nodes};
param Trans_SOC_initial{Trans_Nodes};
param Trans_Type{Trans_Nodes} symbolic;
param Trans_Inst_Cap{Trans_Nodes};


param Trans_Wind_Inj{RES_T, T, S};
param Trans_Solar_Inj{RES_T, T, S};
param RES_type_t{RES_T} symbolic;
param RES_Node_t{RES_T};
param RES_Cap_t{RES_T};
param Wind_min_t{RES_T};
param Wind_max_t{RES_T};
param Wind_nom_t{RES_T};





param Dist_Load{N, T, S};
param LS_Node{LS};
param Dist_Shift_Max{LS};
param Dist_Shift_Min{LS};
param ESS_Node{ESS};
param Dist_P_ESS_min{ESS};
param Dist_P_ESS_max{ESS};
param Dist_SOC_max{ESS};
param Dist_SOC_min{ESS};
param Dist_SOC_initial{ESS};
param Dist_Type{N} symbolic;
param Dist_Inst_Cap{N};



param Dist_Wind_Inj{RES_D, T, S};
param Dist_Solar_Inj{RES_D, T, S};
param RES_type_d{RES_D} symbolic;
param RES_Node_d{RES_D};
param RES_Cap_d{RES_D};
param Wind_min_d{RES_D};
param Wind_max_d{RES_D};
param Wind_nom_d{RES_D};



param Trans_Incidencia{Trans_Nodes,Trans_Lines};
param Trans_Reactance{Trans_Lines};
param Trans_Capacity{Trans_Lines};
param Trans_Status{Trans_Lines};
param From{Trans_Lines};
param To{Trans_Lines};


param Dist_X{L};
param Dist_R{L};
param Dist_Z{L};
param Dist_Status{L};



param Pmax_t{G_T};
param Pmin_t{G_T};
param a_t{G_T};
param b_t{G_T};
param c_t{G_T};
param G_Node{G_T};
param G_Owner_t{G_T};
param Carbon_Cost_t{G_T};

param Pmax_d{G_D};
param Pmin_d{G_D};
param a_d{G_D};
param b_d{G_D};
param c_d{G_D};
param G_Owner_d{G_D};
param G_LDA{G_D};
param G_LDA_Node{G_D};
param Carbon_Cost_d{G_D};

param Carbon_Limit_Trans{T};
param Carbon_Limit_Dist{T};
param Emission_Weighted_Average;
#---- DECLARE VARIABLES, WITH UPPER AND LOWER BOUNDS.

# Primal variables
var VM {N, T, S};
var I  {L, T, S};

var Trans_Shift{Trans_Nodes, T, S};
var Trans_P_ESS{Trans_Nodes, T, S};
var Trans_Flow{Trans_Lines, T, S};
var Trans_Theta{Trans_Nodes, T, S};


var Dist_Shift{LS, T, S};
var Dist_P_ESS{ESS, T, S};
var Dist_SOC{ESS,T, S};
var Dist_Flow {L, T, S};
var Dist_Pfm{L, T, S};
var Dist_Pto{L, T, S};
var Dist_Qfm{L, T, S};
var Dist_Qto{L, T, S};
var P{L, T, S};
var Q{L, T, S};


var P_thermal_dist{G_D, T, S} >= 0;
var Q_thermal_dist{G_D, T, S};
var P_thermal_trans{G_T, T, S} >= 0;
var P_DSO{T, S};
var Q_DSO{T,S};


var Footprint_trans{G_T,T,S};
var Footprint_dist{G_D,T,S};
var Carbon_T{G_T,T,S};
var Carbon_SE{T,S} ;

#Dual Variables
var lambda{Trans_Nodes, T, S};
var omega_U{Trans_Lines, T, S};
var omega_L{Trans_Lines, T, S};	

var kappa_U{T, S};
var kappa_L{T, S};
var alpha_U{G_T, T, S};
var alpha_L{G_T, T, S};
var ramp_up{G_T, T, S};
var ramp_low{G_T, T, S};

var mu{Trans_Lines, T, S};
var phi{G_T, T, S};
var psi{T, S};
var psi_L{T, S};
var epsilon{G_T, T, S};
var Xi{T, S};
var eta_U{T, S};
var eta_L{T, S};




var Bid{T,S} <= 500;
var Carbon_Price{T, S} <= 5;


var DSO_Revenue{T, S};
var GenCosts_dist{G_D, T, S} >= 0;
var GenCosts_trans{G_T, T, S} >= 0;

# eq 1
minimize DSO_Costs: sum{s in S, t in T} (1 * (
					(Bid[t,s] * P_DSO[t,s]) 
					+ sum{g in G_D}(GenCosts_dist[g,t,s])
					+(Carbon_Price[t,s] * Carbon_SE[t,s])
					));

s.t. Calculate_GenCosts_Dist {t in T, g in G_D, s in S}: 
	GenCosts_dist[g,t,s] = a_d[g]*P_thermal_dist[g,t,s]^2 + b_d[g]*P_thermal_dist[g,t,s] + c_d[g];
	
s.t. Calculate_GenCosts_Trans {t in T, g in G_T, s in S}: 
	GenCosts_trans[g,t,s] = a_t[g]*P_thermal_trans[g,t,s]^2 + b_t[g]*P_thermal_trans[g,t,s] + c_t[g];


#s.t. Calculate_DSO_Revenue{t in T}: DSO_Revenue[t] = lambda[5,t] * P_DSO[t] 
#								+ sum{g in G_T, n in LDA} (Carbon_Price[t] * Carbon_SE[n,g,t]);

s.t. ACTIVE_POWER_BALANCE{n in N, t in T, s in S}: #eq
#soma dos geradores nos nós G_LDA_Node[g] para os g no set G_D quando g é igual a n

	sum{g in G_D:G_LDA_Node[g] == n}(P_thermal_dist[g,t,s]) # pot geradores

	-sum{(n,m) in L}(P[n,m,t,s] # pot saindo
	+R[n,m]*I[n,m,t,s])

	+sum{(l,n) in L}(P[l,n,t,s]) # pot entrando

	+(if 1 == n then P_DSO[t,s]) # pot intercambio

	-sum{i in LS:LS_Node[i] == n} (Dist_Shift[i,t,s]) # load shifting 

	+sum{i in ESS:ESS_Node[i] == n} (Dist_P_ESS[i,t,s]) # pot ESS
	
	= Dist_Load[n,t,s];


#s.t. REACTIVE_POWER_BALANCE{n in N, t in T, s in S}: 
#	sum{g in G_D:G_LDA_Node[g] == n}(Q_thermal_dist[g,t,s])
#	-sum{(n,m) in L}(Q[n,m,t,s] + X[n,m]*I[n,m,t,s])
#	+sum{(l,n) in L}(Q[l,n,t,s]) 
#	+(if 1 == n then Q_DSO[t,s])
#	-sum{i in LS:LS_Node[i] == n} (Dist_Shift[i,t,s] * tan(acos(.8))) 
#	= Dist_Load[n,t,s] * tan(acos(.8));
	
	
#s.t. new {(i,j) in L, t in T, s in S}:
#	VM[i,t,s] - VM[j,t,s] = 2*(R[i,j]*P[i,j,t,s] + X[i,j]*Q[i,j,t,s]) + Z[i,j]^2*I[i,j,t,s];
	
#s.t. new2{(i,j) in L, t in T, s in S}:
#	P[i,j,t,s]^2 + Q[i,j,t,s]^2 = I[i,j,t,s]*VM[j,t,s];


s.t. LOSS_REAL_POWER {(i,j) in L, t in T, s in S}: 
	(Dist_Pfm[i,j,t,s] + Dist_Pto[i,j,t,s]) = R[i,j] *I[i,j,t,s]^2;

s.t. REAL_POWER_FLOW {(i,j) in L, t in T, s in S}:
	(Dist_Pfm[i,j,t,s] - Dist_Pto[i,j,t,s]) = R[i,j]/Z[i,j]^2 * (VM[i,t,s]^2 - VM[j,t,s]^2);

#s.t. LOSS_IMAG_POWER {(i,j) in L, t in T, s in S}: 
#	(Dist_Qfm[i,j,t,s] + Dist_Qto[i,j,t,s]) = X[i,j] *I[i,j,t,s]^2;

#s.t. IMAG_POWER_FLOW {(i,j) in L, t in T, s in S}:
#	(Dist_Qfm[i,j,t,s] - Dist_Qto[i,j,t,s]) = X[i,j]/Z[i,j]^2 * (VM[i,t,s]^2 - VM[j,t,s]^2);


#---- CURRENT FLOW EQUATION
s.t. CURRENT_FLOW {(i,j) in L, t in T, s in S}:
	I[i,j,t,s] = (VM[i,t,s]-VM[j,t,s])/Z[i,j];


	
#---- LINE RATING LIMITS
s.t. MINIMUM_CURRENT_FLOW {(n,m) in L, t in T, s in S}: 
	I[n,m,t,s] >= -Imax[n,m];
  
s.t. MAXIMUM_CURRENT_FLOW {(n,m) in L, t in T, s in S}: 
    I[n,m,t,s] <= Imax[n,m];#^2;												

#s.t. MINIMUM_ACTIVE_FLOW {(n,m) in L, t in T, s in S}: 
#	P[n,m,t,s] >= -Imax[n,m]*Vmax[m];
  
#s.t. MAXIMUM_ACTIVE_FLOW {(n,m) in L, t in T, s in S}: 
#    P[n,m,t,s] <= Imax[n,m]*Vmax[m];	
    
#s.t. MINIMUM_REACTIVE_FLOW {(n,m) in L, t in T, s in S}: 
#	Q[n,m,t,s] >= -Imax[n,m]*Vmax[m];
  
#s.t. MAXIMUM_REACTIVE_FLOW {(n,m) in L, t in T, s in S}: 
#    Q[n,m,t,s] <= Imax[n,m]*Vmax[m];	

#---- VOLTAGE MAGNITUDE LIMITS
  
s.t. MAXIMUM_VOLTAGE_MAGNITUDE {t in T,n in N, s in S}:
    VM[n,t,s] <= Vmax[n];#^2;													
  
s.t. MINIMUM_VOLTAGE_MAGNITUDE {t in T,n in N, s in S}:
  	VM[n,t,s] >= Vmin[n];#^2;													
  
#---- ACTIVE GENERATION LIMITS

#s.t. MAXIMUM_DG_APARENT_POWER {t in T, g in G_D, s in S}:
#    P_thermal_dist[g,t,s]^2 + Q_thermal_dist[g,t,s]^2 <=   Pmax_d[g]^2;											

s.t. MAXIMUM_DG_ACTIVE_POWER {t in T, g in G_D, s in S}:
    P_thermal_dist[g,t,s] <=   Pmax_d[g];											

s.t. MINIMUM_DG_ACTIVE_POWER {t in T, g in G_D, s in S}:
    P_thermal_dist[g,t,s] >=   Pmin_d[g];											

#s.t. MAXIMUM_DG_REACTIVE_POWER {t in T, g in G_D, s in S}:
#    Q_thermal_dist[g,t,s] <=   Pmax_d[g];											

#s.t. MINIMUM_DG_REACTIVE_POWER {t in T, g in G_D, s in S}:
#    Q_thermal_dist[g,t,s] >=  -Pmax_d[g];	

#---- LOAD SHIFTING LIMITS
s.t. MAXIMUM_LS_LIMIT {t in T, n in LS, s in S}:
	Dist_Shift[n,t,s] <= Dist_Shift_Max[n] *Dist_Load[LS_Node[n],t,s]; 			
	
s.t. MINIMUM_LS_LIMIT {t in T, n in LS, s in S}:
	-Dist_Shift[n,t,s] <= -Dist_Shift_Min[n] *Dist_Load[LS_Node[n],t,s]; 		

s.t. ZERO_SUM_LS {n in LS, s in S}:
	sum{t in T} Dist_Shift[n,t,s] = 0; 										
	
	
#---- ENERGY STORAGE SYSTEMS LIMITS
s.t. Dist_ESS_Limits  {n in ESS, t in T, s in S}: Dist_P_ESS[n,t,s] >= Dist_P_ESS_min[n]; 					
s.t. Dist_ESS_Limits2 {n in ESS, t in T, s in S}: Dist_P_ESS[n,t,s] <= Dist_P_ESS_max[n]; 					

s.t. Dispatch_Dist_SOCLimits  {n in ESS, t in T, s in S}: Dist_SOC[n,t,s] >= Dist_SOC_min[n];  			
s.t. Dispatch_Dist_SOCLimits2 {n in ESS, t in T, s in S}: Dist_SOC[n,t,s] <= Dist_SOC_max[n]; 				

s.t. Dispatch_Dist_CalculateSOC {n in ESS, t in T, s in S}: (if t = first(T) then Dist_SOC_initial[n] 
else Dist_SOC[n,prev(t),s]) - Dist_P_ESS[n,t,s] = Dist_SOC[n,t,s]; 										

s.t. Dispatch_Dist_SameSOC {n in ESS, s in S}: Dist_SOC_initial[n] = Dist_SOC[n, last(T),s]; 				




#---- CARBON RELATED CONSTRAINTS
s.t. Calculate_Footprint_Dist {t in T, g in G_D, s in S}:
	Footprint_dist[g,t,s] = Carbon_Cost_d[g] * P_thermal_dist[g,t,s] ;
	
	
#subject to Dist_Carbon_Trade_Limit1 {t in T, s in S}:
# 	-sum{g in G_D} (Footprint_dist[g,t,s]) + Carbon_SE[t,s] <= 0;


s.t. Dist_Carbon_Trade_Limit2 {t in T, s in S}:
	sum{g in G_D} (Footprint_dist[g,t,s]) - Carbon_SE[t,s] <= Carbon_Limit_Dist[t] ;	
	










##									TRANSMISSION PRIMAL CONSTRAINTS
# MIN(GenCosts - P_DSO*Bid)
s.t. Calculate_Trans_PowerFlow {l in Trans_Lines, t in T, s in S}: 
	Trans_Flow[l,t,s] - Trans_Status[l] * 1/Trans_Reactance[l] * 
	sum{n in Trans_Nodes}(Trans_Incidencia[n,l]*Trans_Theta[n,t,s]) = 0;					#mu


s.t. Trans_Power_Balance{n in Trans_Nodes, t in T, s in S}: 
	sum{g in G_T : G_Node[g] == n}(P_thermal_trans[g,t,s])
	-sum{l in Trans_Lines}(Trans_Incidencia[n,l]*Trans_Flow[l,t,s])
	-(if 5 == n then P_DSO[t,s]) 
	= Trans_Load[n,t,s];																	#lambda


s.t. Trans_Flow_Upper{l in Trans_Lines, t in T, s in S}:  Trans_Flow[l,t,s] <= Trans_Capacity[l]; #omega_U
s.t. Trans_Flow_Lower{l in Trans_Lines, t in T, s in S}: -Trans_Flow[l,t,s] <= Trans_Capacity[l]; #omega_L

s.t. Trans_Gen_Upper{g in G_T, t in T, s in S}:  P_thermal_trans[g,t,s] <=  Pmax_t[g]; #alpha_U
s.t. Trans_Gen_Lower{g in G_T, t in T, s in S}: -P_thermal_trans[g,t,s] <= -Pmin_t[g]; #alpha_L

s.t. Trans_Gen_Upper_Ramp{g in G_T, t in T, s in S: t<>first(T)}:  P_thermal_trans[g,t-1,s] - P_thermal_trans[g,t,s] <=  0.25 * Pmax_t[g]; #alpha_U
s.t. Trans_Gen_Lower_Ramp{g in G_T, t in T, s in S: t<>first(T)}: -P_thermal_trans[g,t-1,s] + P_thermal_trans[g,t,s] <=  0.25 * Pmax_t[g]; #alpha_L

s.t. MAXIMUM_SE_LIMIT {t in T, s in S}:    P_DSO[t,s] <=   SE_Capacity; #kappa_U											
s.t. MINIMUM_SE_LIMIT {t in T, s in S}:   -P_DSO[t,s] <=  SE_Capacity;  #kappa_L





subject to Calculate_Footprint_Trans {g in G_T, t in T, s in S}:
	Footprint_trans[g,t,s] = Carbon_Cost_t[g] * P_thermal_trans[g,t,s];		#phi
	 
#subject to Generator_Footprint_Limit {g in G_T, t in T, s in S}:
# 	-Footprint_trans[g,t,s] - Carbon_T[g,t,s] <= 0;							#epsilon

s.t. Carbon_Balance_Trans {t in T, s in S}:
	sum{g in G_T} (Carbon_T[g,t,s]) - Carbon_SE[t,s] = 0;					#psi

s.t. Carbon_Limits_Trans {t in T, s in S}:									
	sum{g in G_T} (Footprint_trans[g,t,s] + Carbon_T[g,t,s]) <= Carbon_Limit_Trans[t] ;		#Xi

#s.t. MAXIMUM_CARBON_LIMIT {t in T, s in S}:    Carbon_SE[t,s] <=  Carbon_Limit_Trans[t]; #eta_U											
#s.t. MINIMUM_CARBON_LIMIT {t in T, s in S}:   -Carbon_SE[t,s] <=  Carbon_Limit_Trans[t]; #eta_L

##									TRANSMISSION DUAL CONSTRAINTS
s.t. Deriv_Potencia_with_Ramp {g in G_T, t in T, s in S}:
	(2*a_t[g] * P_thermal_trans[g,t,s] + b_t[g]) - lambda[G_Node[g],t,s] 
	-phi[g,t,s]*Carbon_Cost_t[g] 
	- alpha_L[g,t,s] + alpha_U[g,t,s] 
	+(if t<>last(T) then (ramp_up[g,t+1,s] - ramp_low[g,t+1,s])*0.25 * Pmax_t[g])
	+(if t<>first(T) then (-ramp_up[g,t,s] + ramp_low[g,t,s])*0.25 * Pmax_t[g])
	= 0;
	
s.t. Deriv_Potencia_without_Ramp {g in G_T, t in T, s in S}:
	(2*a_t[g] * P_thermal_trans[g,t,s] + b_t[g]) - lambda[G_Node[g],t,s] 
	-phi[g,t,s]*Carbon_Cost_t[g] 
	- alpha_L[g,t,s] + alpha_U[g,t,s] 
	= 0;

s.t. Deriv_Teta{n in Trans_Nodes, t in T, s in S}: -sum{l in Trans_Lines}(
	-Trans_Status[l] * (1/Trans_Reactance[l]) * Trans_Incidencia[n,l]*mu[l,t,s]) = 0;

s.t. Deriv_Fluxo{l in Trans_Lines, t in T, s in S}: mu[l,t,s] - omega_L[l,t,s] + omega_U[l,t,s] 
	-sum{n in Trans_Nodes}(Trans_Incidencia[n,l] * lambda[n,t,s]) = 0;

s.t. Deriv_PDSO{t in T, s in S}:  
	-Bid[t,s] + lambda[5,t,s] + kappa_U[t,s] - kappa_L[t,s] = 0;

s.t. Deriv_Footprint{g in G_T, t in T, s in S}:
	+phi[g,t,s]  #- epsilon[g,t,s] 
	+ Xi[t,s] = 0;

s.t. Deriv_Carbon_T {g in G_T, t in T, s in S}:
	#-epsilon[g,t,s] 
	+ psi[t,s] + Xi[t,s] = 0;
	
s.t. Deriv_Carbon_SE {t in T, s in S}:
	 -Carbon_Price[t,s] -psi[t,s] 
#	 + eta_U[t,s] - eta_L[t,s] 
	 = 0;


s.t. Omega_Lower{l in Trans_Lines, t in T, s in S}: omega_L[l,t,s] >= 0;
s.t. Omega_Upper{l in Trans_Lines, t in T, s in S}: omega_U[l,t,s] >= 0;

s.t. Alpha_Lower{g in G_T, t in T, s in S}: alpha_L[g,t,s] >= 0;
s.t. Alpha_Upper{g in G_T, t in T, s in S}: alpha_U[g,t,s] >= 0;

s.t. Ramp_Lower{g in G_T, t in T, s in S}: ramp_low[g,t,s] >= 0;
s.t. Ramp_Upper{g in G_T, t in T, s in S}: ramp_up[g,t,s] >= 0;

s.t. Kappa_Lower{t in T, s in S}: kappa_L[t,s] >= 0;
s.t. Kappa_Upper{t in T, s in S}: kappa_U[t,s] >= 0;

#s.t. Epsilon_Limit{g in G_T, t in T, s in S}: epsilon[g,t,s] >=0;
s.t. Xi_Limit{t in T, s in S}: Xi[t,s] >=0;

#s.t. Eta_Lower{t in T, s in S}: eta_L[t,s] >= 0;
#s.t. Eta_Upper{t in T, s in S}: eta_U[t,s] >= 0;

##					LINEAR KKT COMPLEMENTARY SLACKNESS
s.t. slack_NL2 {l in Trans_Lines, t in T, s in S}: omega_U[l,t,s]*(-Trans_Capacity[l] + Trans_Flow[l,t,s]) = 0;
s.t. slack_NL1 {l in Trans_Lines, t in T, s in S}: omega_L[l,t,s]*(-Trans_Flow[l,t,s] - Trans_Capacity[l]) = 0;
s.t. slack_NL4 {g in G_T, t in T, s in S}: alpha_U[g,t,s]*(P_thermal_trans[g,t,s] - Pmax_t[g]) = 0;
s.t. slack_NL3 {g in G_T, t in T, s in S}: alpha_L[g,t,s]*(-P_thermal_trans[g,t,s] + Pmin_t[g]) = 0;
s.t. slack_NL5 {t in T, s in S}: kappa_U[t,s]*(-SE_Capacity + P_DSO[t,s]) = 0;
s.t. slack_NL6 {t in T, s in S}: kappa_L[t,s]*(-P_DSO[t,s] - SE_Capacity) = 0;
s.t. slack_NL7 {t in T, s in S} :Xi[t,s]*(sum{g in G_T} (Footprint_trans[g,t,s] + Carbon_T[g,t,s]) - Carbon_Limit_Trans[t] ) = 0;
#s.t. slack_NL8 {g in G_T, t in T, s in S} :epsilon[g,t,s]*(-Footprint_trans[g,t,s] - Carbon_T[g,t,s]) = 0;
#s.t. slack_NL9 {t in T, s in S}: eta_U[t,s]*(-Carbon_Limit_Trans[t] + Carbon_SE[t,s]) = 0;
#s.t. slack_NL10{t in T, s in S}: eta_L[t,s]*(-Carbon_SE[t,s] - Carbon_Limit_Trans[t]) = 0;


s.t. slack_NL11{g in G_T, t in T, s in S : t<>first(T)}: ramp_up[g,t,s]*(P_thermal_trans[g,t-1,s] - P_thermal_trans[g,t,s] -0.25 * Pmax_t[g]) = 0;
s.t. slack_NL12{g in G_T, t in T, s in S : t<>first(T)}: ramp_low[g,t,s]*(-P_thermal_trans[g,t-1,s] + P_thermal_trans[g,t,s] -0.25 * Pmax_t[g]) = 0;









