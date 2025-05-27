#-----------------------------------------------------------------------
#  DECLARE SET OF BUS
#-----------------------------------------------------------------------
param VBASE  := 24.90;   # base voltage magnitude (kV) 
param SBASE  := 100 ;   # base apparent power (MVA)
param IBASE := (SBASE*1e3/VBASE);
param PWmax_g:= 5;       # Maximun discretization steps for generator costs linearization
param M := 10;

param SE_Capacity;
let SE_Capacity := 50/SBASE;

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
set B := 1..PWmax_g ordered;          # set block offers


#-----------------------------------------------------------------------
# DECLARE ITEMS TO BE READ FROM THE .bus DATA FILE.
#-----------------------------------------------------------------------


param P_block_t{G_T, B};
param Price_block_t{G_T, B};


param P_block_d{G_D, B};
param Price_block_d{G_D, B};




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


param Trans_Load{Trans_Nodes, T};
param Trans_Shift_Max{Trans_Nodes};
param Trans_Shift_Min{Trans_Nodes};
param Trans_P_ESS_min{Trans_Nodes};
param Trans_P_ESS_max{Trans_Nodes};
param Trans_SOC_max{Trans_Nodes};
param Trans_SOC_min{Trans_Nodes};
param Trans_SOC_initial{Trans_Nodes};
param Trans_Type{Trans_Nodes} symbolic;
param Trans_Inst_Cap{Trans_Nodes};


param Trans_Wind_Inj{RES_T, T};
param Trans_Solar_Inj{RES_T, T};
param RES_type_t{RES_T} symbolic;
param RES_Node_t{RES_T};
param RES_Cap_t{RES_T};
param Wind_min_t{RES_T};
param Wind_max_t{RES_T};
param Wind_nom_t{RES_T};





param Dist_Load{N, T};
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



param Dist_Wind_Inj{RES_D, T};
param Dist_Solar_Inj{RES_D, T};
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

param scenario;

#---- DECLARE VARIABLES, WITH UPPER AND LOWER BOUNDS.

# Primal variables
var VM {N, T};
var I  {L, T};

var Trans_Flow{Trans_Lines, T};
var Trans_Theta{Trans_Nodes, T};


var Dist_Shift{LS, T};
var Dist_P_ESS{ESS, T};
var Dist_SOC{ESS,T};
var Dist_Flow {L, T};
var Dist_Pfm{L, T};
var Dist_Pto{L, T};
#var Dist_Qfm{L, T};
#var Dist_Qto{L, T};
var P{L, T};
var Q{L, T};


var P_thermal_dist{G_D, T} >= 0;
var Q_thermal_dist{G_D, T};
var P_renewable_dist{RES_D, T};
var Q_renewable_dist{RES_D, T};
var P_thermal_trans{G_T, T} >= 0;

var Block_Dispatch_t{G_T, B, T};
var Block_Dispatch_d{G_D, B, T};


var P_DSO{T};
var Q_DSO{T};


var Footprint_trans{G_T,T};
var Footprint_dist{G_D,T};
var Carbon_T{G_T,T};
var Carbon_D{G_D,T};
var Carbon_SE{T} ;

#Dual Variables
var lambda{Trans_Nodes, T};
var omega_U{Trans_Lines, T};
var omega_L{Trans_Lines, T};	

var kappa_U{T};
var kappa_L{T};
var alpha_U{G_T, T};
var alpha_L{G_T, T};
var ramp_up{G_T, T};
var ramp_low{G_T, T};
var rho_U{G_T, B, T};
var rho_L{G_T, B, T};
var varphi{G_T,T};

var mu{Trans_Lines, T};
var phi{G_T, T};
var psi{T};
#var psi_L{T};
#var epsilon{G_T, T};
var Xi_U{T};
var Xi_L{T};
#var eta_U{T};
#var eta_L{T};




var z_kappa_U{T} binary;
var z_kappa_L{T} binary;
var z_omega_U{Trans_Lines, T} binary;
var z_omega_L{Trans_Lines, T} binary;
var z_alpha_U{G_T,T} binary;
var z_alpha_L{G_T,T} binary;
var z_ramp_U{G_T,T} binary;
var z_ramp_L{G_T,T} binary;
var z_rho_U{G_T,B,T} binary;
var z_rho_L{G_T,B,T} binary;
var z_Xi_U{T} binary;
var z_Xi_L{T} binary;


var Bid{T};
var Carbon_Price{T};# >= 1e-3;
var DSO_MKT{T};


#var DSO_Revenue{T};
var GenCosts_dist{G_D, T} >= 0;
var GenCosts_trans{G_T, T} >= 0;

minimize DSO_Costs: + sum{t in T} (probability[scenario,t] * (
					+ sum{g in G_D}(GenCosts_dist[g,t])
					+DSO_MKT[t]
#					+ sum{g in G_T}(GenCosts_trans[g,t])
					
					#+(Bid[t] * P_DSO[t])					
					#+(Carbon_Price[t] * Carbon_SE[t])
					));


subj to Calculate_GenCosts_Dist {t in T, g in G_D}: 
	GenCosts_dist[g,t] = sum{b in B}(Block_Dispatch_d[g,b,t] * Price_block_d[g,b]);
	
subj to Calculate_GenCosts_Trans {t in T, g in G_T}: 
	GenCosts_trans[g,t] = sum{b in B}(Block_Dispatch_t[g,b,t] * Price_block_t[g,b]);
	

subj to Linearize_Market_Revenues{t in T}:
	DSO_MKT[t] =
				-(
					sum{n in Trans_Nodes}(lambda[n,t]*(
											+Trans_Load[n,t] * load2[scenario,t]
											-sum{g in RES_T:RES_Node_t[g] == n and RES_type_t[g]=="SOLAR"} (Trans_Solar_Inj[g,t] )
											-sum{g in RES_T:RES_Node_t[g] == n and RES_type_t[g]=="WIND"} (Trans_Wind_Inj[g,t])
										  )) 
					+ sum{l in Trans_Lines}(-Trans_Capacity[l]*(omega_U[l,t]+omega_L[l,t]))
					+ sum{g in G_T}(-alpha_U[g,t]*Pmax_t[g]#-alpha_L[g,t]*Pmin_t[g] 
									)
					+ sum{g in G_T, b in B} (-rho_U[g,b,t] * P_block_t[g,b])
					- Xi_U[t] * Carbon_Limit_Trans[t] 
					- SE_Capacity * (kappa_U[t] + kappa_L[t])
				)
				+sum{g in G_T} (GenCosts_trans[g,t])
	;



#subj to Linearize_Market_Revenues_without_CarbonLimits{t in T}:
#	DSO_MKT[t] =
#				-(
#					sum{n in Trans_Nodes}(-lambda[n,t]*Trans_Load[n,t] * load2[scenario,t]
#											+sum{g in RES_T:RES_Node_t[g] == n and RES_type_d[g]=="SOLAR"} (Trans_Solar_Inj[g,t] ) *lambda[n,t]
#											+sum{g in RES_T:RES_Node_t[g] == n and RES_type_d[g]=="WIND"} (Trans_Wind_Inj[g,t]) *lambda[n,t]
#										  ) 
#					+ sum{l in Trans_Lines}(-Trans_Capacity[l]*(omega_U[l,t]+omega_L[l,t]))
#					+ sum{g in G_T}(-alpha_U[g,t]*Pmax_t[g]#-alpha_L[g,t]*Pmin_t[g] 
#									-  0.25 * Pmax_t[g] *(ramp_up[g,t]+ramp_low[g,t])
#									)
#					- SE_Capacity * (kappa_U[t] + kappa_L[t])
#				)
#				+sum{g in G_T} (GenCosts_trans[g,t])
#				-SE_Capacity * (kappa_U[t] + kappa_L[t])
#	;


#s.t. Calculate_DSO_Revenue{t in T}: DSO_Revenue[t] = lambda[5,t] * P_DSO[t] 
#								+ sum{g in G_T, n in LDA} (Carbon_Price[t] * Carbon_SE[n,g,t]);

s.t. ACTIVE_POWER_BALANCE{n in N, t in T}: 
	sum{g in G_D:G_LDA_Node[g] == n}(P_thermal_dist[g,t])
	
	-sum{(n,m) in L}(P[n,m,t] + R[n,m]*I[n,m,t])
	+sum{(l,n) in L}(P[l,n,t]) 
	
	+(if 1 == n then P_DSO[t])
	-sum{i in LS:LS_Node[i] == n} (Dist_Shift[i,t]) 
	+sum{i in ESS:ESS_Node[i] == n} (Dist_P_ESS[i,t])
	
	+sum{g in RES_D:RES_Node_d[g] == n} (P_renewable_dist[g,t] )
	
	= Dist_Load[n,t] * load1[scenario,t];
	
	
s.t. REACTIVE_POWER_BALANCE{n in N, t in T}: 
	sum{g in G_D:G_LDA_Node[g] == n}(Q_thermal_dist[g,t])
	
	-sum{(n,m) in L}(Q[n,m,t] + X[n,m]*I[n,m,t])
	+sum{(l,n) in L}(Q[l,n,t]) 
	
	+(if 1 == n then Q_DSO[t])
	-sum{i in LS:LS_Node[i] == n} (Dist_Shift[i,t] * tan(acos(0.9))) 
	
	+sum{g in RES_D:RES_Node_d[g] == n} (Q_renewable_dist[g,t] )
	
	= Dist_Load[n,t] * load1[scenario,t] * tan(acos(.9));
	
	
	
s.t. new {(i,j) in L, t in T}:
	VM[i,t] - VM[j,t] = 2*(R[i,j]*P[i,j,t] +X[i,j]*Q[i,j,t]) + Z[i,j]^2*I[i,j,t];
	
s.t. new2{(i,j) in L, t in T}:
	P[i,j,t]^2 + Q[i,j,t]^2 <= I[i,j,t]*VM[j,t];
	
#s.t. new2{(i,j) in L, t in T}:
#	P[i,j,t]^2  <= I[i,j,t]*VM[j,t];


#---- CURRENT FLOW EQUATION
#s.t. CURRENT_FLOW {(i,j) in L, t in T}:
#	I[i,j,t] = (VM[i,t]-VM[j,t])/Z[i,j];


	
#---- LINE RATING LIMITS
s.t. MINIMUM_CURRENT_FLOW {(n,m) in L, t in T}: 
	I[n,m,t] >= 0;#-Imax[n,m];
  
s.t. MAXIMUM_CURRENT_FLOW {(n,m) in L, t in T}: 
    I[n,m,t] <= Imax[n,m]^2;												

#---- VOLTAGE MAGNITUDE LIMITS
  
s.t. MAXIMUM_VOLTAGE_MAGNITUDE {t in T,n in N}:
    VM[n,t] <= Vmax[n]^2;													
  
s.t. MINIMUM_VOLTAGE_MAGNITUDE {t in T,n in N}:
  	VM[n,t] >= Vmin[n]^2;													
  
#---- ACTIVE GENERATION LIMITS

s.t. Calculate_Block_Upper {g in G_D, b in B, t in T}:
	Block_Dispatch_d[g,b,t] <= P_block_d[g,b];
	
s.t. Calculate_Block_Lower {g in G_D, b in B, t in T}:
	Block_Dispatch_d[g,b,t] >= 0;
	
s.t. Calculate_Total_Dispatch {g in G_D, t in T}:
	P_thermal_dist[g,t] = sum{b in B} Block_Dispatch_d[g,b,t];

s.t. MAXIMUM_DG_ACTIVE_POWER {t in T, g in G_D}:
    P_thermal_dist[g,t] <=   Pmax_d[g];											

s.t. MINIMUM_DG_ACTIVE_POWER {t in T, g in G_D}:
    P_thermal_dist[g,t] >=   Pmin_d[g];		
    
s.t. MAXIMUM_DG_REACTIVE_POWER {t in T, g in G_D}:
    Q_thermal_dist[g,t] <=   Pmax_d[g];											

s.t. MINIMUM_DG_REACTIVE_POWER {t in T, g in G_D}:
    Q_thermal_dist[g,t] >=   Pmin_d[g];
    

s.t. MAXIMUM_DG_POWER {t in T, g in G_D}:
    P_thermal_dist[g,t]^2 + Q_thermal_dist[g,t]^2 <= Pmax_d[g];	
    
    
    
s.t. MAXIMUM_DG_ACTIVE_POWER_RENEWABLE {t in T, g in RES_D}:
    P_renewable_dist[g,t] <=   if RES_type_d[g]=="SOLAR" then Dist_Solar_Inj[g,t] else Dist_Wind_Inj[g,t];											

s.t. MINIMUM_DG_ACTIVE_POWER_RENEWABLE {t in T, g in RES_D}:
    P_renewable_dist[g,t] >=   0;		
    
s.t. MAXIMUM_DG_REACTIVE_POWER_RENEWABLE {t in T, g in RES_D}:
    Q_renewable_dist[g,t] <=   if RES_type_d[g]=="SOLAR" then Dist_Solar_Inj[g,t] else Dist_Wind_Inj[g,t];											

s.t. MINIMUM_DG_REACTIVE_POWER_RENEWABLE {t in T, g in RES_D}:
    Q_renewable_dist[g,t] >=   if RES_type_d[g]=="SOLAR" then -Dist_Solar_Inj[g,t] else -Dist_Wind_Inj[g,t];	
								

#---- LOAD SHIFTING LIMITS
s.t. MAXIMUM_LS_LIMIT {t in T, n in LS}:
	Dist_Shift[n,t] <= Dist_Shift_Max[n] *Dist_Load[LS_Node[n],t] * load1[scenario,t]; 			
	
s.t. MINIMUM_LS_LIMIT {t in T, n in LS}:
	-Dist_Shift[n,t] <= -Dist_Shift_Min[n] *Dist_Load[LS_Node[n],t] * load1[scenario,t]; 		

s.t. ZERO_SUM_LS {n in LS}:
	sum{t in T} Dist_Shift[n,t] = 0; 										
	
	
#---- ENERGY STORAGE SYSTEMS LIMITS
s.t. Dist_ESS_Limits  {n in ESS, t in T}: Dist_P_ESS[n,t] >= Dist_P_ESS_min[n]; 					
s.t. Dist_ESS_Limits2 {n in ESS, t in T}: Dist_P_ESS[n,t] <= Dist_P_ESS_max[n]; 					

s.t. Dispatch_Dist_SOCLimits  {n in ESS, t in T}: Dist_SOC[n,t] >= Dist_SOC_min[n];  			
s.t. Dispatch_Dist_SOCLimits2 {n in ESS, t in T}: Dist_SOC[n,t] <= Dist_SOC_max[n]; 				

s.t. Dispatch_Dist_CalculateSOC {n in ESS, t in T}: (if t = first(T) then Dist_SOC_initial[n] 
else Dist_SOC[n,prev(t)]) - Dist_P_ESS[n,t] = Dist_SOC[n,t]; 										

s.t. Dispatch_Dist_SameSOC {n in ESS}: Dist_SOC_initial[n] = Dist_SOC[n, last(T)]; 				




#---- CARBON RELATED CONSTRAINTS
s.t. Calculate_Footprint_Dist {g in G_D, t in T}:
	Footprint_dist[g,t] = Carbon_Cost_d[g] * P_thermal_dist[g,t] ;
	
s.t. Dist_Carbon_Limit {t in T}:
	 sum{g in G_D}(Footprint_dist[g,t] + Carbon_D[g,t]) <= Carbon_Limit_Dist[t] ;	
	
#s.t. Dist_Carbon_Trade_Limit1 {g in G_D, t in T}:
#	Footprint_dist[g,t] + Carbon_D[g,t] >=0;

s.t. Dist_Carbon_SE {t in T}:
	Carbon_SE[t] = -sum{g in G_D}Carbon_D[g,t];








##									TRANSMISSION PRIMAL CONSTRAINTS
# MIN(GenCosts - P_DSO*Bid)
s.t. Calculate_Trans_PowerFlow {l in Trans_Lines, t in T}: 
	Trans_Flow[l,t] - Trans_Status[l] * 1/Trans_Reactance[l] * 
	sum{n in Trans_Nodes}(Trans_Incidencia[n,l]*Trans_Theta[n,t]) = 0;					#mu


s.t. Trans_Power_Balance{n in Trans_Nodes, t in T}: 
	sum{g in G_T : G_Node[g] == n}(P_thermal_trans[g,t])
	-sum{l in Trans_Lines}(Trans_Incidencia[n,l]*Trans_Flow[l,t])
	-(if 5 == n then P_DSO[t]) 
	
	+sum{g in RES_T:RES_Node_t[g] == n and RES_type_t[g]=="SOLAR"} (Trans_Solar_Inj[g,t] )
	+sum{g in RES_T:RES_Node_t[g] == n and RES_type_t[g]=="WIND"} (Trans_Wind_Inj[g,t])
	
	= Trans_Load[n,t] * load2[scenario,t];																	#lambda


s.t. Trans_Flow_Upper{l in Trans_Lines, t in T}:  Trans_Flow[l,t] <= Trans_Capacity[l]; #omega_U
s.t. Trans_Flow_Lower{l in Trans_Lines, t in T}: -Trans_Flow[l,t] <= Trans_Capacity[l]; #omega_L


s.t. Calculate_Block_Upper_Trans {g in G_T, b in B, t in T}: Block_Dispatch_t[g,b,t] <= P_block_t[g,b];				#rho_U
s.t. Calculate_Block_Lower_Trans {g in G_T, b in B, t in T}: -Block_Dispatch_t[g,b,t] <= 0;							#rho_L
	
s.t. Calculate_Total_Dispatch_Trans {g in G_T, t in T}: P_thermal_trans[g,t] = sum{b in B} Block_Dispatch_t[g,b,t];	#varphi

s.t. Trans_Gen_Upper{g in G_T, t in T}:  P_thermal_trans[g,t] <=  Pmax_t[g]; #alpha_U
s.t. Trans_Gen_Lower{g in G_T, t in T}: -P_thermal_trans[g,t] <= -Pmin_t[g]; #alpha_L

#s.t. Trans_Gen_Upper_Ramp{g in G_T, t in T: t<>first(T)}:  P_thermal_trans[g,t-1] - P_thermal_trans[g,t] <=  0.25 * Pmax_t[g]; #alpha_U
#s.t. Trans_Gen_Lower_Ramp{g in G_T, t in T: t<>first(T)}: -P_thermal_trans[g,t-1] + P_thermal_trans[g,t] <=  0.25 * Pmax_t[g]; #alpha_L

s.t. MAXIMUM_SE_LIMIT {t in T}:    P_DSO[t] <=   SE_Capacity; #kappa_U											
s.t. MINIMUM_SE_LIMIT {t in T}:   -P_DSO[t] <=  SE_Capacity;  #kappa_L


subject to Calculate_Footprint_Trans {g in G_T, t in T}:
	Footprint_trans[g,t] = Carbon_Cost_t[g] * P_thermal_trans[g,t];		#phi
	 
s.t. Carbon_Balance_Trans {t in T}:
	sum{g in G_T} (Carbon_T[g,t]) = Carbon_SE[t];					#psi

s.t. Carbon_Limits_Trans {t in T}:									
	sum{g in G_T} (Footprint_trans[g,t] + Carbon_T[g,t]) <= Carbon_Limit_Trans[t] ;		#Xi_U





##									TRANSMISSION DUAL CONSTRAINTS
s.t. Deriv_Potencia_with_Ramp {g in G_T, t in T}:
	#(2*a_t[g] * P_thermal_trans[g,t] + b_t[g]) 
	- lambda[G_Node[g],t] 
	-phi[g,t]*Carbon_Cost_t[g] 
	- alpha_L[g,t] + alpha_U[g,t] 
	+varphi[g,t]
	= 0;
	
s.t. Deriv_Blocos_Potencia {g in G_T, b in B, t in T}:
	Price_block_t[g,b] -varphi[g,t] + rho_U[g,b,t] - rho_L[g,b,t] = 0;
	
s.t. Deriv_Teta{n in Trans_Nodes, t in T}: 
	-sum{l in Trans_Lines}(-Trans_Status[l] * (1/Trans_Reactance[l]) * Trans_Incidencia[n,l]*mu[l,t]) = 0;

s.t. Deriv_Fluxo{l in Trans_Lines, t in T}: 
	mu[l,t] - omega_L[l,t] + omega_U[l,t] 
	-sum{n in Trans_Nodes}(Trans_Incidencia[n,l] * lambda[n,t]) = 0;

s.t. Deriv_PDSO{t in T}:  
	-Bid[t] + lambda[5,t] + kappa_U[t] - kappa_L[t] = 0;

s.t. Deriv_Footprint{g in G_T, t in T}:
	+phi[g,t] + Xi_U[t] 
	= 0;

s.t. Deriv_Carbon_T {g in G_T, t in T}:
	+ psi[t] + Xi_U[t]
	 = 0;
	
s.t. Deriv_Carbon_SE {t in T}:
	 -Carbon_Price[t] -psi[t] = 0;




s.t. Omega_Lower{l in Trans_Lines, t in T}: omega_L[l,t] >= 0;
s.t. Omega_Upper{l in Trans_Lines, t in T}: omega_U[l,t] >= 0;

s.t. Alpha_Lower{g in G_T, t in T}: alpha_L[g,t] >= 0;
s.t. Alpha_Upper{g in G_T, t in T}: alpha_U[g,t] >= 0;

s.t. Rho_Lower{g in G_T, b in B, t in T}: rho_L[g,b,t] >= 0;
s.t. Rho_Upper{g in G_T, b in B, t in T}: rho_U[g,b,t] >= 0;

s.t. Kappa_Lower{t in T}: kappa_L[t] >= 0;
s.t. Kappa_Upper{t in T}: kappa_U[t] >= 0;

s.t. Xi_Upper_Limit{t in T}: Xi_U[t] >=0;


##					LINEAR KKT COMPLEMENTARY SLACKNESS

#---- Lin alpha_U
subj to slack1 {g in G_T, t in T}: alpha_U[g,t] <= M * z_alpha_U[g,t];
subj to slack2 {g in G_T, t in T}: (-P_thermal_trans[g,t] +Pmax_t[g]) <= M * (1 - z_alpha_U[g,t]);

#---- Lin alpha_L
subj to slack3 {g in G_T, t in T}: alpha_L[g,t] <= M * z_alpha_L[g,t];
subj to slack4 {g in G_T, t in T}: (P_thermal_trans[g,t] -Pmin_t[g]) <= M * (1 - z_alpha_L[g,t]);

#---- Lin kappa_U
subj to slack5 {t in T}: kappa_U[t] <= M * z_kappa_U[t];
subj to slack6 {t in T}: (-P_DSO[t] +SE_Capacity) <= M * (1 - z_kappa_U[t]);

#---- Lin kappa_L
subj to slack7 {t in T}: kappa_L[t] <= M * z_kappa_L[t];
subj to slack8 {t in T}: (P_DSO[t] +SE_Capacity) <= M * (1 - z_kappa_L[t]);

#---- Lin omega_U
subj to slack9 {l in Trans_Lines, t in T}: omega_U[l,t] <= M * z_omega_U[l,t];
subj to slack10 {l in Trans_Lines, t in T}: (-Trans_Flow[l,t] +Trans_Capacity[l]) <= M * (1 - z_omega_U[l,t]);

#---- Lin omega_L
subj to slack11 {l in Trans_Lines, t in T}: omega_L[l,t] <= M * z_omega_L[l,t];
subj to slack12 {l in Trans_Lines, t in T}: (Trans_Flow[l,t] +Trans_Capacity[l]) <= M * (1 - z_omega_L[l,t]);

#---- Lin Xi_U
subj to slack13 {t in T}: Xi_U[t] <= M * z_Xi_U[t];
subj to slack14 {t in T}: (Carbon_Limit_Trans[t] - sum{g in G_T} (Footprint_trans[g,t] + Carbon_T[g,t])) <= M * (1 - z_Xi_U[t]);



#---- Lin alpha_U
subj to slack21 {g in G_T, b in B, t in T}: rho_U[g,b,t] <= M * 100*z_rho_U[g,b,t];
subj to slack22 {g in G_T, b in B, t in T}: (-Block_Dispatch_t[g,b,t] +P_block_t[g,b]) <= 100* M * (1 - z_rho_U[g,b,t]);

#---- Lin alpha_L
subj to slack23 {g in G_T, b in B, t in T}: rho_L[g,b,t] <= M * 100*z_rho_L[g,b,t];
subj to slack24 {g in G_T, b in B, t in T}: (Block_Dispatch_t[g,b,t]) <= 100*M * (1 - z_rho_L[g,b,t]);