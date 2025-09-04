# Restrições primais e duais associadas ao problema do TSO (Nível Inferior)
Calculate_Trans_PowerFlow: Dual variable is mu[l,t]          ==
Trans_Power_Balance: Dual variable is lambda[n,t]            ==
Trans_Flow_Upper: Dual variable is omega_U[l,t]              ineq
Trans_Flow_Lower: Dual variable is omega_L[l,t]              ineq
Calculate_Block_Upper_Trans: Dual variable is rho_U[g,b,t]   ineq
Calculate_Block_Lower_Trans: Dual variable is rho_L[g,b,t]   ineq
Calculate_Total_Dispatch_Trans: Dual variable is varphi[g,t] ==
Trans_Gen_Upper: Dual variable is alpha_U[g,t]               ineq
Trans_Gen_Lower: Dual variable is alpha_L[g,t]               ineq
MAXIMUM_SE_LIMIT: Dual variable is kappa_U[t]                ineq
MINIMUM_SE_LIMIT: Dual variable is kappa_L[t]                ineq
Calculate_Footprint_Trans: Dual variable is phi[g,t]         ==
Carbon_Balance_Trans: Dual variable is psi[t]                ==
Carbon_Limits_Trans: Dual variable is Xi_U[t]                ineq

#Derivadas da lagrangeana em relação às primais 
Deriv_Potencia_with_Ramp
Deriv_Blocos_Potencia
Deriv_Teta
Deriv_Fluxo
Deriv_PDSO
Deriv_Footprint
Deriv_Carbon_T
Deriv_Carbon_SE


# Lagrangian of the TSO's problem (Lower-Level):

# L = sum{t in T, g in G_T, b in B} (Price_block_t[g,b] * Block_Dispatch_t[g,b,t])
#
#     + sum{l in Trans_Lines, t in T} mu[l,t] * (Trans_Flow[l,t] - Trans_Status[l] * 1/Trans_Reactance[l] * sum{n in Trans_Nodes}(Trans_Incidencia[n,l]*Trans_Theta[n,t]))
#
#     + sum{n in Trans_Nodes, t in T} lambda[n,t] * (
#           Trans_Load[n,t] * load2[scenario,t]
#           - (sum{g in G_T : G_Node[g] == n}(P_thermal_trans[g,t])
#           - sum{l in Trans_Lines}(Trans_Incidencia[n,l]*Trans_Flow[l,t])
#           - (if 5 == n then P_DSO[t])
#           + sum{g in RES_T:RES_Node_t[g] == n and RES_type_t[g]=="SOLAR"} (Trans_Solar_Inj[g,t])
#           + sum{g in RES_T:RES_Node_t[g] == n and RES_type_t[g]=="WIND"} (Trans_Wind_Inj[g,t]))
#       )
#
#     + sum{l in Trans_Lines, t in T} omega_U[l,t] * (Trans_Flow[l,t] - Trans_Capacity[l])
#     + sum{l in Trans_Lines, t in T} omega_L[l,t] * (-Trans_Flow[l,t] - Trans_Capacity[l])
#
#     + sum{g in G_T, b in B, t in T} rho_U[g,b,t] * (Block_Dispatch_t[g,b,t] - P_block_t[g,b])
#     + sum{g in G_T, b in B, t in T} rho_L[g,b,t] * (-Block_Dispatch_t[g,b,t])
#
#     + sum{g in G_T, t in T} varphi[g,t] * (P_thermal_trans[g,t] - sum{b in B} Block_Dispatch_t[g,b,t])
#
#     + sum{g in G_T, t in T} alpha_U[g,t] * (P_thermal_trans[g,t] - Pmax_t[g])
#     + sum{g in G_T, t in T} alpha_L[g,t] * (-P_thermal_trans[g,t] + Pmin_t[g])
#
#     + sum{t in T} kappa_U[t] * (P_DSO[t] - SE_Capacity)
#     + sum{t in T} kappa_L[t] * (-P_DSO[t] - SE_Capacity)
#
#     + sum{g in G_T, t in T} phi[g,t] * (Footprint_trans[g,t] - Carbon_Cost_t[g] * P_thermal_trans[g,t])
#
#     + sum{t in T} psi[t] * (sum{g in G_T} (Carbon_T[g,t]) - Carbon_SE[t])
#
#     + sum{t in T} Xi_U[t] * (sum{g in G_T} (Footprint_trans[g,t] + Carbon_T[g,t]) - Carbon_Limit_Trans[t])
#