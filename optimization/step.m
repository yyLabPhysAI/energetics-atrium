function [new_data, new_y0] = step(data, y0, eta)

new_data = data;

new_data.V_uni_max = new_data.V_uni_max*exp(eta*randn(1)); 
new_data.psi0 = new_data.psi0*exp(eta*randn(1)); 
new_data.K_act = new_data.K_act*exp(eta*randn(1)); 
new_data.K_trans = new_data.K_trans*exp(eta*randn(1)); 
new_data.L = new_data.L*exp(eta*randn(1)); 
new_data.n_a = new_data.n_a*exp(eta*randn(1)); 
new_data.V_NaCa_max = new_data.V_NaCa_max*exp(eta*randn(1)); 
new_data.b = new_data.b*exp(eta*randn(1)); 
new_data.K_Na = new_data.K_Na*exp(eta*randn(1)); 
new_data.K_Ca = new_data.K_Ca*exp(eta*randn(1)); 
new_data.n = new_data.n*exp(eta*randn(1)); 
new_data.delta = new_data.delta*exp(eta*randn(1)); 
new_data.r_a = new_data.r_a*exp(eta*randn(1)); 
new_data.r_c1 = new_data.r_c1*exp(eta*randn(1)); 
new_data.r_c2 = new_data.r_c2*exp(eta*randn(1)); 
new_data.r_1 = new_data.r_1*exp(eta*randn(1)); 
new_data.r_2 = new_data.r_2*exp(eta*randn(1)); 
new_data.r_3 = new_data.r_3*exp(eta*randn(1)); 
new_data.rho_res = new_data.rho_res*exp(eta*randn(1)); 
new_data.K_res = new_data.K_res*exp(eta*randn(1)); 
new_data.rho_resF = new_data.rho_resF*exp(eta*randn(1)); 
new_data.psi_B = new_data.psi_B*exp(eta*randn(1)); 
new_data.G = new_data.G*exp(eta*randn(1)); 
new_data.K_resF = new_data.K_resF*exp(eta*randn(1)); 
new_data.r_b = new_data.r_b*exp(eta*randn(1)); 
new_data.FADH2 = new_data.FADH2*exp(eta*randn(1)); 
new_data.FADH = new_data.FADH*exp(eta*randn(1)); 
new_data.p_a = new_data.p_a*exp(eta*randn(1)); 
new_data.p_b = new_data.p_b*exp(eta*randn(1)); 
new_data.p_c1 = new_data.p_c1*exp(eta*randn(1)); 
new_data.p_c2 = new_data.p_c2*exp(eta*randn(1)); 
new_data.p_1 = new_data.p_1*exp(eta*randn(1)); 
new_data.p_2 = new_data.p_2*exp(eta*randn(1)); 
new_data.p_3 = new_data.p_3*exp(eta*randn(1)); 
new_data.KCaATP = new_data.KCaATP*exp(eta*randn(1)); 
new_data.rho_F1 = new_data.rho_F1*exp(eta*randn(1)); 
new_data.K_F1 = new_data.K_F1*exp(eta*randn(1)); 
new_data.Pi = new_data.Pi*exp(eta*randn(1)); 
new_data.C_A = new_data.C_A*exp(eta*randn(1)); 
new_data.V_ant_max = new_data.V_ant_max*exp(eta*randn(1)); 
new_data.h_ANT = new_data.h_ANT*exp(eta*randn(1)); 
new_data.g_H = new_data.g_H*exp(eta*randn(1)); 
new_data.delta_pH = new_data.delta_pH*exp(eta*randn(1)); 
new_data.C_PN = new_data.C_PN*exp(eta*randn(1)); 
new_data.C_mito = new_data.C_mito*exp(eta*randn(1)); 
new_data.Cc = new_data.Cc*exp(eta*randn(1)); 
new_data.C_AcCoA = new_data.C_AcCoA*exp(eta*randn(1)); 
new_data.k_cat_cs = new_data.k_cat_cs*exp(eta*randn(1)); 
new_data.E_T_cs = new_data.E_T_cs*exp(eta*randn(1)); 
new_data.K_M_AcCoA = new_data.K_M_AcCoA*exp(eta*randn(1)); 
new_data.K_M_OAA = new_data.K_M_OAA*exp(eta*randn(1)); 
new_data.C_K_int = new_data.C_K_int*exp(eta*randn(1)); 
new_data.k_f_ACO = new_data.k_f_ACO*exp(eta*randn(1)); 
new_data.k_E_ACO = new_data.k_E_ACO*exp(eta*randn(1)); 
new_data.K_ADP_a = new_data.K_ADP_a*exp(eta*randn(1)); 
new_data.K_Ca_a = new_data.K_Ca_a*exp(eta*randn(1)); 
new_data.K_i_NADH = new_data.K_i_NADH*exp(eta*randn(1)); 
new_data.k_cat_IDH = new_data.k_cat_IDH*exp(eta*randn(1)); 
new_data.E_T_IDH = new_data.E_T_IDH*exp(eta*randn(1)); 
new_data.C_H = new_data.C_H*exp(eta*randn(1)); 
new_data.k_h_1 = new_data.k_h_1*exp(eta*randn(1)); 
new_data.k_h_2 = new_data.k_h_2*exp(eta*randn(1)); 
new_data.K_M_ISOC = new_data.K_M_ISOC*exp(eta*randn(1)); 
new_data.N_i = new_data.N_i*exp(eta*randn(1)); 
new_data.K_M_NAD = new_data.K_M_NAD*exp(eta*randn(1)); 
new_data.K_M_Mg = new_data.K_M_Mg*exp(eta*randn(1)); 
new_data.K_M_Ca = new_data.K_M_Ca*exp(eta*randn(1)); 
new_data.E_T_KGDH = new_data.E_T_KGDH*exp(eta*randn(1)); 
new_data.k_cat_KGDH = new_data.k_cat_KGDH*exp(eta*randn(1)); 
new_data.K_M_aKG = new_data.K_M_aKG*exp(eta*randn(1)); 
new_data.K_M_NAD_new = new_data.K_M_NAD_new*exp(eta*randn(1)); 
new_data.n_aKG = new_data.n_aKG*exp(eta*randn(1)); 
new_data.C_Mg = new_data.C_Mg*exp(eta*randn(1)); 
new_data.k_f_SL = new_data.k_f_SL*exp(eta*randn(1)); 
new_data.k_E_SL = new_data.k_E_SL*exp(eta*randn(1)); 
new_data.C_CoA = new_data.C_CoA*exp(eta*randn(1)); 
new_data.k_cat_SDH = new_data.k_cat_SDH*exp(eta*randn(1)); 
new_data.E_T_SDH = new_data.E_T_SDH*exp(eta*randn(1)); 
new_data.K_M_SUC = new_data.K_M_SUC*exp(eta*randn(1)); 
new_data.K_i_FUM = new_data.K_i_FUM*exp(eta*randn(1)); 
new_data.K_i_sdh_OAA = new_data.K_i_sdh_OAA*exp(eta*randn(1)); 
new_data.k_f_FH = new_data.k_f_FH*exp(eta*randn(1)); 
new_data.K_E_FH = new_data.K_E_FH*exp(eta*randn(1)); 
new_data.k_h1 = new_data.k_h1*exp(eta*randn(1)); 
new_data.k_h2 = new_data.k_h2*exp(eta*randn(1)); 
new_data.k_h3 = new_data.k_h3*exp(eta*randn(1)); 
new_data.k_h4 = new_data.k_h4*exp(eta*randn(1)); 
new_data.k_offset = new_data.k_offset*exp(eta*randn(1)); 
new_data.k_cat_MDH = new_data.k_cat_MDH*exp(eta*randn(1)); 
new_data.E_T_MDH = new_data.E_T_MDH*exp(eta*randn(1)); 
new_data.K_M_MAL = new_data.K_M_MAL*exp(eta*randn(1)); 
new_data.K_i_OAA = new_data.K_i_OAA*exp(eta*randn(1)); 
new_data.K_M_NAD_mdh = new_data.K_M_NAD_mdh*exp(eta*randn(1)); 
new_data.C_GLU = new_data.C_GLU*exp(eta*randn(1)); 
new_data.k_f_AAT = new_data.k_f_AAT*exp(eta*randn(1)); 
new_data.K_E_AAT = new_data.K_E_AAT*exp(eta*randn(1)); 
new_data.k_ASP = new_data.k_ASP*exp(eta*randn(1)); 
new_data.C_T = new_data.C_T*exp(eta*randn(1)); 
new_data.k_CK_cyto = new_data.k_CK_cyto*exp(eta*randn(1)); 
new_data.k_CK_mito = new_data.k_CK_mito*exp(eta*randn(1)); 
new_data.k_tr_Cr = new_data.k_tr_Cr*exp(eta*randn(1)); 
new_data.K_EQ = new_data.K_EQ*exp(eta*randn(1)); 
new_data.V_ATPase_cyto = new_data.V_ATPase_cyto*exp(eta*randn(1)); 
new_data.V_myo = new_data.V_myo*exp(eta*randn(1)); 
new_data.V_mito = new_data.V_mito*exp(eta*randn(1)); 
new_data.A_cap = new_data.A_cap*exp(eta*randn(1)); 
new_data.P_Ca = new_data.P_Ca*exp(eta*randn(1)); 
new_data.Z_Ca = new_data.Z_Ca*exp(eta*randn(1)); 
new_data.alpha_m = new_data.alpha_m*exp(eta*randn(1)); 
new_data.alpha_e = new_data.alpha_e*exp(eta*randn(1)); 
new_data.V_NC = new_data.V_NC*exp(eta*randn(1)); 
new_data.Na_e = new_data.Na_e*exp(eta*randn(1)); 
new_data.Na_m = new_data.Na_m*exp(eta*randn(1)); 
new_data.beta_Ca = new_data.beta_Ca*exp(eta*randn(1)); 
new_data.K_D_Ca = new_data.K_D_Ca*exp(eta*randn(1)); 
new_data.K_D_Mg = new_data.K_D_Mg*exp(eta*randn(1)); 


new_y0 = y0;

new_y0(36) = new_y0(36)*exp(eta*randn(1));
new_y0(37) = new_y0(37)*exp(eta*randn(1));
new_y0(38) = new_y0(38)*exp(eta*randn(1));
new_y0(39) = new_y0(39)*exp(eta*randn(1));
new_y0(40) = new_y0(40)*exp(eta*randn(1));
new_y0(41) = new_y0(41)*exp(eta*randn(1));
new_y0(42) = new_y0(42)*exp(eta*randn(1));
new_y0(43) = new_y0(43)*exp(eta*randn(1));
new_y0(44) = new_y0(44)*exp(eta*randn(1));
new_y0(45) = new_y0(45)*exp(eta*randn(1));
new_y0(46) = new_y0(46)*exp(eta*randn(1));
new_y0(47) = new_y0(47)*exp(eta*randn(1));
new_y0(48) = new_y0(48)*exp(eta*randn(1));
new_y0(49) = new_y0(49)*exp(eta*randn(1));
new_y0(50) = new_y0(50)*exp(eta*randn(1));
new_y0(51) = new_y0(51)*exp(eta*randn(1));

end 