function [...
    dC_ISOC, dC_aKG, dC_SCoA, dC_Suc, dC_FUM, dC_MAL, dC_OAA, V_SL, V_IDH, V_KGDH, V_MDH, V_SDH...
    ] = TCA_cycle(C_ISOC,C_aKG,C_SCoA,C_Suc,C_FUM,C_MAL,C_OAA, C_NAD, C_ADP_m, Ca_m, C_NADH, C_ATP_m, data)


C_AcCoA=data.C_AcCoA;k_cat_cs=data.k_cat_cs;
E_T_cs=data.E_T_cs;K_M_AcCoA=data.K_M_AcCoA;K_M_OAA=data.K_M_OAA;C_K_int=data.C_K_int;
k_f_ACO=data.k_f_ACO;k_E_ACO=data.k_E_ACO;K_ADP_a=data.K_ADP_a;K_Ca_a=data.K_Ca_a;
K_i_NADH=data.K_i_NADH;k_cat_IDH=data.k_cat_IDH;E_T_IDH=data.E_T_IDH;C_H=data.C_H;k_h_1=data.k_h_1;
k_h_2=data.k_h_2;K_M_ISOC=data.K_M_ISOC;K_M_NAD=data.K_M_NAD;
E_T_KGDH=data.E_T_KGDH;k_cat_KGDH=data.k_cat_KGDH;K_M_aKG=data.K_M_aKG;
K_M_NAD_new=data.K_M_NAD_new;n_aKG=data.n_aKG;C_Mg=data.C_Mg;k_f_SL=data.k_f_SL;
k_E_SL=data.k_E_SL;C_CoA=data.C_CoA;k_cat_SDH=data.k_cat_SDH;E_T_SDH=data.E_T_SDH;
K_M_SUC=data.K_M_SUC;K_i_FUM=data.K_i_FUM;K_i_sdh_OAA=data.K_i_sdh_OAA;k_f_FH=data.k_f_FH;
K_E_FH=data.K_E_FH;k_h1=data.k_h1;k_h2=data.k_h2;k_h3=data.k_h3;k_h4=data.k_h4;k_offset=data.k_offset;
k_cat_MDH=data.k_cat_MDH;E_T_MDH=data.E_T_MDH;K_M_MAL=data.K_M_MAL;K_i_OAA=data.K_i_OAA;
K_M_NAD_mdh=data.K_M_NAD_mdh;C_GLU=data.C_GLU;k_f_AAT=data.k_f_AAT;K_E_AAT=data.K_E_AAT;
k_ASP=data.k_ASP;K_D_Ca=data.K_D_Ca;K_D_Mg=data.K_D_Mg;

% Citrate
C_CIT = C_K_int - (C_ISOC + C_aKG + C_SCoA + C_Suc + C_FUM + C_MAL + C_OAA);

% Citrate synthase (CS)
V_CS = k_cat_cs*E_T_cs*(1 + (K_M_AcCoA/C_AcCoA) + (K_M_OAA/C_OAA) + (K_M_AcCoA/C_AcCoA) + (K_M_OAA/C_OAA))^(-1);

% Aspartate amino transferase (AAT)
V_AAT = k_f_AAT*C_OAA*C_GLU*(k_ASP*K_E_AAT)/(k_ASP*K_E_AAT + C_aKG*k_f_AAT);

% Aconitase (ACO)
V_ACO = k_f_ACO*(C_CIT - (C_ISOC/k_E_ACO));

% Isocitrate dehydrogenase (IDH)
Ni = 2;

f_a_IDH = ((1 + (C_ADP_m/K_ADP_a))*(1 + (Ca_m/K_Ca_a)))^(-1);
f_i_IDH = 1 + (C_NADH/K_i_NADH);
V_IDH = k_cat_IDH*E_T_IDH*(1 + (C_H/k_h_1) + (k_h_2/C_H) + f_i_IDH*(K_M_NAD/C_NAD) + f_a_IDH*(K_M_ISOC/C_ISOC)^Ni + f_a_IDH*f_i_IDH*(K_M_NAD/C_NAD)*(K_M_ISOC/C_ISOC)^Ni)^-1;

dC_ISOC = V_ACO - V_IDH;


% Alpha-ketoglutarate dehydrogenase (KGDH)
f_a_KGDH = ((1 + (C_Mg/K_D_Mg))*(1 + (Ca_m/K_D_Ca)))^(-1);
V_KGDH = (k_cat_KGDH*E_T_KGDH)/(1 + f_a_KGDH*(K_M_aKG/C_aKG)^(n_aKG)+f_a_KGDH*((K_M_NAD_new)/C_NAD)); % ?????? why not K_M_NAD??

dC_aKG = V_IDH - V_KGDH + V_AAT;

% Succinyl CoA lyase (SL)
V_SL = k_f_SL*((C_SCoA*C_ADP_m - (C_Suc*C_ATP_m*C_CoA)/(k_E_SL)));
dC_SCoA = V_KGDH - V_SL;

% Succinate dehydrogenase (SDH)
V_SDH = (k_cat_SDH*E_T_SDH)/...
        (1 + (K_M_SUC/C_Suc)*(1 + C_OAA/K_i_sdh_OAA)*(1 + C_FUM/K_i_FUM));

dC_Suc = V_SL - V_SDH;

% Fumarate hydratase (FH)
V_FH = k_f_FH*(C_FUM - C_MAL/(K_E_FH));
dC_FUM = V_SDH - V_FH;

% Malate dehydrogenase (MDH)
f_h_a = (1 + C_H/k_h1 + C_H^2/(k_h1*k_h2))^-1 + k_offset;
f_h_i = (1 + k_h3/C_H + (k_h3*k_h4)/C_H^2)^-2;
V_MDH = (k_cat_MDH*E_T_MDH*f_h_a*f_h_i)/(1 + (K_M_MAL/C_MAL)*(1 + (C_OAA/K_i_OAA)) + (K_M_NAD_mdh/C_NAD) + ((K_M_MAL/C_MAL)*(1 + (C_OAA/K_i_OAA))*(K_M_NAD_mdh))/(C_NAD)); % ???????? Why different K_M_NAD_mdh????

dC_MAL = V_FH - V_MDH;

% Oxaloacetate
dC_OAA = V_MDH - V_CS - V_AAT;


end