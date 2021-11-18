function [dO_c, dO_TnCa, dO_TnMgCa, dO_TnMgMg, ...
          dO_Calse, phi_ca_i, dCa_up, dCa_rel, ...
          dF_1, dF_2, dF_3, dfca, I_up, I_rel] = ...
          SR_calcium_handling(ATP_i, Ca_i, Ca_up, ...
          Ca_rel, F_1, F_2, F_3, O_c, O_TnCa, O_TnMgCa, ...
          O_TnMgMg, O_Calse, V, data, fca)

F = data.F;
V_up = data.V_up;
V_rel = data.V_rel;
I_upMax = 2800/1000;
tau_tr = 0.01;
alpha_rel = 200000;
K_cyCa = data.K_cyCa;
K_srCa = data.K_srCa;
K_xcs = data.K_xcs;

% Uptake current of Ca++ to the SR
I_up = I_upMax*(ATP_i/7.977)*((Ca_i/K_cyCa - K_xcs.^2*Ca_up/K_srCa)/...
                             (((Ca_i + K_cyCa)/K_cyCa) + (K_xcs*(Ca_up + K_srCa)/K_srCa)));

% Transfer current between the uptake and release compartments of SR
I_tr = (Ca_up - Ca_rel)*2*F*V_up/tau_tr;

% Release current of Ca++ from the SR
I_rel = alpha_rel*(F_2/(F_2 + 0.25))^2*(Ca_rel - Ca_i);

% calculation of the Ca++ concentration in the different compartments
Mg_i = 2.5; % [mM]

dO_c     = 200000*Ca_i*(1 - O_c) - 476*O_c;
dO_TnCa  = 78400*Ca_i*(1 - O_TnCa) - 392*O_TnCa;
dO_TnMgCa  = 200000*Ca_i*(1 - O_TnMgCa - O_TnMgMg) - 6.6*O_TnMgCa;
dO_TnMgMg  = 2000*Mg_i*(1 - O_TnMgCa - O_TnMgMg) - 666*O_TnMgMg;
dO_Calse = 480*Ca_rel*(1 - O_Calse) - 400*O_Calse;
phi_ca_i = 0.08*dO_TnCa + 0.16*dO_TnMgCa + 0.045*dO_c;

dCa_up = (I_up - I_tr)/(2*F*V_up);
dCa_rel = (I_tr - I_rel)/(2*F*V_rel) - 31*dO_Calse;

% opening and closing of SR Ca++ release channels
K_Mrel = 0.0003;
k_recov = 0.815;
k_act  = 240*exp((V - 20)/12.5) + 203.8*(Ca_i/(Ca_i + K_Mrel))^4;
k_inact = 33.96 + 339.6*(Ca_i/(Ca_i + K_Mrel))^4;

dF_1 = k_recov*F_3 - k_act*F_1;
dF_2 = k_act*F_1 - k_inact*F_2;
dF_3 = k_inact*F_2 - k_recov*F_3;

% variable not really used in the paper - implemented for future use
dfca = 0*fca;


end