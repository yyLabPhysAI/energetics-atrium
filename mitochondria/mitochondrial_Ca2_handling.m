function [dCa_m, J_uni, J_NaCa] = mitochondrial_Ca2_handling(delta_Psi_m, Ca_m, Ca_i, data)

P_Ca=data.P_Ca;
Z_Ca=data.Z_Ca;
alpha_m=data.alpha_m;
alpha_e=data.alpha_e;
V_NC=data.V_NC;
Na_e=data.Na_e;
Na_m=data.Na_m;
beta_Ca=data.beta_Ca*10;
K_Na=data.K_Na;
K_Ca = data.K_Ca;
R=data.R;
T=data.T;
F=data.F;

% The Ca2+ uniporter
J_uni = P_Ca.*(Z_Ca.*delta_Psi_m.*F./(R.*T)).*(alpha_m.*Ca_m.*exp(-Z_Ca.*delta_Psi_m.*F./(R.*T)) - alpha_e.*Ca_i)./...
    (exp(-Z_Ca.*delta_Psi_m.*F./(R.*T)) - 1);

% The Na+./Ca2+ exchanger
J_NaCa = V_NC.*(exp(0.5.*delta_Psi_m.*F./(R.*T)).*Na_e.^3.*Ca_m./(K_Na.^3.*K_Ca) - exp(-0.5.*delta_Psi_m.*F./(R.*T)).*(Na_m.^3.*Ca_i)./(K_Na.^3.*K_Ca))./...
    (1 + Na_e.^3./K_Na.^3 + Ca_m./K_Ca + Na_e.^3.*Ca_m./(K_Na.^3.*K_Ca) + Na_m.^3./K_Na.^3 + Ca_i./K_Ca + Na_m.^3.*Ca_i./(K_Na.^3.*K_Ca));

dCa_m = beta_Ca.*(J_uni - J_NaCa);

end