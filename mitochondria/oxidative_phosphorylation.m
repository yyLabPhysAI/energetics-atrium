function [...
    V_He, dC_FLV, V_He_F, dC_NADH, V_Hu, V_H_Leak, dC_ADP_m, V_ANT, V_O2]...
    = oxidative_phosphorylation(V_SL, V_IDH, V_KGDH, V_MDH, V_SDH, ...
    delta_Psi_m, C_NADH, C_NAD, C_ATP_m, C_ADP_m, Ca_m, ATP_i, ...
    ADP_i, dC_AcCoA,  data)

r_a = data.r_a;
r_c1 = data.r_c1;
r_c2 = data.r_c2;
r_1 = data.r_1;
r_2 = data.r_2;
r_3 = data.r_3;
rho_res = data.rho_res;
K_res = data.K_res;
rho_resF = data.rho_resF;
psi_B = data.psi_B;
K_resF = data.K_resF;
r_b = data.r_b;
FADH2 = data.FADH2;
FADH = data.FADH;
p_a = data.p_a;
p_b = data.p_b;
p_c1 = data.p_c1;
p_c2 = data.p_c2;
p_1 = data.p_1;
p_2 = data.p_2;
p_3 = data.p_3;
KCaATP = data.KCaATP;
rho_F1 = data.rho_F1;
K_F1 = data.K_F1;
V_ant_max = data.V_ant_max;
h_ANT = data.h_ANT;
g_H = data.g_H;
delta_pH = data.delta_pH;
R = data.R;
T = data.T;
F = data.F;

% The respiration-driven proton pump
g = 1/2; % voltage correction factor
delta_mu_h = -2.303.*(R.*T./F).*delta_pH + delta_Psi_m;
A_res = (R.*T./F).*log(K_res.*sqrt(C_NADH./C_NAD));
V_O2 = 9.1979*0.2515*7.0810*0.5.*rho_res.*((r_a + r_c1.*exp(6.*F.*psi_B./(R.*T))).*exp(A_res.*F./(R.*T)) - r_a.*exp((g.*6.*F.*delta_mu_h)./(R.*T)) + r_c2.*exp(A_res.*F./(R.*T)).*exp(g.*6.*F.*delta_mu_h./(R.*T)))./...
                   ((1 + r_1.*exp(A_res.*F./(R.*T))) .* exp(6.*F.*psi_B./(R.*T)) + (r_2 + r_3.*exp(A_res.*F./(R.*T))).*exp(g.*6.*F.*delta_mu_h./(R.*T)));
V_He = 6.*rho_res.*(r_a.*exp(A_res.*F./(R.*T)) - (r_a + r_b).*exp(g.*6.*F.*delta_mu_h./(R.*T)))./...
                 ((1 + r_1.*exp(A_res.*F./(R.*T))) .* exp(6.*F.*psi_B./(R.*T)) + (r_2 + r_3.*exp(A_res.*F./(R.*T))).*exp(g.*6.*F.*delta_mu_h./(R.*T)));
dC_FLV = V_SDH - V_O2;
A_res_F = (R.*T./F).*log(K_resF.*sqrt(FADH2./FADH)); % 
V_He_F = 4.*rho_resF.*(r_a.*exp(A_res_F.*F./(R.*T)) - (r_a + r_b).*exp(g.*6.*F.*delta_mu_h./(R.*T)))./...
                    ((1 + r_1.*exp(A_res_F.*F./(R.*T))).*exp(6.*F.*psi_B./(R.*T)) + (r_2 + r_3.*exp(A_res_F.*F./(R.*T))).*exp(g.*6.*F.*delta_mu_h./(R.*T)));
dC_NADH = -V_O2 + V_IDH + V_KGDH + V_MDH + dC_AcCoA;

% The F_1F_0-ATPase
A_F1 = (R.*T./F).*log(K_F1.*C_ATP_m./(C_ADP_m.*data.Pi));
V_ATPase = -rho_F1.*((10.^2.*p_a + p_c1.*exp(3.*F.*psi_B./(R.*T))).*exp(F.*A_F1./(R.*T)) - (p_a.*exp(3.*F.*delta_mu_h./(R.*T)) + p_c2.*exp(F.*A_F1./(R.*T)).*exp(3.*F.*delta_mu_h./(R.*T))))./...
                   ((1 + p_1.*exp(A_F1.*F./(R.*T))).*exp(3.*F.*psi_B./(R.*T)) + (p_2 + p_3.*exp(A_F1.*F./(R.*T))).*exp(3.*F.*delta_mu_h./(R.*T))).*(1 - exp(-Ca_m./KCaATP));
V_Hu = -3.*rho_F1.*(10.^2.*p_a.*(1 + exp(A_F1./(R.*T))) - (p_a + p_b).*exp(3.*F.*delta_mu_h./(R.*T)))./...
                  ((1 + p_1.*exp(A_F1.*F./(R.*T))).*exp(3.*F.*psi_B./(R.*T)) + (p_2 + p_3.*exp(A_F1.*F./(R.*T))).*exp(3.*F.*delta_mu_h./(R.*T)));

% Adenine nucleotide translocator (ANT) and proton leak
V_H_Leak = g_H.*delta_mu_h;
V_ANT = V_ant_max.*(0.75.*(1 - 1e-4*((0.25.*ATP_i.*0.45.*C_ADP_m)./(0.17.*ADP_i.*0.025.*C_ATP_m))).*exp(-(F.*delta_Psi_m)./(R.*T)))./...
                  ((1 + ((0.25.*ATP_i)./(0.225.*ADP_i)).*exp((-h_ANT.*F.*delta_Psi_m)./(R.*T))).*(1 + ((0.45.*C_ADP_m)./(0.025.*C_ATP_m ))));
dC_ADP_m = V_ANT - V_ATPase - V_SL;

end