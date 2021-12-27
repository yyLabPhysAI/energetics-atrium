function dX = model(t, X, data, f_stim, start_stim, end_stim)

% Names for the state vector values
V           = X(1);
P_a         = X(2);
P_i         = X(3);
n           = X(4);
r           = X(5);
s_1         = X(6);
s_2         = X(7);
s_3         = X(8);
m           = X(9);
h_1         = X(10);
h_2         = X(11);
d_L         = X(12);
f_L         = X(13);
d_T         = X(14);
f_T         = X(15);
Na_i        = X(16);
Ca_up       = X(17);
Ca_rel      = X(18);
Ca_i        = X(19);
O_c         = X(20);
O_TnCa      = X(21);
O_TnMgCa    = X(22);
O_TnMgMg    = X(23);
O_Calse     = X(24);
K_o         = X(25);
K_i         = X(26);
F_1         = X(27);
F_2         = X(28);
F_3         = X(29);
% fca       = X(30), not used
SL          = X(31);
A           = X(32);
TT          = X(33);
U           = X(34);
V_e         = X(35);
ATP_i       = X(36);
Ca_m        = X(37);
% C_ATP_ic  = X(38), not used
% C_CrP_i   = X(39), not used
% C_CrP_ic  = X(40), not used
C_ADP_m     = X(41);
C_NADH      = X(42);
delta_Psi_m = X(43);
C_ISOC      = X(44);
C_aKG       = X(45);
C_SCoA      = X(46);
C_Suc       = X(47);
C_FUM       = X(48);
C_MAL       = X(49);
C_OAA       = X(50);
C_FLV       = X(51);


%% Unpack constants to simplefy the code
F=data.F;

%% External stimulation
Is=data.Is;
Cm=data.Cm;
impulseFactor = 100;

if ((mod(t,1./f_stim) >= start_stim && mod(t,1./f_stim) < 0.001*impulseFactor + start_stim) && (t <= end_stim))
    I_stim=Is./impulseFactor;
else
    I_stim=0;
end

%% Action potential, membrane currents and intracellular Ca2+

% Nernst potentials
Na_o = data.Na_o; Ca_o = data.Ca_o;

E_k   = nernst(K_i, K_o, 1, data);
E_Na  = nernst(Na_i, Na_o, 1, data);
E_Ca  = nernst(Ca_i, Ca_o, 2, data);

% Fast delayed rectifier K+ current
[I_Kr, dP_a, dP_i] = fast_delayed_rectifier_k(V, P_a, P_i, E_k);

% Slow delayed rectifier K+ current
[I_Ks, dn] = slow_delayed_rectifier_k(V, n, E_k);

% Inward rectifier K+ channel
I_k1 = inward_rectifier_k(V, E_k, K_o, data);

% Transient outward K+ channel
[I_Kto, dr, ds_1, ds_2, ds_3] = transient_outward_k(r, s_1, s_2, s_3, V, E_k);

% Na+-K+ ATPase
I_NaK = NaK_ATPase(K_o, Na_i, V);

% Na+-Ca(2+) exchanger (NCX)
I_NaCa = NCX(Na_i, Ca_i, V, data);

% Fast Na+ channels
[I_Na, dm, dh_1, dh_2] = fast_na(V, m, h_1, h_2, E_Na, data);

% Membranal Ca(2+) pump
I_Cap_max = 9.509;
I_Cap = I_Cap_max.*(Ca_i./(Ca_i + 0.0002));

% L-type Ca(2+) channels
[I_CaL, df_L, dd_L] = L_type_Ca(V, d_L, f_L);

% T-Type calcium channels
[dd_T , df_T, I_CaT] = T_Type_calcium(V, d_T, f_T);

% Background currents
g_NaB = 0.03;
g_CaB = 0.03;

I_NaB = g_NaB.*(V - E_Na);
I_CaB = g_CaB.*(V - E_Ca);

% Calculate the total current
I_tot = I_Kr + I_Ks + I_k1 + I_Kto + I_NaK + I_NaCa + I_Na + I_NaB + I_CaL + I_CaT + I_Cap + I_CaB;

% Sodium and potassium concentration in cytoplasm
V_i = data.V_i;
V_Ca = data.V_Ca;
V_c = data.V_c;

dNa_i = (-3.*I_NaK - 3.*I_NaCa - I_NaB - I_Na)./(F.*V_i);
dK_o = (-2.*I_NaK + I_Kr + I_Ks + I_Kto + I_k1)./(F.*V_c);
dK_i = (2.*I_NaK - I_Kr - I_Ks - I_Kto - I_k1)./(F.*V_i);

% Change in membrane potential
dV = -(I_tot + I_stim)./...
    Cm;

%% Sarcoplasmic Reticulum (SR) and Calcium Handling
[dO_c, dO_TnCa, dO_TnMgCa, dO_TnMgMg, dO_Calse, phi_ca_i, ...
    dCa_up, dCa_rel, dF_1, dF_2, dF_3, I_up, I_rel] ...
    = ...
    SR_calcium_handling(ATP_i, Ca_i, Ca_up, Ca_rel, F_1, F_2, ...
    F_3, O_c, O_TnCa, O_TnMgCa, O_TnMgMg, O_Calse, V, data);

dCa_i = ((2.*I_NaCa - I_CaL - I_CaT - I_Cap - I_CaB - I_up + I_rel)./(2.*V_Ca.*F) - phi_ca_i);

%% Mitochondrial energy metabolism, Ca2+ dynamics and oxygen consumption

% Mitochondrial energetics and EC coupling
C_A_m = data.C_A_m;
C_A_i = data.C_A_i;
ADP_i = C_A_i - ATP_i;
C_ATP_m = C_A_m - C_ADP_m;

% Force generation and energy consumption
[dSL, dA, dTT, dU, dV_e, Force]...
    = ...
    force_generation(V_e, SL, TT, A, U, data, Ca_i, ATP_i, ADP_i);
[V_AM, ATP_XB] = force_energy_consumption(Force, ATP_i, A, data);

% Oxidation state
C_PN = data.C_PN;
C_NAD = C_PN - C_NADH;

% The TCA cycle
[dC_ISOC, dC_aKG, dC_SCoA, dC_Suc, dC_FUM, dC_MAL, dC_OAA, V_SL, V_IDH, ...
    V_KGDH, V_MDH, V_SDH] ...
    = ...
    TCA_cycle(C_ISOC,C_aKG,C_SCoA,C_Suc,C_FUM,C_MAL,...
    C_OAA, C_NAD, C_ADP_m, Ca_m, C_NADH, C_ATP_m, data);

% Oxidative phosphorylation
[V_He, dC_FLV, V_He_F, dC_NADH, V_Hu, V_H_Leak, dC_ADP_m, V_ANT, ~]...
    = ...
    oxidative_phosphorylation(V_SL, V_IDH, V_KGDH, V_MDH, V_SDH, ...
    delta_Psi_m, C_NADH, C_NAD, C_ATP_m, C_ADP_m, Ca_m, ATP_i, ...
    ADP_i, data);

V_myo = data.V_myo;V_mito=data.V_mito;A_cap=data.A_cap;
dATP_i = V_ANT.*V_mito./V_myo - 0.5.*I_up - (I_NaK + I_Cap).*A_cap./(V_myo.*F) - V_AM;

% Mitochondrial Ca2+ handling processes
[dCa_m, J_uni, J_NaCa] = mitochondrial_Ca2_handling(delta_Psi_m, Ca_m, Ca_i, data);

% Mitochondrial membrane voltage
C_mito = data.C_mito;
ddelta_Psi_m = (V_He + V_He_F - V_Hu - V_ANT - V_H_Leak - J_NaCa - 2.*J_uni)./...
    C_mito;

% Update the state vector
dX=zeros(51,1);

dX(1)   = dV;
dX(2)   = dP_a;
dX(3)   = dP_i;
dX(4)   = dn;
dX(5)   = dr;
dX(6)   = ds_1;
dX(7)   = ds_2;
dX(8)   = ds_3;
dX(9)   = dm;
dX(10)  = dh_1;
dX(11)  = dh_2;
dX(12)  = dd_L;
dX(13)  = df_L;
dX(14)  = dd_T;
dX(15)  = df_T;
dX(16)  = dNa_i;
dX(17)  = dCa_up;
dX(18)  = dCa_rel;
dX(19)  = dCa_i;
dX(20)  = dO_c;
dX(21)  = dO_TnCa;
dX(22)  = dO_TnMgCa;
dX(23)  = dO_TnMgMg;
dX(24)  = dO_Calse;
dX(25)  = dK_o;
dX(26)  = dK_i;
dX(27)  = dF_1;
dX(28)  = dF_2;
dX(29)  = dF_3;
%   dX(30) = dfca; - not used
dX(31)  = dSL;
dX(32)  = dA;
dX(33)  = dTT;
dX(34)  = dU;
dX(35)  = dV_e;
dX(36)  = dATP_i;
dX(37)  = dCa_m;
%   dX(38) = dC_ATP_ic - not used
%   dX(39) = dC_CrP_i  - not used
%   dX(40) = dC_CrP_ic - not used
dX(41)  = dC_ADP_m;
dX(42)  = dC_NADH;
dX(43)  = ddelta_Psi_m;
dX(44)  = dC_ISOC;
dX(45)  = dC_aKG;
dX(46)  = dC_SCoA;
dX(47)  = dC_Suc;
dX(48)  = dC_FUM;
dX(49)  = dC_MAL;
dX(50)  = dC_OAA;
dX(51)  = dC_FLV;

end



