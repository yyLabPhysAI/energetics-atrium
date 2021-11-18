function dX=differential_equations_model_C (t,X,data,f_stim,start_stim,end_stim)
% Names for the state vector values
V      = X(1);
P_a     = X(2);
P_i     = X(3);
n      = X(4);
r      = X(5);
s_1     = X(6);
s_2     = X(7);
s_3     = X(8);
m      = X(9);
h_1     = X(10);
h_2     = X(11);
d_L     = X(12);
f_L     = X(13);
d_T     = X(14);
f_T     = X(15);
Na_i    = X(16);
Ca_up   = X(17);
Ca_rel  = X(18);
Ca_i    = X(19);
O_c    = X(20);
O_TnCa   = X(21);
O_TnMgCa = X(22);
O_TnMgMg = X(23);
O_Calse = X(24);
K_o     = X(25);
K_i     = X(26);
F_1     = X(27);
F_2     = X(28);
F_3     = X(29);
fca    = X(30);
SL     = X(31);
A      = X(32);
TT     = X(33);
U      = X(34);
V_e     = X(35);
ATP_i   = X(36);
Ca_m   = X(37);
C_ATP_ic=X(38);
C_CrP_i= X(39);
C_CrP_ic=X(40);
C_ADP_m= X(41);
C_NADH=  X(42);
delta_Psi_m=X(43);
C_ISOC=  X(44);
C_aKG =  X(45);
C_SCoA=  X(46);
C_Suc=   X(47);
C_FUM=   X(48);
C_MAL=   X(49);
C_OAA=   X(50);
C_FLV=   X(51);


% Unpack the constants to simplefy the code
Is=data.Is;
R=data.R;
T=data.T;
Cm=data.Cm;
F=data.F;
V_i=data.V_i;
V_Ca=data.V_Ca;
V_c=data.V_c;

K_Ca = data.K_Ca;
C_PN=data.C_PN;
C_mito = data.C_mito;
Ff=data.F_f;
MaxATP=data.MaxATP;
K_Na=data.K_Na;r_a=data.r_a;r_c1=data.r_c1;
r_c2=data.r_c2;r_1=data.r_1;r_2=data.r_2;r_3=data.r_3;rho_res=data.rho_res;K_res=data.K_res;
rho_resF=data.rho_resF;psi_B=data.psi_B;K_resF=data.K_resF;r_b=data.r_b;
FADH2=data.FADH2;FADH=data.FADH;p_a=data.p_a;p_b=data.p_b;p_c1=data.p_c1;p_c2=data.p_c2;
p_1=data.p_1;p_2=data.p_2;p_3=data.p_3;KCaATP=data.KCaATP;rho_F1=data.rho_F1;K_F1=data.K_F1;
C_A=data.C_A;V_ant_max=data.V_ant_max;h_ANT=data.h_ANT;g_H=data.g_H;delta_pH=data.delta_pH;
V_myo=data.V_myo;V_mito=data.V_mito;A_cap=data.A_cap;
P_Ca=data.P_Ca;Z_Ca=data.Z_Ca;alpha_m=data.alpha_m;alpha_e=data.alpha_e;V_NC=data.V_NC;Na_e=data.Na_e;
Na_m=data.Na_m;Beta_Ca=data.Beta_Ca;
F_f=data.F_f;
K_M_ATP=data.K_M_ATP;Max_ATP=data.Max_ATP;CATPi=data.CATPi;
K_M_ADP = data.K_M_ADP;
Na_o = data.Na_o; Ca_o = data.Ca_o;

% Nernst potentials
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
I_Cap = I_Cap_max*(Ca_i/(Ca_i + 0.0002));

% L-type Ca(2+) channels
[I_CaL, df_L, dd_L] = L_type_Ca(V, d_L, f_L);

% T-Type calcium channels
[dd_T , df_T, I_CaT] = T_Type_calcium(V, d_T, f_T);

% Background currents
g_NaB = 0.03;
g_CaB = 0.03;

INaB = g_NaB*(V - E_Na);
ICaB = g_CaB*(V - E_Ca);

% Calculate the total current
Itot = I_Kr + I_Ks + I_k1 + I_Kto + I_NaK + I_NaCa + I_Na + INaB + I_CaL + I_CaT + I_Cap + ICaB;

% Sodium and potassium concentration in cytoplasm
dNa_i = (-3*I_NaK - 3*I_NaCa - INaB - I_Na)/(F*V_i);
dK_o = (-2*I_NaK + I_Kr + I_Ks + I_Kto + I_k1)/(F*V_c);
dK_i = (2*I_NaK - I_Kr - I_Ks - I_Kto - I_k1)/(F*V_i);

% Sarcoplasmic Reticulum (SR) and Calcium Handling 
[dO_c, dO_TnCa, dO_TnMgCa, dO_TnMgMg, dO_Calse, phi_ca_i, ...
dCa_up, dCa_rel, dF_1, dF_2, dF_3, dfca, I_up, I_rel] = ...
          SR_calcium_handling(ATP_i, Ca_i, Ca_up, Ca_rel, F_1, F_2, ...
          F_3, O_c, O_TnCa, O_TnMgCa, O_TnMgMg, O_Calse, V, data, fca);

dCa_i = ((2*I_NaCa - I_CaL - I_CaT - I_Cap - ICaB - I_up + I_rel)/(2*V_Ca*F) - phi_ca_i);

%% Mitochondrial energy metabolism, Ca2+ dynamics and oxygen consumption

% Mitochondrial energetics and EC coupling
ADP_i = C_A - ATP_i;

% Force generation
[dSL, dA, dTT, dU, dV_e, Force] = ...
force_generation(V_e, SL, TT, A, U, data, Ca_i, ATP_i, ADP_i);

%energy metabolites reaction rates
C_ATP_m = C_A - C_ADP_m;

% Oxidation state
C_NAD = C_PN - C_NADH;

% The TCA cycle
[dC_ISOC, dC_aKG, dC_SCoA, dC_Suc, dC_FUM,...
    dC_MAL, dC_OAA, V_SL, V_IDH, V_KGDH, V_MDH, V_SDH] ...
    = TCA_cycle(C_ISOC,C_aKG,C_SCoA,C_Suc,C_FUM,C_MAL,...
    C_OAA, C_NAD, C_ADP_m, Ca_m, C_NADH, C_ATP_m, data);

%C_ADP_ic=C_A-C_ATP_ic;
A_res = ((R*T)/(F))*log(K_res*sqrt(C_NADH./C_NAD));
delta_mu_h = -2.303*(R*T/F)*delta_pH + delta_Psi_m;

A_res_F = ((R*T)/(F))*log(K_resF*sqrt(FADH2/FADH));
A_F1 = ((R*T)/F)*log((K_F1*C_ATP_m)/(C_ADP_m*P_i));

V_He = 6*rho_res*(r_a*exp((A_res*F)/(R*T)) - (r_a + r_b)*exp((12*F*delta_mu_h)/(R*T)))/...
                 ((1 + r_1*exp((A_res*F)/(R*T))) * exp((6*F*psi_B)/(R*T)) + (r_2 + r_3 * exp((A_res*F)/(R*T)))* exp((12*F*delta_mu_h)/(R*T)));

V_He_F = 4*rho_resF*(r_a*exp((A_res_F*F)/(R*T)) - (r_a + r_b)*exp((12*F*delta_mu_h)/(R*T)))/...
                    ((1+r_1* exp((A_res_F*F)/(R*T)))* exp((6*F*psi_B)/(R*T))+(r_2+r_3* exp((A_res_F*F)/(R*T)))* exp((12*F*delta_mu_h)/(R*T)));

Juni = P_Ca*((Z_Ca*delta_Psi_m*F)/(R*T))*(alpha_m*Ca_m*exp((-Z_Ca*delta_Psi_m*F)/(R*T))-(alpha_e*Ca_i))/...
                                         (exp((-Z_Ca*delta_Psi_m*F)/(R*T))-1);

Jnc = V_NC*(exp(0.5*delta_Psi_m*F/(R*T))*(Na_e^3*Ca_m)/(K_Na^3*K_Ca) - exp(-0.5*delta_Psi_m*F/(R*T))*(Na_m^3*Ca_i)/(K_Na^3*K_Ca))/...
           (1 + (Na_e^3/K_Na^3) + (Ca_m/K_Ca) + (Na_e^3*Ca_m)/(K_Na^3*K_Ca) + (Na_m^3/K_Na^3) + (Ca_i/K_Ca) + (Na_m^3*Ca_i)/(K_Na^3*K_Ca));

V_Hu = -3*rho_F1*(10^2*p_a*(1+exp(A_F1/(R*T)))-(p_a+p_b)*exp((3*F*delta_mu_h)/(R*T)))/...
                  ((1 + p_1*exp((A_F1*F)/(R*T)))*exp((3*F*psi_B)/(R*T)) + (p_2 + p_3*exp((A_F1*F)/(R*T)))*exp((3*F*delta_mu_h)/(R*T)));

V_H_Leak = g_H*delta_mu_h;

V_ANT = V_ant_max*(0.75*(1 - ((0.25*ATP_i*0.45*C_ADP_m)/(0.17*ADP_i*0.025*C_ATP_m)))*exp(-(F*delta_Psi_m)/(R*T)))/...
                  ((1 + ((0.25*ATP_i)/(0.225*ADP_i))*exp((-h_ANT*F*delta_Psi_m)/(R*T)))*(1 + ((0.45*C_ADP_m)/(0.025*C_ATP_m ))));


Force_ATP = 1.02/...
            (1 + (K_M_ATP/ATP_i)*(1 + (CATPi - ATP_i)/K_M_ADP));

V_AM = Max_ATP*F_f*A*Force_ATP;
ATP_XB = MaxATP*Ff*A.*Force; 

dATPi = V_ANT*(V_mito/V_myo) - (I_NaK+I_Cap)*A_cap/(V_myo*F) - 0.5*I_up - V_AM;


V_ATPase = -rho_F1*((10^2*p_a + p_c1*exp((3*F*psi_B)/(R*T)))*exp(F*A_F1/(R*T)) - (p_a*exp((3*F*delta_mu_h)/(R*T)) + p_c2*exp(F*A_F1/(R*T))*exp((3*F*delta_mu_h)/(R*T))))/...
                   ((1 + p_1*exp((A_F1*F)/(R*T)))*exp((3*F*psi_B)/(R*T)) + (p_2 + p_3*exp((A_F1*F)/(R*T)))*exp((3*F*delta_mu_h)/(R*T)))*(1 - exp(-(Ca_m/KCaATP)));

dC_ADP_m = V_ANT - V_ATPase - V_SL;

V_O2 = 0.5*rho_res*((r_a + r_c1*exp((6*F*psi_B)/(R*T)))*exp((A_res*F)/(R*T)) - r_a*exp((12*F*delta_mu_h)/(R*T)) + r_c2*exp((A_res*F)/(R*T))*exp((12*F*delta_mu_h)/(R*T)))/...
                   ((1 + r_1*exp((A_res*F)/(R*T)))*exp((6*F*psi_B)/(R*T)) + (r_2 + r_3*exp((A_res*F)/(R*T)))*exp((12*F*delta_mu_h)/(R*T)));


dC_NADH = -V_O2 + V_IDH + V_KGDH + V_MDH;


dC_FLV = V_SDH - V_O2;

dCa_m = Beta_Ca.*(Juni - Jnc);

ddelta_Psi_m = (V_He + V_He_F - V_Hu - V_ANT - V_H_Leak - Jnc - 2*Juni)/...
                C_mito;

%%
impulseFactor = 100;

if ((mod(t,1/f_stim) >= start_stim && mod(t,1/f_stim) < 0.001*impulseFactor + start_stim) && (t <= end_stim))
    I_stim=Is/impulseFactor;
else 
    I_stim=0;
end

dV = -(Itot + I_stim)/...
       Cm;



% Update the state vector
dX=zeros(51,1);

dX(1)=(dV);
dX(2)=(dP_a);
dX(3)=(dP_i);
dX(4)=(dn);
dX(5)=(dr);
dX(6)=(ds_1);
dX(7)=(ds_2);
dX(8)=(ds_3);
dX(9)=(dm);
dX(10)=(dh_1);
dX(11)=(dh_2);
dX(12)=(dd_L);
dX(13)=(df_L);
dX(14)=(dd_T);
dX(15)=(df_T);
dX(16)=(dNa_i);
dX(17)=(dCa_up);
dX(18)=(dCa_rel);
dX(19)=(dCa_i);
dX(20)=(dO_c);
dX(21)=(dO_TnCa);
dX(22)=(dO_TnMgCa);
dX(23)=(dO_TnMgMg);
dX(24)=(dO_Calse);
dX(25)=(dK_o);
dX(26)=(dK_i);
dX(27)=(dF_1);
dX(28)=(dF_2);
dX(29)=(dF_3);
dX(30)=(dfca);
dX(31)=(dSL);
dX(32)=(dA);
dX(33)=(dTT);
dX(34)=(dU);
dX(35)=(dV_e);
dX(36)=(dATPi);
dX(37)=(dCa_m);
%dX(38)=(dC_ATP_ic);
%dX(39)=(dC_CrP_i);
%dX(40)=(dC_CrP_ic);
dX(41)=(dC_ADP_m);
dX(42)=(dC_NADH);
dX(43)=(ddelta_Psi_m);
dX(44)=(dC_ISOC);
dX(45)=(dC_aKG);
dX(46)=(dC_SCoA);
dX(47)=(dC_Suc);
dX(48)=(dC_FUM);
dX(49)=(dC_MAL);
dX(50)=(dC_OAA);
dX(51)=(dC_FLV);
end



