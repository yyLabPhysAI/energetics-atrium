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
h1     = X(10);
h2     = X(11);
dL     = X(12);
fL     = X(13);
dT     = X(14);
fT     = X(15);
Na_i    = X(16);
Caup   = X(17);
Carel  = X(18);
Ca_i    = X(19);
fac    = X(20);
faTc   = X(21);
faTmgc = X(22);
faTmgm = X(23);
faCalse= X(24);
K_o     = X(25);
K_i     = X(26);
F1     = X(27);
F2     = X(28);
F3     = X(29);
fca    = X(30);
SL     = X(31);
A      = X(32);
TT     = X(33);
U      = X(34);
Ve     = X(35);
ATPi   = X(36);
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
RTONF=data.RTONF;
V_i=data.V_i;
V_Ca=data.V_Ca;
V_c=data.V_c;
V_up=data.V_up;
V_rel=data.V_rel;

Kcy_ca=data.Kcy_ca;Kxcs=data.Kxcs;
Ksr_ca=data.Ksr_ca;Ff=data.Ff;SL_0=data.SL_0;
F_k05=data.F_k05;F_go=data.F_go;F_XB=data.F_XB;MaxATP=data.MaxATP;
k_XBATP=data.k_XBATP;
k_XBADP=data.k_XBADP;
K_Na=data.K_Na;r_a=data.r_a;r_c1=data.r_c1;
r_c2=data.r_c2;r_1=data.r_1;r_2=data.r_2;r_3=data.r_3;rho_res=data.rho_res;K_res=data.K_res;
rho_resF=data.rho_resF;psi_B=data.psi_B;K_resF=data.K_resF;r_b=data.r_b;
FADH2=data.FADH2;FADH=data.FADH;p_a=data.p_a;p_b=data.p_b;p_c1=data.p_c1;p_c2=data.p_c2;
p_1=data.p_1;p_2=data.p_2;p_3=data.p_3;KCaATP=data.KCaATP;rho_F1=data.rho_F1;K_F1=data.K_F1;
C_A=data.C_A;V_ant_max=data.V_ant_max;h_ANT=data.h_ANT;g_H=data.g_H;delta_pH=data.delta_pH;
C_PN=data.C_PN;C_mito=data.C_mito;C_AcCoA=data.C_AcCoA;k_cat_cs=data.k_cat_cs;
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
k_ASP=data.k_ASP;V_myo=data.V_myo;V_mito=data.V_mito;A_cap=data.A_cap;
P_Ca=data.P_Ca;Z_Ca=data.Z_Ca;alpha_m=data.alpha_m;alpha_e=data.alpha_e;V_NC=data.V_NC;Na_e=data.Na_e;
Na_m=data.Na_m;Beta_Ca=data.Beta_Ca;K_D_Ca=data.K_D_Ca;K_D_Mg=data.K_D_Mg;N_c=data.N_c;
F_k0=data.F_k0;F_k1=data.F_k1;FN=data.FN;F_kl=data.F_kl;F_f=data.F_f;
F_gl=data.F_gl;K_M_ATP=data.K_M_ATP;Max_ATP=data.Max_ATP;CATPi=data.CATPi;
K_M_ADP = data.K_M_ADP;
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
if(abs(V + 44.4) < 0.0001)
    Am = 460*12.673;
else
    Am = -460*(V + 44.4)/(exp(-(V + 44.4)/12.673) - 1);
end

Bm = 18400*exp(-(V + 44.4)/12.673);
dm = Am*(1 - m) - Bm*m;

Ah = 44.9*exp(-(V + 66.9)/5.57);
Bh = 1491.0/(1 + 323.3*exp(-(V + 94.6)/12.9));

th1 = 0.03/(1 + exp((V+40)/6)) + 0.00015; %0.00035 - Lindblad
th2 = 0.12/(1 + exp((V+60)/2)) + 0.00045;  %0.00295 - Lindblab

hm  = Ah/(Ah + Bh);
dh1 = (hm - h1)/th1;
dh2 = (hm - h2)/th2;

if(abs(V) > 0.0001)
    INa = 0.0014*m^3*(0.635*h1 + 0.365*h2)*140*V...
        *(F/RTONF)*(exp((V-E_Na)/RTONF) - 1)/(exp(V/RTONF) - 1); 
else
    INa = 0.0014*m^3*(0.635*h1 + 0.365*h2)*140*(F)*(exp((V-E_Na)/RTONF)-1);     
end

% calculate icap
ICap = 9.509*(Ca_i/(Ca_i + 0.0002));

%% L-type

% calculate ical
V_CaL = V + 10;

Adl = -16.72*(V_CaL + 35)/(exp(-(V_CaL + 35)/2.5)-1) - 50.0*V_CaL/(exp(-V_CaL/4.808) - 1);
if(abs(V_CaL+35) < 0.0001)
     Adl = 16.72*2.5 - 50.0*V_CaL/(exp(-V_CaL/4.808) - 1);
end
if(abs(V_CaL) < 0.0001)
     Adl = -16.72*(V_CaL + 35)/(exp(-(V + 35)/2.5) - 1) + 50*4.808;
end

if(abs(V_CaL - 5) < 0.0001)
    Bdl = 4.48*2.5;
else
    Bdl = 4.48*(V_CaL-5)/(exp((V_CaL - 5)/2.5) - 1);
end

tdl = 1/(Adl + Bdl);

dlm = 1/(1 + exp(-(V_CaL+0.95)/6.6));
ddL=(dlm - dL)/tdl;

V_CaL = V_CaL - 10 - 10; % updated JB

if(abs(V_CaL + 28) < 0.0001)
    Afl = 8.49*4;
else
    Afl = 8.49*(V_CaL + 28)/(exp((V_CaL + 28)/4) - 1);
end

Bfl = 67.922/(1 + exp(-(V_CaL + 28)/4));

flm = Afl/(Afl + Bfl);
tfl = 1/(Afl + Bfl);
dfL = (flm - fL)/tfl;

ICaL = 2.1*4*(dL*fL + 1.0/(1 + exp(-(V - 23.0)/12.0)))*(V - 50);      % CT - 1.8, PM - 2.1

% calculace icat
Adt = 674.173*exp((V + 23.3)/30);
Bdt = 674.173*exp(-(V + 23.3)/30);
tdt = 1/(Adt + Bdt);
dtm = 1/(1 + exp(-(V + 23)/6.1));
ddT = (dtm - dT)/tdt;

Aft = 9.637*exp(-(V + 75)/83.3);
Bft = 9.637*exp((V + 75)/15.38);
tft = 1/(Aft + Bft);
ftm = Aft/(Aft + Bft);
dfT = (ftm - fT)/tft;

ICaT = 6*dT*fT*(V - 38);

% calculate ibna
IbNa = 0.03*(V - E_Na);  % 0.02 - CT, 0.03 - PM

% calculate ibca
IbCa = 0.03*(V - E_Ca);  % 0.02 - CT, 0.03 - PM


% Calculate the total current
Itot = I_Kr + I_Ks + I_k1 + I_Kto + I_NaK + I_NaCa + INa + IbNa + ICaL + ICaT + ICap + IbCa;

%% LINDBLAD Ca++ handling - SR Ca++

% Na+ concentration change
dNai = (-3*I_NaK-3*I_NaCa-IbNa-INa)/(F*V_i);

% Uptake current of Ca++ to the SR
Iup = (2800/1000)*(ATPi/7.977)*(((Ca_i/Kcy_ca) - (Kxcs*Kxcs*Caup/Ksr_ca))/(((Ca_i + Kcy_ca)/Kcy_ca) + (Kxcs*(Caup + Ksr_ca)/Ksr_ca)));
% Release current of Ca++ from the SR
Irel = 200000*((F2/(F2 + 0.25))^2)*(Carel - Ca_i);

% transfer current of Ca++ between the SR compartments
Itr = (Caup - Carel)*2*F*V_up/0.01;

% calculation of the Ca++ concentration in the different compartments
dfac     = 200000*Ca_i*(1 - fac) - 476*fac;
dfaTc    = 78400*Ca_i*(1 - faTc) - 392*faTc;
dfaTmgc  = 200000*Ca_i*(1 - faTmgc - faTmgm) - 6.6*faTmgc;
dfaTmgm  = 2000*2.5*(1 - faTmgc - faTmgm) - 666*faTmgm;
dfaCalse = 480*Carel*(1 - faCalse) - 400*faCalse;
dfab     = 0.08*dfaTc + 0.16*dfaTmgc + 0.045*dfac;

Faraday = F;

dCaup = (Iup - Itr)/(2*Faraday*V_up);
dCarel = ((Itr - Irel)/(2*Faraday*V_rel) - 31*(480*Carel*(1 - faCalse) - 400*faCalse));
dCai =((2*I_NaCa - ICaL - ICaT - ICap - IbCa - Iup + Irel)/(2*V_Ca*Faraday) - dfab);

% calculation of the K+ concentration in the different compartments
dKc = (-2*I_NaK + I_Kr + I_Ks + I_Kto + I_k1)/(Faraday*V_c);
dKi = (2*I_NaK - I_Kr - I_Ks - I_Kto - I_k1)/(Faraday*V_i);

% opening and closing of SR Ca++ release channels
ract  = 240*exp((V - 20)/12.5) + 203.8*(Ca_i/(Ca_i+0.0003))^4;
rinact= 33.96 + 339.6*(Ca_i/(Ca_i+0.0003))^4;

dF1 = 0.815*F3 - ract*F1;
dF2 = ract*F1 - rinact*F2;
dF3 = rinact*F2 - 0.815*F3;

% variable not really used in the paper - implemented for future use
dfca = 0*fca;

%% Force generation

dSL = -Ve;
N_XB = (1e-6)*(SL - SL_0)*N_c*(TT + U)*1000/2;
K_Ca = F_k0 + F_k1*(N_XB^FN)/(F_k05^FN + N_XB^FN);
K_minus1 = F_kl/K_Ca;
dA = F_kl*Ca_i*(1 - A - TT - U) - A*(Ff + K_minus1) + TT*(F_go + F_gl*Ve);
dTT = Ff*A - TT*(F_go + F_gl*Ve + K_minus1) + F_kl*Ca_i*U;
dU = K_minus1*TT - (F_go + F_gl*Ve + F_kl*Ca_i)*U;

ADPi=C_A-ATPi;
N_XB_ATP=1.02/(1+(k_XBATP/ATPi)*(1+ADPi/k_XBADP));
Force=F_XB*N_XB_ATP*((SL-SL_0)/2)*(TT+U)*N_c;

dVe = 0; % We assume here isometric contraction, change here for moving models

%% Energy consumption

%energy metabolites reaction rates
C_ATP_m = C_A - C_ADP_m;

%C_ADP_ic=C_A-C_ATP_ic;
C_NAD = C_PN - C_NADH;
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

V_ANT = V_ant_max*(0.75*(1 - ((0.25*ATPi*0.45*C_ADP_m)/(0.17*ADPi*0.025*C_ATP_m)))*exp(-(F*delta_Psi_m)/(R*T)))/...
                  ((1 + ((0.25*ATPi)/(0.225*ADPi))*exp((-h_ANT*F*delta_Psi_m)/(R*T)))*(1 + ((0.45*C_ADP_m)/(0.025*C_ATP_m ))));


Force_ATP = 1.02/...
            (1 + (K_M_ATP/ATPi)*(1 + (CATPi - ATPi)/K_M_ADP));

V_AM = Max_ATP*F_f*A*Force_ATP;
ATP_XB = MaxATP*Ff*A.*Force; 

dATPi = V_ANT*(V_mito/V_myo) - (I_NaK+ICap)*A_cap/(V_myo*F) - 0.5*Iup - V_AM;


V_ATPase = -rho_F1*((10^2*p_a + p_c1*exp((3*F*psi_B)/(R*T)))*exp(F*A_F1/(R*T)) - (p_a*exp((3*F*delta_mu_h)/(R*T)) + p_c2*exp(F*A_F1/(R*T))*exp((3*F*delta_mu_h)/(R*T))))/...
                   ((1 + p_1*exp((A_F1*F)/(R*T)))*exp((3*F*psi_B)/(R*T)) + (p_2 + p_3*exp((A_F1*F)/(R*T)))*exp((3*F*delta_mu_h)/(R*T)))*(1 - exp(-(Ca_m/KCaATP)));

V_SL = k_f_SL*((C_SCoA*C_ADP_m - (C_Suc*C_ATP_m*C_CoA)/(k_E_SL )));

dC_ADP_m = V_ANT - V_ATPase - V_SL;

V_O2 = 0.5*rho_res*((r_a + r_c1*exp((6*F*psi_B)/(R*T)))*exp((A_res*F)/(R*T)) - r_a*exp((12*F*delta_mu_h)/(R*T)) + r_c2*exp((A_res*F)/(R*T))*exp((12*F*delta_mu_h)/(R*T)))/...
                   ((1 + r_1*exp((A_res*F)/(R*T)))*exp((6*F*psi_B)/(R*T)) + (r_2 + r_3*exp((A_res*F)/(R*T)))*exp((12*F*delta_mu_h)/(R*T)));

f_a_IDH = ((1 + (C_ADP_m/K_ADP_a))*(1 + (Ca_m/K_Ca_a)))^(-1);

f_i_IDH = 1 + (C_NADH/K_i_NADH);

V_IDH = k_cat_IDH*E_T_IDH*(1 + (C_H/k_h_1) + (k_h_2/C_H) + f_i_IDH*(K_M_NAD/C_NAD) + f_a_IDH*(K_M_ISOC/C_ISOC)^2 + f_a_IDH*f_i_IDH*(K_M_ISOC/C_ISOC)^2*(K_M_NAD/C_NAD))^-1;

f_a_KGDH = ((1 + (C_Mg/K_D_Mg))*(1 + (Ca_m/K_D_Ca)))^(-1);

V_KGDH = (k_cat_KGDH*E_T_KGDH)/(1 + f_a_KGDH*((K_M_aKG)/C_aKG)^(n_aKG)+f_a_KGDH*((K_M_NAD_new)/C_NAD));

f_h_a = 1/...
        (1 + C_H/k_h1 + C_H^2/(k_h1*k_h2)) + k_offset;

f_h_i = (1 + k_h3/C_H + (k_h3*k_h4)/(C_H^2))^-2;

V_MDH = (k_cat_MDH*E_T_MDH *f_h_a*f_h_i)/(1 + (K_M_MAL/C_MAL)*(1 + (C_OAA/K_i_OAA)) + (K_M_NAD_mdh/C_NAD) + ((K_M_MAL/C_MAL)*(1 + (C_OAA/K_i_OAA))*(K_M_NAD_mdh))/(C_NAD));

dC_NADH = -V_O2 + V_IDH + V_KGDH + V_MDH;

C_CIT = C_K_int - (C_ISOC + C_aKG + C_SCoA + C_Suc + C_FUM + C_MAL + C_OAA);

V_ACO = k_f_ACO*(C_CIT - (C_ISOC/k_E_ACO));

dC_ISOC = V_ACO - V_IDH;

V_AAT = k_f_AAT*C_OAA*C_GLU*(k_ASP*K_E_AAT)/...
                            (k_ASP*K_E_AAT+C_aKG*k_f_AAT);

dC_aKG = V_IDH - V_KGDH + V_AAT;

dC_SCoA = V_KGDH - V_SL;

V_SDH = (k_cat_SDH*E_T_SDH)/...
        (1 + ((K_M_SUC)/C_Suc)*(1 + (C_OAA/K_i_sdh_OAA))*((1 + (C_FUM/(K_i_FUM)))));

dC_Suc = V_SL - V_SDH;

V_FH = k_f_FH*(C_FUM - C_MAL/(K_E_FH));

dC_FUM = V_SDH - V_FH;

dC_MAL = V_FH - V_MDH;

V_CS = k_cat_cs*E_T_cs*(1 + (K_M_AcCoA/C_AcCoA) + (K_M_OAA/C_OAA) + (K_M_AcCoA/C_AcCoA) + (K_M_OAA/C_OAA))^(-1);

dC_OAA = V_MDH - V_CS - V_AAT;

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
dX(10)=(dh1);
dX(11)=(dh2);
dX(12)=(ddL);
dX(13)=(dfL);
dX(14)=(ddT);
dX(15)=(dfT);
dX(16)=(dNai);
dX(17)=(dCaup);
dX(18)=(dCarel);
dX(19)=(dCai);
dX(20)=(dfac);
dX(21)=(dfaTc);
dX(22)=(dfaTmgc);
dX(23)=(dfaTmgm);
dX(24)=(dfaCalse);
dX(25)=(dKc);
dX(26)=(dKi);
dX(27)=(dF1);
dX(28)=(dF2);
dX(29)=(dF3);
dX(30)=(dfca);
dX(31)=(dSL);
dX(32)=(dA);
dX(33)=(dTT);
dX(34)=(dU);
dX(35)=(dVe);
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



