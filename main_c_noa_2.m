clear 
clc

%% Constant decleration and loading
model_constants_energetics
F_XB = data.F_XB;
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
%% Choice of stimulation parameters
start_stim = 0; % [s] when to start the current stimulation
end_stim   = 100;    % [s] when to end it
f_stim     = 1;    % [Hz] at what frequency do you want the current pulses to be?
TMAX       = 100;    % [s] until what time to calculate the simulation?

%% Run the model
[t,x] = ode23tb(@(t,x)differential_equations_model_C(t,x,data,f_stim,start_stim,end_stim), [0 TMAX] ,y0',options);
V = x(:,1);
P_a = x(:,2);
P_i = x(:,3);
n = x(:,4);
r = x(:,5);
s_1 = x(:,6);
s_2 = x(:,7);
s_3 = x(:,8);
m = x(:,9);
h_1 = x(:,10);
h_2 = x(:,11);
d_L = x(:,12);
f_L = x(:,13);
d_T = x(:,14);
f_T = x(:,15);
Na_i = x(:,16);
Ca_up = x(:,17);
Ca_rel = x(:,18);
Ca_i = x(:,19);
O_c = x(:,20);
O_TnCa = x(:,21);
O_TnMgCa = x(:,22);
O_TnMgMg = x(:,23);
O_Calse = x(:,24);
K_o = x(:,25);
K_i = x(:,26);
F_1 = x(:,27);
F_2 = x(:,28);
F_3 = x(:,29);
fca = x(:,30);
SL = x(:,31);
A = x(:,32);
TT = x(:,33);
U = x(:,34);
Ve = x(:,35);
ATP_i = x(:,36);
Ca_m = x(:,37);
C_ATP_ic = x(:,38);
C_CrP_i = x(:,39);
C_CrP_ic = x(:,40);
C_ADP_m = x(:,41);
C_NADH = x(:,42);
delta_Psi_m = x(:,43);
C_ISOC = x(:,44);
C_aKG = x(:,45);
C_SCoA = x(:,46);
C_Suc = x(:,47);
C_FUM = x(:,48);
C_MAL = x(:,49);
C_OAA = x(:,50);
C_FLV = x(:,51);

Ek   = data.RTONF*log(K_o./K_i);
ENa  = data.RTONF*log(140./Na_i);
ECa  = (data.RTONF/2)*log(2.5./Ca_i);

% calculation of ikf
IKf = 3.5*P_a.*P_i.*(V-Ek);

Apa = 9.0*exp(V/25.371);
Bpa = 1.3*exp(-V./13.026);
pam = 1./(1+exp(-(V+5.1)./7.4));
tpa = 1./(Apa+Bpa);
dPa = (pam-P_a)./tpa;

Api = 100*exp(-V./54.645);
Bpi = 656*exp(V./106.157);
pim = 1./(1+exp((V+47.3921)./18.6603));
tpi = 1./(Api+Bpi);
dPi = (pim-P_i)./tpi;

% calculation of iks
IKs = 2.5*n.*(V-Ek);

An = 1.66*exp(V./69.452);
Bn = 0.3*exp(-V./21.826);
nm = 1./(1+exp(-(V-0.9)./13.8));
tn = 1./(An+Bn) + 0.06;
dn = (nm-n)./tn;

% calculation of ik1
E0 = V-Ek+3.6;

Ik1 = 2.5*5.08*(K_o./(K_o+0.59)).^3.*(V-Ek)./(1+exp(1.393*E0./data.RTONF)); %CT - 2.0, PM - 2.5

% calculation of ito
Ito = 0.35*50.02*r.*(0.590*s_1.^3 + 0.410*s_2.^3).*(0.600*s_3.^6 + 0.4).*(V - Ek); % CT - 0.2, PM - 0.35

Isus = 2.4*(V+70);  %CT - 1.4, PM - 2.4
Ito  = Ito + Isus;

Ar = 386.6*exp(V./12);
Br = 8.011*exp(-V./7.2);
rm = 1./(1 + exp(-(V+15)./5.633));
tr = 1./(Ar+Br) + 0.0004;
dr = (rm-r)./tr;

s1m = 1./(1 + exp((V+28.29)./7.06));
ts1 = 0.5466./(1 + exp((V+32.8)./0.1))+0.0204;
ds1 = (s1m-s_1)./ts1;

s2m = 1./(1 + exp((V+28.29)./7.06));
ts2 = 5.75./(1 + exp((V+32.8)./0.1)) + 0.45./(1 + exp(-(V-13.54)./13.97));
ds2 = (s2m-s_2)./ts2;

s3m = ((1./(1 + exp((V+50.67)./27.38))) + 0.666)./1.666;
ts3 = (7.5./(1 + exp((V+23.0)./0.5))) + 0.5;
ds3 = (s3m-s_3)./ts3;

% calculate ip <------- This is the Na-K ATPase current!!!!!!!!!!!!!
Ip = 64.41*(K_o./(K_o+1)).*(Na_i.^1.5./(Na_i.^1.5 + 11.^1.5)).*(1.6./(1.5 + exp(-(V+60)./40)));

% calculate inaca
INaCa = 0.02*(Na_i.^3*2.5.*exp(0.450*V./data.RTONF) - 140^3*Ca_i.*exp(V.*(0.45-1)./data.RTONF))./(1 + 0.0003*(Ca_i*140^3 + 2.5*Na_i.^3));



%calculate ina
if(abs(V+44.4) < 0.0001)
    Am = 460*12.673;
else
    Am = -460*(V + 44.4)./(exp(-(V+44.4)./12.673) - 1);
end

Bm = 18400*exp(-(V+44.4)./12.673);
dm = Am.*(1-m) - Bm.*m;

Ah = 44.9*exp(-(V+66.9)./5.57);
Bh = 1491.0./(1+323.3*exp(-(V+94.6)./12.9));

th1 = 0.03./(1 + exp((V+40)/6)) + 0.00015; %0.00035 - Lindblad
th2 = 0.12./(1 + exp((V+60)/2)) + 0.00045;  %0.00295 - Lindblab

hm  = Ah./(Ah+Bh);
dh1 = (hm-h_1)./th1;
dh2 = (hm-h_2)./th2;

if(abs(V) > 0.0001)
    INa = 0.0014.*m.^3.*(0.635.*h_1 + 0.365.*h_2)*140.*V...
        .*(data.F/data.RTONF).*(exp((V-ENa)./data.RTONF) - 1)./(exp(V./data.RTONF) - 1); 
else
    INa = 0.0014*m.^3.*(0.635*h_1 + 0.365*h_2).*140*(data.F).*(exp((V-ENa)./data.RTONF)-1);     
end

% calculate icap
ICap = 9.509*(Ca_i./(Ca_i + 0.0002));



% calculate ical
V_CaL = V + 10;

Adl = -16.72*(V_CaL + 35)./(exp(-(V_CaL+35)./2.5)-1) - 50.0*V_CaL./(exp(-V_CaL./4.808) - 1);
if(abs(V_CaL+35) < 0.0001)
     Adl = 16.72*2.5 - 50.0*V_CaL./(exp(-V_CaL./4.808) - 1);
end
if(abs(V_CaL) < 0.0001)
     Adl = -16.72*(V_CaL + 35)./(exp(-(V+35)./2.5)-1) + 50*4.808;
end

if(abs(V_CaL-5) < 0.0001)
    Bdl = 4.48*2.5;
else
    Bdl = 4.48*(V_CaL-5)./(exp((V_CaL-5)./2.5) - 1);
end

tdl = 1./(Adl + Bdl);

dlm = 1./(1+exp(-(V_CaL+0.95)./6.6));
ddL=(dlm-d_L)./tdl;

V_CaL = V_CaL-10-10; % updated JB

if(abs(V_CaL+28) < 0.0001)
    Afl = 8.49*4;
else
    Afl = 8.49*(V_CaL + 28)./(exp((V_CaL + 28)./4) - 1);
end

Bfl = 67.922./(1 + exp(-(V_CaL + 28)./4));

flm = Afl./(Afl + Bfl);
tfl = 1./(Afl + Bfl);
dfL = (flm - f_L)./tfl;

ICaL = 2.1*4*(d_L.*f_L + 1.0./(1 + exp(-(V - 23.0)./12.0))).*(V - 50);      % CT - 1.8, PM - 2.1

% calculace icat
Adt = 674.173*exp((V + 23.3)./30);
Bdt = 674.173*exp(-(V + 23.3)./30);
tdt = 1./(Adt + Bdt);
dtm = 1./(1 + exp(-(V + 23)./6.1));
ddT = (dtm - d_T)./tdt;

Aft = 9.637*exp(-(V + 75)./83.3);
Bft = 9.637*exp((V + 75)./15.38);
tft = 1./(Aft + Bft);
ftm = Aft./(Aft + Bft);
dfT = (ftm - f_T)./tft;

ICaT = 6*d_T.*f_T.*(V - 38);

% calculate ibna
IbNa = 0.03*(V - ENa);  % 0.02 - CT, 0.03 - PM

% calculate ibca
IbCa = 0.03*(V  - ECa);  % 0.02 - CT, 0.03 - PM


% Calculate the total current
Itot = IKf + IKs + Ik1 +Ito+Ip+INaCa+INa+IbNa+ICaL+ICaT+ICap+IbCa;


% Na+ concentration change
dNai = (-3*Ip-3*INaCa-IbNa-INa)./(data.F*data.V_i);

% Uptake current of Ca++ to the SR
Iup = (2800/1000).*(ATP_i/7.977).*(((Ca_i./data.K_cyCa) - (data.K_xcs*data.K_xcs*Ca_up./data.K_srCa))./(((Ca_i + data.K_cyCa)./data.K_cyCa) + (data.K_xcs.*(Ca_up + data.K_srCa)./data.K_srCa)));
% Release current of Ca++ from the SR
Irel = 200000*((F_2./(F_2 + 0.25)).^2).*(Ca_rel - Ca_i);

% transfer current of Ca++ between the SR compartments
Itr = (Ca_up - Ca_rel)*2*data.F*data.V_up./0.01;

% calculation of the Ca++ concentration in the different compartments
dfac     = 200000*Ca_i.*(1 - O_c) - 476*O_c;
dfaTc    = 78400*Ca_i.*(1 - O_TnCa) - 392*O_TnCa;
dfaTmgc  = 200000*Ca_i.*(1 - O_TnMgCa - O_TnMgMg) - 6.6*O_TnMgCa;
dfaTmgm  = 2000*2.5*(1 - O_TnMgCa - O_TnMgMg) - 666*O_TnMgMg;
dfaCalse = 480*Ca_rel.*(1 - O_Calse) - 400*O_Calse;
dfab     = 0.08*dfaTc + 0.16*dfaTmgc + 0.045*dfac;

Faraday = data.F;
Vup     = data.V_up;
Vrel    = data.V_rel;
VCa     = data.V_Ca;
Vc      = data.V_c;
Vi      = data.V_i;

dCaup = (Iup - Itr)./(2*Faraday*Vup);
dCarel = ((Itr - Irel)./(2*Faraday.*Vrel) - 31*(480*Ca_rel.*(1 - O_Calse) - 400*O_Calse));
dCai =((2*INaCa - ICaL - ICaT - ICap - IbCa - Iup + Irel)./(2*VCa*Faraday) - dfab);

% calculation of the K+ concentration in the different compartments
dKc = (-2*Ip + IKf + IKs + Ito + Ik1)./(Faraday*Vc);
dKi = (2*Ip - IKf - IKs - Ito - Ik1)./(Faraday*Vi);

% opening and closing of SR Ca++ release channels
ract  = 240*exp((V - 20)./12.5) + 203.8*(Ca_i./(Ca_i+0.0003)).^4;
rinact= 33.96 + 339.6*(Ca_i./(Ca_i+0.0003)).^4;

dF1 = 0.815*F_3 - ract.*F_1;
dF2 = ract.*F_1 - rinact.*F_2;
dF3 = rinact.*F_2 - 0.815.*F_3;

% variable not really used in the paper - implemented for future use
dfca = 0*fca;
% %%


SL_0 = data.SL_0;
N_c = data.N_c;
F_k0 = data.F_k0;
F_k1 = data.F_k1;
F_k05 = data.F_k05;
F_kl = data.F_kl;
Ff = data.F_f;
F_go = data.F_go; 
F_gl = data.F_gl;
FN = data.FN;

dSL = -Ve;
N_XB = (1e-6)*(SL - SL_0).*N_c.*(TT + U).*1000/2;
K_Ca = F_k0 + F_k1.*(N_XB.^FN)./(F_k05.^FN + N_XB.^FN);
K_minus1 = F_kl./K_Ca;
dA = F_kl.*Ca_i.*(1 - A - TT - U) - A.*(Ff + K_minus1) + TT.*(F_go + F_gl.*Ve);
dTT = Ff.*A - TT.*(F_go + F_gl.*Ve + K_minus1) + F_kl.*Ca_i.*U;
dU = K_minus1.*TT - (F_go + F_gl.*Ve + F_kl.*Ca_i).*U;
ADPi=data.C_A-ATP_i;
N_XB_ATP=1.02./(1+(data.k_XBATP./ATP_i).*(1+ADPi./data.k_XBADP));
Force=data.F_XB.*N_XB_ATP.*((SL-SL_0)./2).*(TT+U).*data.N_c;

dVe = 0; % We assume here isometric contraction, change here for moving models


%dATPi = 0*ATP_i; % assuming const. ATP concentration  - implemented for future use
%energy metabolites reaction rates
%ATP_i=7.977;
C_ATP_m=data.C_A-C_ADP_m;
C_ADP_ic=data.C_A-C_ATP_ic;
C_NAD=data.C_PN-C_NADH;
A_res=((data.R*data.T)/(data.F))*log(data.K_res.*sqrt(C_NADH./C_NAD));
delta_mu_h=-2.303*(data.R*data.T/data.F).*data.delta_pH+delta_Psi_m;

A_res_F=((data.R*data.T)/(data.F)).*log(data.K_resF.*sqrt(data.FADH2./data.FADH));

A_F1=((data.R*data.T)/data.F).*log((data.K_F1.*C_ATP_m)./(C_ADP_m.*data.Pi));

V_He= 6*data.rho_res*(data.r_a.*exp((A_res*data.F)./(data.R*data.T))-(data.r_a+data.r_b).*exp((12*data.F.*delta_mu_h)./(data.R*data.T)))./((1+data.r_1.*exp((A_res*data.F)./(data.R*data.T))).*exp((6*data.F*data.psi_B)./(data.R*data.T))+(data.r_2+data.r_3.* exp((A_res*data.F)./(data.R*data.T))).* exp((12*data.F*delta_mu_h)./(data.R*data.T)));

V_He_F=4*data.rho_resF*(data.r_a.*exp((A_res_F*data.F)./(data.R*data.T))-(data.r_a+data.r_b).*exp((12*data.F.*delta_mu_h)./(data.R*data.T)))./((1+data.r_1.* exp((A_res_F*data.F)./(data.R*data.T))).* exp((6*data.F*data.psi_B)./(data.R*data.T))+(data.r_2+data.r_3.* exp((A_res_F*data.F)./(data.R*data.T))).* exp((12*data.F*delta_mu_h)./(data.R*data.T)))                                                               ;

Juni=data.P_Ca*((data.Z_Ca*delta_Psi_m.*data.F)./(data.R*data.T)).*(data.alpha_m*Ca_m.*exp((-data.Z_Ca*delta_Psi_m.*data.F)./(data.R*data.T))-(data.alpha_e.*Ca_i))./(exp((-data.Z_Ca*delta_Psi_m.*data.F)./(data.R*data.T))-1);

Jnc=data.V_NC*(exp(0.5*delta_Psi_m.*data.F./(data.R*data.T)).*(data.Na_e^3.*Ca_m)./(data.K_Na^3.*data.K_Ca)-exp(-0.5*delta_Psi_m.*data.F/(data.R*data.T)).*(data.Na_m^3*Ca_i)./(data.K_Na^3*data.K_Ca))./(1+(data.Na_e^3/data.K_Na^3)+(Ca_m./data.K_Ca)+(data.Na_e^3.*Ca_m)./(data.K_Na^3*data.K_Ca)+(data.Na_m^3/data.K_Na^3)+(Ca_i/data.K_Ca)+(data.Na_m^3*Ca_i)./(data.K_Na^3*data.K_Ca));

V_Hu=-3*data.rho_F1*(10^2*data.p_a*(1+exp(A_F1./(data.R*data.T)))-(data.p_a+data.p_b)*exp((3*data.F*delta_mu_h/(data.R*data.T)))./((1+data.p_1* exp((A_F1*data.F)./(data.R*data.T))).*exp((3*data.F*data.psi_B)/(data.R*data.T))+( data.p_2+data.p_3*exp((A_F1.*data.F)./(data.R*data.T))).* exp((3*data.F.*delta_mu_h)./(data.R*data.T))));

V_H_Leak=data.g_H*delta_mu_h;

V_ANT=data.V_ant_max.*0.75.*(1-((0.25.*ATP_i.*0.45.*C_ADP_m)./(0.17.*ADPi.*0.025.*C_ATP_m))).*exp(-(data.F.*delta_Psi_m)./(data.R*data.T))./((1+((0.25.*ATP_i)./(0.225.*ADPi)).*exp((-data.h_ANT.*data.F.*delta_Psi_m)./(data.R*data.T))).*(1+((0.45.*C_ADP_m)./(0.025.*C_ATP_m ))));

ddelta_Psi_m=(V_He+V_He_F-V_Hu-V_ANT-V_H_Leak-Jnc-2*Juni)/data.C_mito;

%C_Cr_i=data.Cc-C_CrP_i;

%C_Cr_ic=data.Cc-C_CrP_ic;

%V_CK_mito=data.k_CK_mito*(ATP_i.*C_Cr_i-(ADPi.*C_CrP_i./data.K_EQ));

Force_ATP=1.02./(1+(data.K_M_ATP./ATP_i).*(1+(data.CATPi-ATP_i)./data.K_M_ADP));

V_AM=data.Max_ATP*data.F_f.*A.*Force_ATP;
ATP_XB = data.MaxATP.*Ff.*A.*Force; 

dATPi=V_ANT*(data.V_mito/data.V_myo)-V_AM-(Ip+ICap).*data.A_cap./(data.V_myo*data.F)-0.5*Iup;

%V_CK_cyto=data.k_CK_cyto*(C_ATP_ic.*C_Cr_ic-(C_ADP_ic.*C_CrP_ic./data.K_EQ));

%dC_ATP_ic=-V_CK_cyto-data.V_ATPase_cyto;
V_tr_CrP=data.k_tr_Cr*(C_CrP_i-C_CrP_ic);

%dC_CrP_i=V_CK_mito-V_tr_CrP;
%dC_CrP_ic=V_tr_CrP+V_CK_cyto;

%V_ATPase=-data.rho_F1*((10^2*data.p_a+data.p_c1.*exp((3*data.F*data.psi_B)/(data.R*data.T))).*exp(data.F*A_F1./(data.R*data.T))-(data.p_a.*exp((3*data.F.*delta_mu_h)./(data.R*data.T))+data.p_c2*exp(data.F.*A_F1./(data.R*data.T)).*exp((3*data.F.*delta_mu_h)./(data.R*data.T))))./((1+data.p_1* exp((A_F1.*data.F)./(data.R*data.T))).* exp((3*data.F*data.psi_B)./(data.R*data.T))+(data.p_2+data.p_3* exp((A_F1.*data.F)./(data.R*data.T))).* exp((3*data.F.*delta_mu_h)./(data.R*data.T)));
V_ATPase=-data.rho_F1*((10^2*data.p_a+data.p_c1.*exp((3*data.F*data.psi_B)./(data.R*data.T))).*exp(data.F*A_F1./(data.R*data.T))-(data.p_a.*exp((3*data.F*delta_mu_h)./(data.R*data.T))+data.p_c2*exp(data.F*A_F1./(data.R*data.T)).*exp((3*data.F*delta_mu_h)./(data.R*data.T))))./((1+data.p_1* exp((A_F1*data.F)./(data.R*data.T))).* exp((3*data.F*data.psi_B)./(data.R*data.T))+(data.p_2+data.p_3* exp((A_F1*data.F)./(data.R*data.T))).* exp((3*data.F*delta_mu_h)./(data.R*data.T))).*(1-exp(-(Ca_m./data.KCaATP)));

V_SL=data.k_f_SL*((C_SCoA.*C_ADP_m-(C_Suc.*C_ATP_m.*data.C_CoA)./(data.k_E_SL )));

dC_ADP_m=V_ANT-V_ATPase-V_SL;

V_O2=0.5*data.rho_res*((data.r_a+data.r_c1*exp((6*data.F*data.psi_B)/(data.R*data.T))).*exp((A_res*data.F)./(data.R*data.T))-data.r_a*exp((12*data.F*delta_mu_h)./(data.R*data.T))+data.r_c2*exp((A_res*data.F)./(data.R*data.T)).*exp((12*data.F*delta_mu_h)./(data.R*data.T)))./((1+data.r_1*exp((A_res*data.F)./(data.R*data.T))).*exp((6*data.F*data.psi_B)/(data.R*data.T))+(data.r_2+data.r_3*exp((A_res*data.F)./(data.R*data.T))).*exp((12*data.F*delta_mu_h)./(data.R*data.T)));

f_a_IDH=((1+(C_ADP_m./data.K_ADP_a)).*(1+(Ca_m./data.K_Ca_a))).^(-1);

f_i_IDH=1+(C_NADH./data.K_i_NADH);

V_IDH=data.k_cat_IDH*data.E_T_IDH*(1+(data.C_H/data.k_h_1)+(data.k_h_2/data.C_H)+f_i_IDH.*(data.K_M_NAD./C_NAD)+f_a_IDH.*(data.K_M_ISOC./C_ISOC).^2+f_a_IDH.*f_i_IDH.*(data.K_M_ISOC./C_ISOC).^2.*(data.K_M_NAD./C_NAD)).^-1;

f_a_KGDH=((1+(data.C_Mg/data.K_D_Mg)).*(1+(Ca_m./data.K_D_Ca))).^(-1);

V_KGDH=(data.k_cat_KGDH*data.E_T_KGDH)./(1+f_a_KGDH.*((data.K_M_aKG)./C_aKG).^(data.n_aKG)+f_a_KGDH.*((data.K_M_NAD_new)./C_NAD));

f_h_a=1./(1+data.C_H/data.k_h1+data.C_H^2./(data.k_h1*data.k_h2))+data.k_offset;

f_h_i=(1+data.k_h3./data.C_H+(data.k_h3*data.k_h4)./(data.C_H^2)).^-2;

V_MDH=(data.k_cat_MDH*data.E_T_MDH .*f_h_a.*f_h_i)./(1+(data.K_M_MAL./C_MAL).*(1+(C_OAA./data.K_i_OAA))+(data.K_M_NAD_mdh./C_NAD)+((data.K_M_MAL./C_MAL).*(1+(C_OAA./data.K_i_OAA)).*(data.K_M_NAD_mdh))./(C_NAD));

dC_NADH=-V_O2+V_IDH+V_KGDH+V_MDH;

C_CIT=data.C_K_int-(C_ISOC+C_aKG+C_SCoA+C_Suc+C_FUM+C_MAL+C_OAA);

V_ACO=data.k_f_ACO*(C_CIT-(C_ISOC./data.k_E_ACO));

dC_ISOC=V_ACO-V_IDH;

V_AAT=data.k_f_AAT.*C_OAA*data.C_GLU.*(data.k_ASP*data.K_E_AAT)./(data.k_ASP*data.K_E_AAT+C_aKG*data.k_f_AAT);

dC_aKG=V_IDH-V_KGDH+V_AAT;

dC_SCoA=V_KGDH-V_SL;

V_SDH=(data.k_cat_SDH*data.E_T_SDH)./(1+((data.K_M_SUC)./C_Suc).*(1+(C_OAA./data.K_i_sdh_OAA)).*((1+(C_FUM./(data.K_i_FUM)))));

dC_Suc=V_SL-V_SDH;

V_FH=data.k_f_FH*(C_FUM-C_MAL./(data.K_E_FH));

dC_FUM=V_SDH-V_FH;

dC_MAL=V_FH-V_MDH;

V_CS=data.k_cat_cs*data.E_T_cs*(1+(data.K_M_AcCoA/data.C_AcCoA)+(data.K_M_OAA./C_OAA)+(data.K_M_AcCoA/data.C_AcCoA)+(data.K_M_OAA./C_OAA)).^(-1);

dC_OAA=V_MDH-V_CS-V_AAT;

dC_FLV=V_SDH-V_O2;

dCa_m=data.Beta_Ca.*(Juni-Jnc);

impulseFactor =100;

if ((mod(t,1/f_stim) >= start_stim & mod(t,1/f_stim)< 0.001*impulseFactor + start_stim) & (t <= end_stim))
    I_stim=data.Is/impulseFactor;
else 
    I_stim=0;
end

dV = -(Itot + I_stim)./data.Cm;


%%

figure()
plot(t,V_O2)
figure()
plot(t,V_ANT)
%%
figure()
plot(t,ATP_i)


