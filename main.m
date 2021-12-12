clear 
clc
close all

%% Constant decleration and loading
constants
options = odeset('RelTol',1e-4,'AbsTol',1e-4);

%% Choice of stimulation parameters
start_stim = 100;     % [s]   When to start the current stimulation
end_stim   = 100;    % [s]   When to end it
f_stim     = 1;     % [Hz]  At what frequency do you want the current pulses to be?
TMAX       = 100;    % [s]   Until what time to calculate the simulation?
TMIN       = 0;     % [s]   Time to start plotting?

%% Run the model
tic;
[t,x] = ode23tb(@(t,x)model(t,x,data,f_stim,start_stim,end_stim), [0 TMAX] ,y0',options);
toc
ind_min = find(t>=TMIN);
t = t(ind_min:end);
V=x(ind_min:end,1); P_a=x(ind_min:end,2); P_i=x(ind_min:end,3); n=x(ind_min:end,4); r=x(ind_min:end,5); s_1=x(ind_min:end,6); s_2=x(ind_min:end,7); s_3=x(ind_min:end,8); m=x(ind_min:end,9); h_1=x(ind_min:end,10); h_2=x(ind_min:end,11); d_L=x(ind_min:end,12); f_L=x(ind_min:end,13); d_T=x(ind_min:end,14); f_T=x(ind_min:end,15); Na_i=x(ind_min:end,16); Ca_up=x(ind_min:end,17); Ca_rel=x(ind_min:end,18); Ca_i=x(ind_min:end,19); O_c=x(ind_min:end,20); O_TnCa=x(ind_min:end,21); O_TnMgCa=x(ind_min:end,22); O_TnMgMg=x(ind_min:end,23); O_Calse=x(ind_min:end,24); K_o=x(ind_min:end,25); K_i=x(ind_min:end,26); F_1=x(ind_min:end,27); F_2=x(ind_min:end,28); F_3=x(ind_min:end,29); SL=x(ind_min:end,31); A=x(ind_min:end,32); TT=x(ind_min:end,33); U=x(ind_min:end,34); Ve=x(ind_min:end,35); ATP_i=x(ind_min:end,36); Ca_m=x(ind_min:end,37); C_ATP_ic=x(ind_min:end,38); C_CrP_i=x(ind_min:end,39); C_CrP_ic=x(ind_min:end,40); C_ADP_m=x(ind_min:end,41); C_NADH=x(ind_min:end,42); delta_Psi_m=x(ind_min:end,43); C_ISOC=x(ind_min:end,44); C_aKG=x(ind_min:end,45); C_SCoA=x(ind_min:end,46); C_Suc=x(ind_min:end,47); C_FUM=x(ind_min:end,48); C_MAL=x(ind_min:end,49); C_OAA=x(ind_min:end,50); C_FLV=x(ind_min:end,51);  

%% Action potential, membrane currents and intracellular Ca2+

% Nernst potentials
Na_o = data.Na_o; Ca_o = data.Ca_o;

E_k   = nernst(K_i, K_o, 1, data);
E_Na  = nernst(Na_i, Na_o, 1, data);
E_Ca  = nernst(Ca_i, Ca_o, 2, data);

% Fast delayed rectifier K+ current
[I_Kr, ~, ~] = fast_delayed_rectifier_k(V, P_a, P_i, E_k);

% Slow delayed rectifier K+ current
[I_Ks, ~] = slow_delayed_rectifier_k(V, n, E_k);

% Inward rectifier K+ channel
I_k1 = inward_rectifier_k(V, E_k, K_o, data);

% Transient outward K+ channel
[I_Kto, ~, ~, ~, ~] = transient_outward_k(r, s_1, s_2, s_3, V, E_k);

% Na+-K+ ATPase
I_NaK = NaK_ATPase(K_o, Na_i, V);

% Na+-Ca(2+) exchanger (NCX)
I_NaCa = NCX(Na_i, Ca_i, V, data);

% Fast Na+ channels
[I_Na, ~, ~, ~] = fast_na(V, m, h_1, h_2, E_Na, data);

% Membranal Ca(2+) pump
I_Cap_max = 9.509;
I_Cap = I_Cap_max.*(Ca_i./(Ca_i + 0.0002));

% L-type Ca(2+) channels
[I_CaL, ~, ~] = L_type_Ca(V, d_L, f_L);

% T-Type calcium channels
[~ , ~, I_CaT] = T_Type_calcium(V, d_T, f_T);

% Background currents
g_NaB = 0.03;
g_CaB = 0.03;

I_NaB = g_NaB.*(V - E_Na);
I_CaB = g_CaB.*(V - E_Ca);

% Calculate the total current
I_tot = I_Kr + I_Ks + I_k1 + I_Kto + I_NaK + I_NaCa + I_Na + I_NaB + I_CaL + I_CaT + I_Cap + I_CaB;


%% Mitochondrial energy metabolism, Ca2+ dynamics and oxygen consumption

C_PN = data.C_PN;
C_NAD = C_PN - C_NADH;
C_A = data.C_A;
ADP_i = C_A - ATP_i;
ATP_m = C_A - C_ADP_m;

[~, ~, ~, ~, ~, ~, ~, V_SL, V_IDH, V_KGDH, V_MDH, V_SDH] = TCA_cycle(...
    C_ISOC, C_aKG, C_SCoA, C_Suc, C_FUM, C_MAL, ...
    C_OAA, C_NAD, C_ADP_m, Ca_m, C_NADH, ATP_m, data);

[~, ~, ~, ~, ~, ~, ~, V_ANT, V_O2]...
    = oxidative_phosphorylation(V_SL, V_IDH, V_KGDH, V_MDH, V_SDH, ...
    delta_Psi_m, C_NADH, C_NAD, ATP_m, C_ADP_m, Ca_m, ATP_i, ...
    ADP_i, data);

%% Membrane potential
figure()

nice_plot(V, t)

%% Membranal currents
figure()

subplot(3, 4, 1); nice_plot(I_Kr, t) 
subplot(3, 4, 2); nice_plot(I_Ks, t) 
subplot(3, 4, 3); nice_plot(I_k1, t) 
subplot(3, 4, 4); nice_plot(I_Kto, t) 
subplot(3, 4, 5); nice_plot(I_NaK, t) 
subplot(3, 4, 6); nice_plot(I_NaCa, t) 
subplot(3, 4, 7); nice_plot(I_Na, t) 
subplot(3, 4, 8); nice_plot(I_NaB, t) 
subplot(3, 4, 9); nice_plot(I_CaL, t) 
subplot(3, 4, 10); nice_plot(I_CaT, t) 
subplot(3, 4, 11); nice_plot(I_Cap, t) 
subplot(3, 4, 12); nice_plot(I_CaB, t)

%% Calcium 
figure()

subplot(2, 2, 1); nice_plot(Ca_up, t)
subplot(2, 2, 2); nice_plot(Ca_rel, t)
subplot(2, 2, 3); nice_plot(Ca_i, t)
subplot(2, 2, 4); nice_plot(Ca_m, t)

%% ATP
figure()

subplot(1, 2, 1); nice_plot(ATP_i, t)
subplot(1, 2, 2); nice_plot(ATP_m, t)

%% V_O2
figure()

nice_plot(V_O2, t)

%% FLV
figure()

nice_plot(C_FLV, t)
