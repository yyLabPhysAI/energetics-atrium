sb
v_data
plot(t, V)
ylim([-70, 35])
xlim([1.01, 1.135])



% Action potential, membrane currents and intracellular Ca2+

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


% figure(9)
% plot(t, V)
% hold on
% plot(t, I_Kr)
% plot(t, I_Ks)
% plot(t, I_k1)
% plot(t, I_Kto)
% plot(t, I_NaK)
% plot(t, I_NaCa)
% plot(t, I_Na)
% plot(t, I_NaB)
% plot(t, I_CaL)
% plot(t, I_CaT)
% plot(t, I_Cap)
% plot(t, I_CaB)
% plot(t, I_tot)
% xlim([0, 0.1])
% 
% legend(["V";
% "I_Kr";
% "I_Ks";
% "I_k1";
% "I_Kto";
% "I_NaK";
% "I_NaCa";
% "I_Na";
% "I_NaB";
% "I_CaL";
% "I_CaT";
% "I_Cap";
% "I_CaB";
% "I_tot"])
% 

%%
%close all

[dO_c, dO_TnCa, dO_TnMgCa, dO_TnMgMg, dO_Calse, phi_ca_i, ...
    dCa_up, dCa_rel, dF_1, dF_2, dF_3, I_up, I_rel] ...
    = ...
    SR_calcium_handling(ATP_i, Ca_i, Ca_up, Ca_rel, F_1, F_2, ...
    F_3, O_c, O_TnCa, O_TnMgCa, O_TnMgMg, O_Calse, V, data);

figure(10)
hold off
plot(t, 2.*I_NaCa)
hold on
plot(t, -I_CaL)
hold on
plot(t, -I_CaT)
hold on
plot(t, -I_Cap)
hold on
plot(t, -I_CaB)
hold on
plot(t, -I_up)
hold on
plot(t, I_rel)

legend("I_NaCa","I_CaL","I_CaT","I_Cap","I_CaB","I_up","I_rel")

%%

[...
dSL, dA, dTT, dU, dV_e, Force...
] = force_generation(Ve, SL, TT, A, U, data, Ca_i, ATP_i, ADP_i);


Ca_ratio = sum(I_rel(500:end)/sum(I_CaL(500:end)))

%% 

[V_AM, ATP_XB] = force_energy_consumption(Force, ATP_i, A, data);

SR_consumption = sum(0.5.*I_up/(data.F*data.V_myo));
pump_consumption = sum((I_NaK + I_Cap).*data.A_cap./(data.V_myo.*data.F));
sarcomere_consumtion = sum(V_AM);

con_vec = [SR_consumption, pump_consumption, sarcomere_consumtion];


figure(11)
labels = {'SR','Membrane','Sarcomere'};
pie(con_vec/sum(con_vec),[1,1,1], labels)

%%

figure(12); plot(t, Force)
legend("force")