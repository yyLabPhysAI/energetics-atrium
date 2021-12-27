function [I_Kto, dr, ds_1, ds_2, ds_3] = transient_outward_k(r, s_1, s_2, s_3, V, E_k)

g_Kto = 50.02*10;
I_Kto = 0.35.*g_Kto.*r.*(0.59.*s_1.^3 + 0.41.*s_2.^3).*(0.6.*s_3.^6 + 0.4).*(V - E_k);

g_sus = 2.4; % Conductance of the sustained current [nS]
E_sus = 70; % Reversal potential of the sustained current [mV]
Isus = g_sus.*(V + E_sus);
I_Kto  = I_Kto + Isus;

alpha_r = 386.6.*exp(V./12);
beta_r = 8.011.*exp(-V./7.2);
r_m = 1./(1 + exp(-(V + 15)./5.633));
tau_r = 1./(alpha_r + beta_r) + 0.0004;
dr = (r_m - r)./tau_r;

s_1m = 1./(1 + exp((V + 28.29)./7.06));
tau_s_1 = 0.5466./(1 + exp((V + 32.8)./0.1)) + 0.0204;
ds_1 = (s_1m - s_1)./tau_s_1;

s_2m = 1./(1 + exp((V + 28.29)./7.06));
tau_s_2 = 5.75./(1 + exp((V + 32.8)./0.1)) + 0.45./(1 + exp(-(V - 13.54)./13.97));
ds_2 = (s_2m - s_2)./tau_s_2;

s_3m = ((1./(1 + exp((V + 50.67)./27.38))) + 0.666)./1.666;
tau_s_3 = 7.5./(1 + exp((V + 23.0)./0.5)) + 0.5;
ds_3 = (s_3m - s_3)./tau_s_3;

end
