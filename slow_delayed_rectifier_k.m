function [I_Ks, dn] = slow_delayed_rectifier_k(V, n, E_k)

g_Ks = 2.5; % ????????? 
I_Ks = g_Ks*n*(V - E_k);

alpha_n = 1.66*exp(V/69.452);
beta_n = 0.3*exp(-V/21.826);
n_m = 1/(1 + exp(-(V - 0.9)/13.8));
tau_n = 1/(alpha_n + beta_n) + 0.06;
dn = (n_m - n)/tau_n;

end
