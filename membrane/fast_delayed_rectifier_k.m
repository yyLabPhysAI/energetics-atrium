function [I_Kr, dP_a, dP_i] = fast_delayed_rectifier_k(V, P_a, P_i, E_k)

g_Kr = 3.5; % Condactance of the fast delayed rectifier k+ [nS]
I_Kr = g_Kr.*P_a.*P_i.*(V - E_k);

alpha_p_a = 9.0.*exp(V./25.371);
beta_p_a = 1.3.*exp(-V./13.026);
p_am = 1./(1 + exp(-(V + 5.1)./7.4));
tau_p_a = 1./(alpha_p_a + beta_p_a);
dP_a = (p_am - P_a)./tau_p_a;

alpha_p_i = 100.*exp(-V./54.645);
beta_p_i = 656.*exp(V./106.157);
p_im = 1./(1 + exp((V + 47.3921)./18.6603));
tau_p_i = 1./(alpha_p_i + beta_p_i);
dP_i = (p_im - P_i)./tau_p_i;

end