function I_k1 = inward_rectifier_k(V, E_k, K_o, data)

E_0 = V - E_k + 3.6;
g_K1 = 5.08;
k_mK1 = 0.59;
I_k1 = 2.5*g_K1 * (K_o/(K_o + k_mK1))^3 * (V - E_k)/(1 + exp(1.393*E_0/data.RTONF));

end
