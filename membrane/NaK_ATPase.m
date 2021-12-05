function I_NaK = NaK_ATPase(K_o, Na_i, V)

I_NaKmax = 64.41;
k_mK = 1;
k_mNa = 11.^1.5;

I_NaK = I_NaKmax .*(K_o./(K_o + k_mK)).*(Na_i.^1.5./(Na_i.^1.5 + k_mNa)).*(1.6./(1.5 + exp(-(V + 60)./40)));

end