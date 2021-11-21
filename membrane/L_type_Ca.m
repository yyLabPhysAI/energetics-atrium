function [I_CaL, df_L, dd_L] = L_type_Ca(V, d_L, f_L)

g_CaL = 8.4;
E_CaL = 50;

if(abs(V + 45) < 0.0001)
    alpha_d_l = 16.72*2.5 - 50.0*(V + 10)/(exp(-(V + 10)/4.808) - 1);
elseif(abs(V + 10) < 0.0001)
    alpha_d_l = -16.72*(V + 45)/(exp(-(V + 45)/2.5) - 1) + 50*4.808;
else
    alpha_d_l = -16.72*(V + 45)/(exp(-(V + 45)/2.5) - 1) + 50.0*(V + 10)/(1 - exp(-(V + 10)/4.808));
end

if(abs(V + 5) < 0.0001)
    beta_d_l = 4.48*2.5;
else
    beta_d_l = - 4.48*(V + 5)/(1 - exp((V + 5)/2.5));
end

tau_d_l = 1/(alpha_d_l + beta_d_l);

d_lm = 1/(1 + exp(-((V + 10.95)/6.6)));
dd_L=(d_lm - d_L)/tau_d_l;

if(abs((V - 10) + 28) < 0.0001)
    alpha_f_l = 8.49*4;
else
    alpha_f_l = -8.49*(V + 18)/(1 - exp((V + 18)/4));
end

beta_f_l = 67.922/(1 + exp(-((V + 18)/4)));

f_lm = alpha_f_l/(alpha_f_l + beta_f_l);
t_f_l = 1/(alpha_f_l + beta_f_l);
df_L = (f_lm - f_L)/t_f_l;

I_CaL = g_CaL*(d_L*f_L + 1/(1 + exp(-(V - 23)/12)))*(V - E_CaL);

end

