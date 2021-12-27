function [dd_T , df_T, I_CaT] = T_Type_calcium(V, d_T, f_T)

g_CaT = 6*10;
E_CaT = 38;

alpha_d_t = 674.173.*exp((V + 23.3)./30);
beta_d_t = 674.173.*exp(-(V + 23.3)./30);
tau_d_t = 1./(alpha_d_t + beta_d_t);
d_tm = 1./(1 + exp(-(V + 23)./6.1));
dd_T = (d_tm - d_T)./tau_d_t;

alpha_f_t = 9.637.*exp(-(V + 75)./83.3);
beta_f_t = 9.637.*exp((V + 75)./15.38);
tau_f_t = 1./(alpha_f_t + beta_f_t);
f_tm = alpha_f_t./(alpha_f_t + beta_f_t);
df_T = (f_tm - f_T)./tau_f_t;

I_CaT = g_CaT.*d_T.*f_T.*(V - E_CaT);

end