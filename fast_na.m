function [I_Na, dm, dh_1, dh_2] = fast_na(V, m, h_1, h_2, E_Na, data)

RTONF = data.RTONF;
F = data.F;
g_Na = 0.0014;
Na_o = data.Na_o;

h = 0.635*h_1 + 0.365*h_2; % ??????????????? diff from doc!!!!!

if(abs(V + 44.4) < 0.0001)
    alpha_m = 460*12.673; % The limit calculated L'Hôpital's rule
else
    alpha_m = -460*(V + 44.4)/(exp(-(V + 44.4)/12.673) - 1);
end

beta_m = 18400*exp(-(V + 44.4)/12.673);
dm = alpha_m*(1 - m) - beta_m*m;

alpha_h = 44.9*exp(-(V + 66.9)/5.57);
beta_h = 1491.0/(1 + 323.3*exp(-(V + 94.6)/12.9));

tau_h_1 = 0.03/(1 + exp((V + 40)/6)) + 0.00015; %0.00035 - Lindblad ????? so which one do you want??
tau_h_2 = 0.12/(1 + exp((V + 60)/2)) + 0.00045;  %0.00295 - Lindblab ????? so which one do you want??

h_m  = alpha_h/(alpha_h + beta_h);
dh_1 = (h_m - h_1)/tau_h_1;
dh_2 = (h_m - h_2)/tau_h_2;

if(abs(V) > 0.0001)
    I_Na = g_Na*m^3*h*Na_o*V...
        *(F/RTONF)*(exp((V - E_Na)/RTONF) - 1)/(exp(V/RTONF) - 1);
else
    I_Na = g_Na*m^3*h*Na_o*(F)*(exp((V-E_Na)/RTONF)-1);
end

end