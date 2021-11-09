function [I_Na, dm, dh1, dh2] = fast_na()

if(abs(V + 44.4) < 0.0001)
    Am = 460*12.673;
else
    Am = -460*(V + 44.4)/(exp(-(V + 44.4)/12.673) - 1);
end

Bm = 18400*exp(-(V + 44.4)/12.673);
dm = Am*(1 - m) - Bm*m;

Ah = 44.9*exp(-(V + 66.9)/5.57);
Bh = 1491.0/(1 + 323.3*exp(-(V + 94.6)/12.9));

th1 = 0.03/(1 + exp((V+40)/6)) + 0.00015; %0.00035 - Lindblad
th2 = 0.12/(1 + exp((V+60)/2)) + 0.00045;  %0.00295 - Lindblab

hm  = Ah/(Ah + Bh);
dh1 = (hm - h1)/th1;
dh2 = (hm - h2)/th2;

if(abs(V) > 0.0001)
    I_Na = 0.0014*m^3*(0.635*h1 + 0.365*h2)*140*V...
        *(F/RTONF)*(exp((V-E_Na)/RTONF) - 1)/(exp(V/RTONF) - 1); 
else
    I_Na = 0.0014*m^3*(0.635*h1 + 0.365*h2)*140*(F)*(exp((V-E_Na)/RTONF)-1);     
end

end