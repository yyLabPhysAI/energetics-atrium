function [dPa, dPi] = fast_delayed_rectifier_k(V, Pa, Pi)

Apa = 9.0*exp(V/25.371);
Bpa = 1.3*exp(-V/13.026);
pam = 1/(1 + exp(-(V + 5.1)/7.4));
tpa = 1/(Apa+Bpa);
dPa = (pam - Pa)/tpa;

Api = 100*exp(-V/54.645);
Bpi = 656*exp(V/106.157);
pim = 1/(1+exp((V+47.3921)/18.6603));
tpi = 1/(Api+Bpi);
dPi = (pim-Pi)/tpi;

end