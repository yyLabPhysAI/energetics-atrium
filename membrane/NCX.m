function I_NaCa = NCX(Na_i, Ca_i, V, data)

Ca_o = data.Ca_o;
Na_o = data.Na_o;
gamma_NaCa = 0.450;
d_NaCa = 0.0003*20;

I_NaCa = 0.02.*(Na_i.^3.*Ca_o.*exp(gamma_NaCa.*V./data.RTONF) - Na_o.^3.*Ca_i.*exp(V.*(gamma_NaCa - 1)./data.RTONF))...
                ./(1 + d_NaCa.*(Ca_i.*Na_o.^3 + Ca_o.*Na_i.^3));

end