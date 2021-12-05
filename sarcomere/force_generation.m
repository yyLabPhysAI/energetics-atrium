function [...
    dSL, dA, dTT, dU, dV_e, Force...
    ] = force_generation(V_e, SL, TT, A, U, data, Ca_i, ATP_i, ADP_i)

F_kl = data.F_kl;
F_f = data.F_f;
F_go = data.F_go;
F_gl = data.F_gl;
SL_0 = data.SL_0;
N_c = data.N_c;
F_k0 = data.F_k0;
F_k1 = data.F_k1;
FN = data.FN;
F_k05 = data.F_k05;
k_XBATP = data.k_XBATP;
k_XBADP = data.k_XBADP;
F_XB = data.F_XB;


%  Cell Contraction
dV_e = 0; % Assume constant contraction velocity
dSL = -V_e;


N_XB = (1e-6).*(SL - SL_0).*N_c.*(TT + U).*1000./2;
K_Ca = F_k0 + F_k1.*(N_XB.^FN)./(F_k05.^FN + N_XB.^FN);
K_minus1 = F_kl./K_Ca;
dA = F_kl.*Ca_i.*(1 - A - TT - U) - A.*(F_f + K_minus1) + TT.*(F_go + F_gl.*V_e);
dTT = F_f.*A - TT.*(F_go + F_gl.*V_e + K_minus1) + F_kl.*Ca_i.*U;
dU = K_minus1.*TT - (F_go + F_gl.*V_e + F_kl.*Ca_i).*U;


N_XB_ATP = 1.02./(1 + (k_XBATP./ATP_i).*(1 + ADP_i./k_XBADP));
Force = F_XB.*N_XB_ATP.*((SL - SL_0)./2).*(TT + U).*N_c;

end
