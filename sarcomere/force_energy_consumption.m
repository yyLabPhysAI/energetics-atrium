function [V_AM, ATP_XB] = force_energy_consumption(Force, ATP_i, A, data)

Ff=data.F_f;
MaxATP=data.MaxATP;
F_f=data.F_f;
K_M_ATP=data.K_M_ATP;
Max_ATP=data.Max_ATP;
CATPi=data.CATPi;
K_M_ADP = data.K_M_ADP;

Force_ATP = 1.02./(1 + (K_M_ATP./ATP_i).*(1 + (CATPi - ATP_i)./K_M_ADP));
V_AM = Max_ATP.*F_f.*A.*Force_ATP;
ATP_XB = MaxATP*Ff*A.*Force;

end