function l = loss(x)

V=x(:,1); P_a=x(:,2); P_i=x(:,3); n=x(:,4); r=x(:,5); s_1=x(:,6); s_2=x(:,7); s_3=x(:,8); m=x(:,9); h_1=x(:,10); h_2=x(:,11); d_L=x(:,12); f_L=x(:,13); d_T=x(:,14); f_T=x(:,15); Na_i=x(:,16); Ca_up=x(:,17); Ca_rel=x(:,18); Ca_i=x(:,19); O_c=x(:,20); O_TnCa=x(:,21); O_TnMgCa=x(:,22); O_TnMgMg=x(:,23); O_Calse=x(:,24); K_o=x(:,25); K_i=x(:,26); F_1=x(:,27); F_2=x(:,28); F_3=x(:,29); SL=x(:,31); A=x(:,32); TT=x(:,33); U=x(:,34); Ve=x(:,35); ATP_i=x(:,36); Ca_m=x(:,37); C_ATP_ic=x(:,38); C_CrP_i=x(:,39); C_CrP_ic=x(:,40); C_ADP_m=x(:,41); C_NADH=x(:,42); delta_Psi_m=x(:,43); C_ISOC=x(:,44); C_aKG=x(:,45); C_SCoA=x(:,46); C_Suc=x(:,47); C_FUM=x(:,48); C_MAL=x(:,49); C_OAA=x(:,50); C_FLV=x(:,51);  

l =   std(C_FLV)				./(abs(mean(C_FLV)) + 1e-10)...
 	+ std(C_OAA)				./(abs(mean(C_OAA)) + 1e-10)...
 	+ std(ATP_i)				./(abs(mean(ATP_i)) + 1e-10)...
 	+ std(Ca_m)				./(abs(mean(Ca_m)) + 1e-10)...
 	+ std(C_ADP_m)			./(abs(mean(C_ADP_m)) + 1e-10)...
 	+ std(C_NADH)			./(abs(mean(C_NADH)) + 1e-10)...
 	+ std(delta_Psi_m)		./(abs(mean(delta_Psi_m)) + 1e-10)...
 	+ std(C_ISOC)			./(abs(mean(C_ISOC)) + 1e-10)...
 	+ std(C_aKG)				./(abs(mean(C_aKG)) + 1e-10)...
 	+ std(C_SCoA)			./(abs(mean(C_SCoA)) + 1e-10)...
 	+ std(C_Suc)				./(abs(mean(C_Suc)) + 1e-10)...
 	+ std(C_FUM)				./(abs(mean(C_FUM)) + 1e-10)...
 	+ std(C_MAL)				./(abs(mean(C_MAL)) + 1e-10)...
 	+ std(C_OAA)				./(abs(mean(C_OAA)) + 1e-10)...
  	+ std(Ca_rel)			./(abs(mean(Ca_rel)) + 1e-10)...
    + std(Ca_up)				./(abs(mean(Ca_up)) + 1e-10)...
    + std(C_FLV)				./(abs(mean(C_FLV)) + 1e-10);

l = abs(l);

end
