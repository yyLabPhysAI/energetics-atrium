%% Constant decleration and loading
constants
options = odeset('RelTol',1e-4,'AbsTol',1e-4);

%% Choice of stimulation parameters
start_stim = 0;     % [s]   When to start the current stimulation
end_stim   = 0;    % [s]   When to end it
f_stim     = 1;     % [Hz]  At what frequency do you want the current pulses to be?
TMAX       = 10;    % [s]   Until what time to calculate the simulation?
TMIN       = 0;     % [s]   Time to start plotting?

%% Run the model
tic;
[t,x] = ode23tb(@(t,x)model(t,x,data,f_stim,start_stim,end_stim), [0 TMAX] ,y0',options);
toc
ind_min = find(t>=TMIN);
t = t(ind_min:end);
V=x(ind_min:end,1); P_a=x(ind_min:end,2); P_i=x(ind_min:end,3); n=x(ind_min:end,4); r=x(ind_min:end,5); s_1=x(ind_min:end,6); s_2=x(ind_min:end,7); s_3=x(ind_min:end,8); m=x(ind_min:end,9); h_1=x(ind_min:end,10); h_2=x(ind_min:end,11); d_L=x(ind_min:end,12); f_L=x(ind_min:end,13); d_T=x(ind_min:end,14); f_T=x(ind_min:end,15); Na_i=x(ind_min:end,16); Ca_up=x(ind_min:end,17); Ca_rel=x(ind_min:end,18); Ca_i=x(ind_min:end,19); O_c=x(ind_min:end,20); O_TnCa=x(ind_min:end,21); O_TnMgCa=x(ind_min:end,22); O_TnMgMg=x(ind_min:end,23); O_Calse=x(ind_min:end,24); K_o=x(ind_min:end,25); K_i=x(ind_min:end,26); F_1=x(ind_min:end,27); F_2=x(ind_min:end,28); F_3=x(ind_min:end,29); SL=x(ind_min:end,31); A=x(ind_min:end,32); TT=x(ind_min:end,33); U=x(ind_min:end,34); Ve=x(ind_min:end,35); ATP_i=x(ind_min:end,36); Ca_m=x(ind_min:end,37); C_ATP_ic=x(ind_min:end,38); C_CrP_i=x(ind_min:end,39); C_CrP_ic=x(ind_min:end,40); C_ADP_m=x(ind_min:end,41); C_NADH=x(ind_min:end,42); delta_Psi_m=x(ind_min:end,43); C_ISOC=x(ind_min:end,44); C_aKG=x(ind_min:end,45); C_SCoA=x(ind_min:end,46); C_Suc=x(ind_min:end,47); C_FUM=x(ind_min:end,48); C_MAL=x(ind_min:end,49); C_OAA=x(ind_min:end,50); C_FLV=x(ind_min:end,51);  
