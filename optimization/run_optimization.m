constants
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
%% Choice of stimulation parameters
start_stim = 10;     % [s]   When to start the current stimulation
end_stim   = 10;    % [s]   When to end it
f_stim     = 1;     % [Hz]  At what frequency do you want the current pulses to be?
TMAX       = 10;    % [s]   Until what time to calculate the simulation?
TMIN       = 0;     % [s]   Time to start plotting?

%% Run the model
l = Inf;
loss_list = [];
eta = 0.1;
for i = 1:1000

    [new_data, new_y0] = step(data, y0, eta);
    
    [t,x] = ode23tb(@(t,x)model(t, x, new_data, f_stim, start_stim, end_stim), [0 TMAX] ,new_y0',options);
    
    current_l = loss(x);
    if (0 < current_l) && (current_l < l) 
        data = new_data;
        y0 = new_y0;
        l = current_l;
        loss_list = [loss_list, [i;l]];
        eta = eta*0.99;
        disp([i, current_l])
    elseif ~mod(i, 100)
        disp(i)
    end
end

figure()
plot(loss_list(1, :), loss_list(2, :))
%%
V=x(:,1); P_a=x(:,2); P_i=x(:,3); n=x(:,4); r=x(:,5); s_1=x(:,6); s_2=x(:,7); s_3=x(:,8); m=x(:,9); h_1=x(:,10); h_2=x(:,11); d_L=x(:,12); f_L=x(:,13); d_T=x(:,14); f_T=x(:,15); Na_i=x(:,16); Ca_up=x(:,17); Ca_rel=x(:,18); Ca_i=x(:,19); O_c=x(:,20); O_TnCa=x(:,21); O_TnMgCa=x(:,22); O_TnMgMg=x(:,23); O_Calse=x(:,24); K_o=x(:,25); K_i=x(:,26); F_1=x(:,27); F_2=x(:,28); F_3=x(:,29); SL=x(:,31); A=x(:,32); TT=x(:,33); U=x(:,34); Ve=x(:,35); ATP_i=x(:,36); Ca_m=x(:,37); C_ATP_ic=x(:,38); C_CrP_i=x(:,39); C_CrP_ic=x(:,40); C_ADP_m=x(:,41); C_NADH=x(:,42); delta_Psi_m=x(:,43); C_ISOC=x(:,44); C_aKG=x(:,45); C_SCoA=x(:,46); C_Suc=x(:,47); C_FUM=x(:,48); C_MAL=x(:,49); C_OAA=x(:,50); C_FLV=x(:,51);  

figure()
subplot(4, 4, 1); nice_plot(C_FLV, t);
subplot(4, 4, 2); nice_plot(C_OAA, t);
subplot(4, 4, 3); nice_plot(ATP_i, t);
subplot(4, 4, 4); nice_plot(Ca_m, t);
subplot(4, 4, 5); nice_plot(C_ADP_m, t);
subplot(4, 4, 6); nice_plot(C_NADH, t);
subplot(4, 4, 7); nice_plot(delta_Psi_m, t);
subplot(4, 4, 8); nice_plot(C_ISOC, t);
subplot(4, 4, 9); nice_plot(C_aKG, t);
subplot(4, 4, 10); nice_plot(C_SCoA, t);
subplot(4, 4, 11); nice_plot(C_Suc, t);
subplot(4, 4, 12); nice_plot(C_FUM, t);
subplot(4, 4, 13); nice_plot(C_MAL, t);
subplot(4, 4, 14); nice_plot(C_FLV, t);
