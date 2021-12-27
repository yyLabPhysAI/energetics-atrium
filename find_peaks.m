function boolim = find_peaks(signal_vec)
    
    means = zeros(length(signal_vec),1);
    for ind = -100:100
        means = means + circshift(signal_vec,ind);
    end
    means = means/201;
    
    SDs = zeros(length(signal_vec),1);
    for ind = -100:100
        SDs = SDs + (circshift(signal_vec,ind)-means).^2;
    end
    SDs = SDs/201;
    SDs = SDs.^0.5;
    
    %dis_vec = abs(signal_vec-means);
    dis_vec = (signal_vec-means);
    
    dis_env = [];
    for ind = -20:20
        dis_env = [dis_env circshift(dis_vec,ind)];
    end
    
    boolim = dis_vec' == max(dis_env') | dis_vec' == min(dis_env');
    %boolim = boolim & (dis_vec>3*SDs)';
    boolim = boolim & (dis_vec>1*SDs)';
end
