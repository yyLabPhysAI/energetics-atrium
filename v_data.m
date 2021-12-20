load('v1.mat')
load('v2.mat')
load('v3.mat')
load('v4.mat')

figure()

tm = (0:1/10000: (length(v1)-1)/10000)';
plot(tm,v1)
hold on
tm = (0:1/10000: (length(v2)-1)/10000)';
plot(tm,v2)
hold on
tm = (0:1/10000: (length(v3)-1)/10000)';
plot(tm,v3)
hold on
tm = (0:1/10000: (length(v4)-1)/10000)';
plot(tm,v4)
hold on

legend('V_1', 'V_2', 'V_3', 'V_4')