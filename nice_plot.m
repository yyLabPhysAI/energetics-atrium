function nice_plot(var, t)
var=real(var);
n = inputname(1);
closes = count(n, "_");
if contains(n, '_')
    n = [replace(n, "_", "_{"), repelem("}", closes)];
end

plot(t, var, "LineWidth", 2)
title(n)
xlabel('Time [sec]')
ylabel(n)
xlim([t(1), t(end)]);
ylim([1.5*min([0, min(var)]), 1.5*max([0, max(var)])]);
set(gca,'FontSize',10)
set(gcf,'color','white')
set(gca, 'FontName', 'Times New Roman');
set(gca,'linewidth',2)
set(gca,'FontWeight','bold')
box off
end
