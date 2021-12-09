function nice_plot(var, t)

n = inputname(1);
if contains(n, '_')
    n = [replace(n, "_", "_{"), "}"];
end

plot(t, var, "LineWidth", 2)
title(n)
xlabel('Time [sec]')
ylabel(n)
xlim([t(1), t(end)]);
set(gca,'FontSize',10)
set(gcf,'color','white')
set(gca, 'FontName', 'Times New Roman');
set(gca,'linewidth',2)
set(gca,'FontWeight','bold')

end
