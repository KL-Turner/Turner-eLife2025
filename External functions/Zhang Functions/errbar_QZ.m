function errbar_QZ(X,Y,E,color)
errbar(X,Y,E,'linestyle','-','color',color,'LineWidth',1);
%% Re-plot the horizontal lines
% Width of the top and bottom lines of errorbar
xlength = 0.04;
% Make horizontal lines with 'line'
for k = 1:length(X)
    x_h = [X(k) - xlength, X(k) + xlength];
    y_h = [Y(k) + E(k), Y(k) + E(k)];
    line(x_h, y_h,'color',color, 'LineWidth',1);
    y_b = [Y(k) - E(k), Y(k) - E(k)];
    line(x_h, y_b,'color',color, 'LineWidth',1);
end
end