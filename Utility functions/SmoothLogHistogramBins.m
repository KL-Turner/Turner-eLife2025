function [curve] = SmoothLogHistogramBins(histData,edges)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%----------------------------------------------------------------------------------------------------------
curveFig = figure;
h = histogram(log10(histData),edges,'Normalization','probability');
% set(gca,'xscale','log')
heights = h.Values;
centers = h.BinEdges + h.BinWidth/2;
centers = centers(1:end-1);
hold on
ax = gca;
ax.XTickLabel = [];
n = length(centers);
w = centers(2) - centers(1);
t = linspace(centers(1) - w/2,centers(end) + w/2,n + 1);
p = fix(n/2);
fill(t([p,p,p + 1,p + 1]),[0 heights([p,p]),0],'w')
plot(centers([p,p]),[0,heights(p)],'r:')
h = text(centers(p) - .2,heights(p)/2,'   h');
dep = -70;
tL = text(t(p),dep,'L');
tR = text(t(p + 1),dep,'R');
hold off
dt = diff(t);
Fvals = cumsum([0,heights.*dt]);
F = spline(t,[0,Fvals,0]);
curve = fnder(F);
h.String = 'h(i)';
tL.String = 't(i)';
tR.String = 't(i+1)';
hold on
fnplt(curve,'r',2)
hold off
ylims = ylim;
ylim([0,ylims(2)]);
close(curveFig)