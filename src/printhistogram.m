% Print histograms for the rapidity gap simulation
%
% mikael.mieskolainen@cern.ch, 2019

function printhistogram(histo, param, cuts)

k = 1;

% Create figure
figure1 = figure('units','normalized','outerposition',[0 0 1 1]);

% ------------------------------------------------------------------------
% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.325932504440498 0.388009049773757 0.368561278863232 0.469457013574661]);
hold(axes1,'on');

% Plot ND,SD,DD
for kk = 1:size(histo.hMC,1)
    stephistedge(histo.binedges, histo.hMC(kk,:), 'linewidth', 1.5); hold on;    
end

% Total sum
hTOTAL = sum(histo.hMC,1);
stephistedge(histo.binedges, histo.hTOTAL, 'k-', 'linewidth', 1.5);

% LHC DATA [x,y,neg,pos]
errorbar(histo.xPOINTS, histo.hDATA, histo.hDATA.*(histo.hDATA_D/100), histo.hDATA.*(histo.hDATA_U/100), ...
    's', 'Color', zeros(3,1), 'MarkerFaceColor', ones(3,1)*0.7, 'MarkerSize', 4, 'CapSize', 0, 'LineWidth', 1.25);

% ------------------------------------------------------------------------

% label
str = 'F';
ylabel(sprintf('$d \\sigma / d \\Delta \\eta^%s$ (mb)', str),'interpreter','latex');

axis square;
set(gca,'yscale','log');
set(gca, 'XTick', 0:1:max(histo.binedges));
axis([0 8 1e-2 1e4]);

%
l = legend(sprintf('ND: %0.1f/%0.1f mb', 0, (1 - sum(param.diff)) * param.sigmainel), ...
           sprintf('SD: %0.1f/%0.1f mb', 0, param.diff(1) * param.sigmainel), ...
           sprintf('DD: %0.1f/%0.1f mb', 0, param.diff(2) * param.sigmainel), ...
           sprintf('Total: %0.1f/%0.1f mb', 0, param.sigmainel), ...
           'ATLAS 7 TeV');

set(l,'interpreter','latex');
legend('boxoff');
%}

if (k == 1)
bstr = 'I';
end
if (k == 2)
bstr = 'II';
end

title(sprintf('Type %s Boundary Conditions', bstr),'interpreter','latex');
text(0.4, 2*10^3, sprintf('$p_T > %0.1f$ GeV',  cuts.ptcut), 'interpreter','latex', 'fontsize', 12);
text(0.4, 1*10^3, sprintf('$| \\eta| <  %0.1f$', cuts.etacut), 'interpreter','latex', 'fontsize', 12);

axis(axes1,'square');
set(axes1,'XTickLabel', []);
set(gca,'box','on');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Interpreter','latex');

% Remove the first y-axis tick due to overlap with underneath plot
%ax = get(gca,'Parent');
yTick = get(gca,'YTickLabel'); yTick{1} = '';
set(gca,'YTickLabel',yTick);

% ------------------------------------------------------------------------
% Create lower axes
axes2 = axes('Parent',figure1,...
    'Position',[0.395 0.256216923529732 0.23 0.126625586073341]);
hold(axes2,'on');

% Plot ratio plot with data uncertainty (fill bars)
plot([0 4 8], ones(3,1), 'k-'); hold on; % Horizontal line
stephistedge(histo.binedges, histo.hTOTAL(:) ./ histo.hDATA(:), 'k', 'linewidth', 1);
stepfilledge(histo.binedges, 1+histo.hDATA_U/100, 1-histo.hDATA_D/100, [1 0 0]*0.5, [1 0 0]*0.5, 0.1);

% Axes properties
set(axes2,'XTick',[0 1 2 3 4 5 6 7 8],'YMinorTick','on','YScale','linear');
axis([0 8 0.5 1.5]);

xlabel(sprintf('$ \\Delta \\eta^%s$', str),'interpreter','latex');
ylabel('MC / data','interpreter','latex');
set(gca,'box','on');

% Print
outputfile = sprintf('../figs/sim%d_pt_%0.0f.pdf', k, cuts.ptcut*1e3);
cmd = sprintf('print -dpdf %s', outputfile);
eval(cmd);
system(sprintf('pdfcrop --margins 10 %s %s', outputfile, outputfile));


end