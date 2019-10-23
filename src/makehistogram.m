% Histogram simulations
%
% mikael.mieskolainen@cern.ch, 2018

function [histo, chi2] = makehistogram(obs, DATA, param)

% Histogram binning setup
histo.binedges = [table2array(DATA.table(:,2)); table2array(DATA.table(end,3))];

% Take the observable
O = obs.observable{1};

% Histogram it
count_ND  = hist1m(O(O(:,2) == 1), histo.binedges);
count_SD  = hist1m(O(O(:,2) == 2), histo.binedges);
count_DD  = hist1m(O(O(:,2) == 3), histo.binedges);
count_tot = count_ND + count_SD + count_DD;

% Normalization [histogram efficiency] x [total fiducial factor] x [total inelastic]
% The histogram efficiency takes into account that not all fiducial events
% belong to the histogram
N = (sum(count_tot) / size(O,1)) * (obs.accepted / obs.trials) * param.sigmainel;

% Binwidth
N = N ./ (histo.binedges(2:end)-histo.binedges(1:end-1));

% Calculate histograms
histo.hMC(1,:) = N .* count_ND / sum(count_tot);
histo.hMC(2,:) = N .* count_SD / sum(count_tot);
histo.hMC(3,:) = N .* count_DD / sum(count_tot);
histo.hTOTAL   = sum(histo.hMC,1);

% GET DATA
histo.xPOINTS = table2array(DATA.table(:,1));
histo.hDATA   = table2array(DATA.table(:,4));

% Get uncertanties (stat, syst, luminosity) sum in quadrature
histo.hDATA_U = zeros(length(histo.xPOINTS),1);
histo.hDATA_D = zeros(length(histo.xPOINTS),1);

% Strip numbers from ascii
hDATA_cell_U = table2array(DATA.table(:,[5 7 9]));
hDATA_cell_D = table2array(DATA.table(:,[6 8 10]));
for ii = 1:size(hDATA_cell_U,1)
   val = zeros(3,1);
   for jj = 1:3
      val(jj) = sscanf(hDATA_cell_U{ii,jj}, '''%f'''); 
   end
   histo.hDATA_U(ii) = sqrt(sum(val.^2)); % quadrature
   
   val = zeros(3,1);
   for jj = 1:3
      val(jj) = sscanf(hDATA_cell_D{ii,jj}, '''%f''');
   end
   histo.hDATA_D(ii) = sqrt(sum(val.^2)); % quadrature
end

% chi^2
chi2 = sum( (1 - histo.hTOTAL(:) ./ histo.hDATA(:)).^2);

end