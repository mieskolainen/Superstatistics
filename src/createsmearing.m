% Calculate smearing matrix based on source -> target mappings
%
% mikael.mieskolainen@cern.ch, 2019

function W = createsmearing(MC_targets_pu, MCtest_targets, C)

% Two-dimensional histogram
bins = 1:size(C,1);

W0 = hist3([MC_targets_pu MCtest_targets(1:length(MC_targets_pu))], {bins, bins})';

% Normalize each row
W = W0;
for i = 1:size(W0,1)
    W(i,:) = W0(i,:) / (sum(W0(i,:)) + 1e-9);
end

% Finally transpose
W = W';

% Fold data (for a test)
% y = W*x;

% W = [0.8400    0       0
%      0         0.84    0
%      0.1600    0.16    1.0000];

end