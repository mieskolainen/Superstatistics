% PDF print cropping
%
% input = [hor,ver] [0 ... 1] scale factors
%
% mikael.mieskolainen@cern.ch, 2019

function pdfcrop(hor, ver)

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3)*hor fig_pos(4)*ver]; % horizontal vertical

end