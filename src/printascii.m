% Print integer vector in ascii format
%
% input: vec = input vector
%
% mikael.mieskolainen@cern.ch, 2019

function printascii(vec)

for i = 1:length(vec)-1
    fprintf('%d,', vec(i));
end
fprintf('%d\n', vec(end));

end