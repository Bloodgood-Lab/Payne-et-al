% Shannon information, sparseness, and selectivity
function [information,sparsity,selectivity] = getMapStats(map,posPDF);
% Written by Geoff Diehl with modifications by Anja Payne
% Last Modified: 10/26/2018

% Incoming Information: 
%   1) The rate map
%   2) The PDF of the time map

% Returns:
%   1) The spatial information for the cell
%   2) The sparsity of the cell
%   3) The selectivity of the cell

n = size(map,1);
meanrate = nansum(nansum( map .* posPDF ));
meansquarerate = nansum(nansum( (map.^2) .* posPDF ));
if meansquarerate == 0
    sparsity = NaN;
else
    sparsity = meanrate^2 / meansquarerate;
end
maxrate = max(max(map));
if meanrate == 0;
    selectivity = NaN;
else
    selectivity = maxrate/meanrate;
end
[i1, i2] = find( (map>0) & (posPDF>0) );  % the limit of x*log(x) as x->0 is 0
if length(i1)>0
    akksum = 0;
    for i = 1:length(i1);
        ii1 = i1(i);
        ii2 = i2(i);
        akksum = akksum + posPDF(ii1,ii2) * (map(ii1,ii2)/meanrate) * log2( map(ii1,ii2) / meanrate );
    end
    information = akksum;
else
    information = NaN;
end