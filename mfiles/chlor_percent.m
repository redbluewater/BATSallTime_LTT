function [I, chl_pct] = chlor_percent(chlor, imaxz)
% returns the index of the chlorophyl maximum in each column of chlor, and creates a property
% chl_pct reflecting percentage of the maximum
% Input 
%    chlor:  [nz x nprof]
%    imaxz: indx of max level to evaluate chlor
% 
%      I  :  index of the dcm  [1 x nprof]
% chl_pct :  array [nz x nprof]
%

chl_pct = ones(size(chlor)).* NaN;
[dcm_val, I] = max(chlor);
for ii= 1:length(I)
    imax = I(ii);
    if isnan(dcm_val(ii))  % check this because max(NaN) returns I=1
        I(ii) = NaN;
    else
        indx = find(~isnan(chlor(:,ii)),1,'last');
        if ~isnan(imax) && ~isempty(indx)
            theMax = chlor(imax,ii);
            chl_pct(:,ii) = chlor(:,ii) ./ theMax;  % range 0 -> 1
            chl_pct(imaxz:end,ii) = 0;
        end
    end
    
end


end %function