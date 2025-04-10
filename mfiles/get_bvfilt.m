function [ bvfilt ] = get_bvfilt(bvfreq,filtwidth  )
%  function [ bvfilt ] = get_bvfilt(bvfreq,filtwidth  )
% Applies a 3rd order butterworth filter to a rectangular array of 
% buoyancy frequency profiles, and returns the filtered profiles.
%  filtwidth is the number of points in the filter window
% code from Ruth (June 2020)
%in Ruth's 2020 sample code she called this as follows: bvfilt = get_bvfilt(bvfrq,5)
[b,a] = butter(3,1/filtwidth);
bvfilt = ones(size(bvfreq)).* NaN;

[~,nprof] = size(bvfreq);
for ii = 1:nprof
    indx = find(~isnan(bvfreq(:,ii)));
    if length(indx) > 10
        bvfilt(indx,ii) = filtfilt(b,a,bvfreq(indx,ii));
    end
end
end

