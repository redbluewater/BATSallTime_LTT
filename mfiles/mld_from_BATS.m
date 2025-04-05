function [index, threashold]=mld(x, dT)

% function [index, threashold] = mld(x, dT)
%
% function to return the MLD index where:
%	index = row index for MLD
%	threashold = calculated value
%	x = x(p,t,s,sigmat)
%	dT =  degrees Celcius
%
% MLD => Sigmat(Zmld) >= Sigmat(avg(5-10m)) + dT*alpha
% alpha = coef thermal expansion
%	= a(s(avg(5-10m)), t(avg(5-10m)), p(0))
%
% Uses:
% calab.m - calculate alpha (Currently *only* for P==0).
% (calls: rhostp.m - calculate the density of seawater )
%
% Created: JCS - 21APR93

% dT=0.5;		% threashold temperature change
			% passed from calling routine for now.

s0 = median(x(1:5,3));
t0 = median(x(1:5,2));
sigmat0 = median(x(1:5,4));

a0 = calab(s0,t0,0);	% calculate the alpha from medians
threashold=[sigmat0+a0*dT*(sigmat0+1000)];

index=min(find(x(:,4) >= threashold));
if isempty(index)     % if no points found set index to zero == TRAPFLAG
      index=[-9.99];
else
 index=x(index,1);
end
