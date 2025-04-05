function [ Xout ] = get_BATS_par_vars(PAR_in, de_in, zrange, zmax )
%  function [ Xout ] = get_BATS_par_vars(PAR_in, de_in, zrange, zmax )
%   
%  Fits an exponential curve to the PAR/depth profile over the interval
%  [zrange](usually upper 100 m) from sea surface down to a depth = zmax.
%  Missing values for BATS PAR  are set to -999.

%  Returns info in fields of Xout:
%    .par_est  :fitted profile
%    .kpar     : scalar value for each profile
%    .par0     : surface value
%    .z_1pcnt : depth of 1% light levels
%    .z_halfpcnt :  0.5%
%    .z_tenthpcnt : 0.1%
%
%  NOTE:  Missing values returned in Xout are set to NaN.
%% 19-Oct-2023 rcurry Updated to catch error in fit() of par profiles
%% KL turned off some warnings 14 June 2024

%% Initialize output struct

Xout = struct();
Xout.par_est = ones(size(PAR_in)).* NaN;
Xout.kpar = NaN;
Xout.par0 = NaN;
Xout.z_par_1pcnt = NaN;
Xout.z_par_halfpcnt = NaN;
Xout.z_par_tenthpcnt = NaN;

%  replace -999 values with NaN
PAR_in(PAR_in < -990) = NaN;
de_in(de_in < -990) = NaN;
notNan = ~isnan(PAR_in) & ~isnan(de_in);
if sum(notNan) < 30
    %disp('Not enough pts to fit PAR')
    return
end
%
indx = find(de_in <= zmax);
if max(PAR_in(notNan)) < 2
    Xout.par_est(indx) = 0;  %set profile to zero
    Xout.par0 = 0;
    return
end
%
igood = find(notNan & de_in >= zrange(1) & de_in <= zrange(2));
if length(igood) < 30
    %disp('Not enough points to fit PAR')
    return
end
%
% exponential fit 
 de_good = de_in(igood);
 PAR_good = PAR_in(igood);
 try
     [f, gof,options] = fit(de_good,PAR_good,'exp1','robust','bisquare');
    
 catch
     %warning('WARNING: fit() for par profile returned an error. Par variables set to NaN'); 
     return
 end
 
 % Compute par values
 Coef=coeffvalues(f);
%
 Xout.par_est(indx)= Coef(1)*exp(Coef(2)* de_in(indx)); 
 Xout.kpar = -Coef(2);
 Xout.par0= Coef(1); % 
%
 xx = 0.01 * Xout.par0;   % 1% light level
 iz = find(Xout.par_est < xx, 1);  
 if ~isempty(iz)
    Xout.z_par_1pcnt = de_in(iz);
 end
%
 xx = 0.005 * Xout.par0;   % 0.5% light level
 iz = find(Xout.par_est < xx, 1); 
 if ~isempty(iz)
    Xout.z_par_halfpcnt = de_in(iz);
 end
    %
 xx = 0.001 * Xout.par0;   % 0.1% light level
 iz = find(Xout.par_est < xx, 1);  
 if ~isempty(iz)
    Xout.z_par_tenthpcnt = de_in(iz);
 end

end


