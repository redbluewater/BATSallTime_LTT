function [ ML ] = get_mld_ctd( CTD )
%compute Mixed Layer Depth for a CTD profile by various methods
% and return a struct with a field for each criterion. 
% 
%  sig10 is computed as average of all values of sig_theta <= 10 meters
%  ML
%     .dens125 = depth where density >  sig10  + 0.125 
%     .bvfrq =  depth where bvfrq >  stddev(bvfrq) & de > 5m  Mojica and Gaube(2020 formulation)
%     .densT2 =  depth where density > sig10 + 0.2 * alpha
%                 (thermal expansion coefficient)
%
% Ruth Curry, BIOS / ASU
% Uploaded for BIOS-SCOPE project 19 October 2023
%

ML = struct();
ML.bvfrq = -999;
ML.dens125 = -999;
ML.densT2 = -999;

iz10 = find(CTD.Depth <= 10,1,'last');

if isempty(iz10)
    return;
end

if ~isfield(CTD,'bvfrq')
   BVFRQ = sw_bfrq(CTD.sa(:),CTD.te(:),CTD.pr(:),CTD.lat);
   CTD.bvfrq = [NaN; BVFRQ];
end

sd_nsq = nanstd(CTD.bvfrq(:));   %ML_bvfrq criteria
sig10 = nanmean(CTD.sig_theta(1:iz10));
alpha = sw_alpha(nanmean(CTD.Salt(1:iz10)),nanmean(CTD.Temp(1:iz10)),nanmean(CTD.Pressure(1:iz10)),'temp');

D125_crit = 0.125;
TE_crit = 0.2 * alpha * (sig10+1000) ;
  indx=find(CTD.bvfrq >= sd_nsq & CTD.Depth > 5,1,'first');
  if ~isempty(indx)
      ML.bvfrq = CTD.Depth(indx);
  end
  
  indx = find(CTD.sig_theta > sig10+D125_crit,1);
  if ~isempty(indx)
      ML.dens125 = CTD.Depth(indx);
  end
  
    indx = find(CTD.sig_theta > sig10+TE_crit,1);
  if ~isempty(indx)
      ML.densT2 = CTD.Depth(indx);
  end

end  %function

