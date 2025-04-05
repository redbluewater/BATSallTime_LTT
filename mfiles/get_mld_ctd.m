function [ ML ] = get_mld_ctd( CTD )
%compute Mixed Layer Depth for a CTD profile by various methods
% and return a struct with a field for each criterion. 
% 
%  sig10 is computed as average of all values of sig_theta <= 10 meters
%  ML
%     .bvfrq =  depth where bvfrq >  stddev(bvfrq) & de > 5m  Mojica and Gaube(2020 formulation)
%     .dens125 = depth where density >  sig10  + 0.125 
%     .dens2 = depth where density >  sig10  + 0.2  ('mld6' on Craig's list)
%     .densT1 =  depth where density > sig10 + 0.1 * alpha 
%     .densT2 =  depth where density > sig10 + 0.2 * alpha
%                 (thermal expansion coefficient); this is Sprintall and
%                 Tomczak (1992); equivalent to MLD_BATS
%     .densT3 =  depth where density > sig10 + 0.3 * alpha 
%     .densGR = depth found by searching for a change of 
%                  0.001 (kg/m3)/m in sigma_t profile ('mld1' on Craig's
%                  list, name densGR because this is searching for a
%                  gradient
%     .te2 = depth where T < SST - 0.2C ('mld5' on Craig's list) 

%
% Ruth Curry, BIOS / ASU
% Uploaded for BIOS-SCOPE project 19 October 2023
% Krista Longnecker editing to add more MLD definitions 24 June 2024

%KL adding options 24 June 2024
ML = struct();
ML.bvfrq = -999;
ML.dens125 = -999;
ML.dens2 = -999;
ML.densT1 = -999;
ML.densT2 = -999;
ML.densT3 = -999;
ML.densGR = -999;
ML.te2 = -999;

%only move forward if we have CTD data from shallow depths (have at least
%one cruise/cast where the data do not start until 2000m (cruise 10002,cast 4)
% if isequal(CTD.BATS_id(1),'10001003')
%     keyboard
% end
iz10 = find(CTD.Depth <= 10,1,'last'); %find the depth where CTD data are shallower than 10 m to do this calculation

% put a check here, cannot calculate variables if sensor is not working
noSalt = 0;
if isnan(CTD.Salt(1:iz10))
    noSalt = 1;
end

if isempty(iz10)
    return;
end

%This calculation has already been done, comment this out KL 6/24/2024
% if ~isfield(CTD,'bvfrq')
%    BVFRQ = sw_bfrq(CTD.sa(:),CTD.te(:),CTD.pr(:),CTD.lat);
%    CTD.bvfrq = [NaN; BVFRQ];
% end

%do the math for the various options for MLD
sd_nsq = nanstd(CTD.bvfrq(:));   %ML_bvfrq criteria
sig10 = nanmean(CTD.sig_theta(1:iz10));
te10 = nanmean(CTD.Temp(1:iz10));
alpha = sw_alpha(nanmean(CTD.Salt(1:iz10)),nanmean(CTD.Temp(1:iz10)),nanmean(CTD.Pressure(1:iz10)),'temp');
%search for a change of 0.001 (kg/m3)/m in sigma_t profile:	
sigDiff = [diff(CTD.sig_theta)./diff(CTD.Depth)]; %need both diffs here (sig_theta and depth)

%these are the values that are searched for in the steps below
D125_crit = 0.125;
D2_crit = 0.2;
TE_crit = 0.2 * alpha * (sig10+1000) ; %KL note: TE is the thermal expansion coefficient
TE_critOne = 0.1 * alpha * (sig10+1000) ; %spell out 'one' to make this easier to read
TE_crit3 = 0.3 * alpha * (sig10+1000) ; 
dTemp = te10 - 0.2;
sDiff = 0.001;

%now go through the various options to find the depth matching each
%criteria
  indx=find(CTD.bvfrq >= sd_nsq & CTD.Depth > 5,1,'first');
  if ~isempty(indx) && ~noSalt
      ML.bvfrq = CTD.Depth(indx);
  end
  
  indx = find(CTD.sig_theta > sig10+D125_crit,1);
  if ~isempty(indx) && ~noSalt
      ML.dens125 = CTD.Depth(indx);
  end
  
  indx = find(CTD.sig_theta > sig10+D2_crit,1);
  if ~isempty(indx) && ~noSalt
      ML.dens2 = CTD.Depth(indx);
  end
  
  indx = find(CTD.sig_theta > sig10+TE_critOne,1);
  if ~isempty(indx) && ~noSalt
      ML.densT1 = CTD.Depth(indx);
  end
  
  indx = find(CTD.sig_theta > sig10+TE_crit,1);
  if ~isempty(indx) && ~noSalt
      ML.densT2 = CTD.Depth(indx);
  end

  indx = find(CTD.sig_theta > sig10+TE_crit3,1);
  if ~isempty(indx) && ~noSalt
      ML.densT3 = CTD.Depth(indx);
  end
  
  indx = find(CTD.Temp < dTemp,1);
  if ~isempty(indx) && ~noSalt
      ML.te2 = CTD.Depth(indx);
  end
   
  indx = find(abs(sigDiff)> sDiff,1,'first');
  if ~isempty(indx) && ~noSalt
      ML.densGR = CTD.Depth(indx);
  end
  
end  %function

