function [Zout] = label_vertical_zones_ctd(CTD,DCM,mld)
% function [Zout] = label_vertical_zones_ctd(CTD,DCM,mld)
%  Returns a column vector that assigns each element to 
%  a vertical zone  [0 .. 10]  for BIOSSCOPE 
%
%  0 :  Mixed layer: from sea surface to value supplied in mld
%  1 :  USEZ (upper stratified euphotic zone) :  from mld to top of DCM layer
%  2 :  DCM layer : defined by chlorophyll (fluorometer) values above 0.35 * max measured value
%  3 :  WMW layer  (Winter Mode Water) :  base of DCM to sigma_theta = 26.32 
%  4 :  sigma_theta 26.32 - 26.50
%  5 :  sigma_theta 26.50 - 26.70
%  6 :  sigma_theta 26.70 - 26.90
%  7 :  sigma_theta 26.90 - 27.10
%  8 :  sigma_theta 27.10 - 27.30  aka the Deep Ox Minimum
%  9 :  sigma_theta 27.30 - 27.50
% 10 :  anything denser
% -999 : no zone identified 
% KL note: can have no zone identified in cases where there is no
% fluorescence data because some zones depend on those data to give a zone
% (e.g., zone 2)
%
%  Note that zones 1 & 2 will be absent if the top of DCM is within the Mixed Layer
%
%  INPUT
%      CTD : struct representing a single profile
%      DCM : struct with fields describing DCM layer:  from get_dcm_layer_ctd()
%      mld : scalar value of mixed layer depth
%  OUTPUT
%      Zout : vector of vertical zones
%
% Ruth Curry, BIOS / ASU
% Uploaded for BIOS-SCOPE project 19 October 2023
%
%%
  Zout = ones(size(CTD.Depth)) .* -999;
  nz = length(Zout);
  if mld < -990 %too many messages, shut off the warning
      %disp(['NO MIXED LAYER defined for  ',num2str(CTD.Cruise(1)),' - ',num2str(CTD.Cast(1))]);
      return
  end
% 
 zone = 0;  % Mixed layer
 indx = find(CTD.Depth <= mld);
 if ~isempty(indx)
     Zout(indx) = zone;
 end
%
 zone = 1;  % USEZ layer
 if ~isnan(DCM.itop)
     indx = find(Zout~=0 & CTD.Depth <= DCM.de_top);
     if ~isempty(indx)
         Zout(indx) = zone;
     end
 end
 %
 zone = 2;  % DCM layer
 indx = find(CTD.Depth >= DCM.de_top & CTD.Depth <= DCM.de_bot);
 if ~isempty(indx)
     Zout(indx) = zone;
 end
 %
 zone = 3;  % WMW layer
 sigmax = 26.32;
 itop = find(CTD.Depth > DCM.de_bot,1);
 ibot = find(CTD.sig_theta <= sigmax,1,'last');
 if ~isempty(itop) && ~isempty(ibot) && ibot >= itop
     Zout(itop:ibot) = zone;
 end
 %
 zone = 4;
 sigmin = 26.32;
 sigmax = 26.50;
 itop = find(CTD.sig_theta > sigmin, 1);
 ibot = find(CTD.sig_theta <= sigmax, 1, 'last');
 if ~isempty(itop) 
      if isempty(ibot)
          ibot = nz;
      end
      Zout(itop:ibot) = zone;
 end
%
 zone = 5;
 sigmin = 26.50;
 sigmax = 26.70;
 itop = find(CTD.sig_theta > sigmin, 1);
 ibot = find(CTD.sig_theta <= sigmax, 1, 'last');
 if ~isempty(itop) 
      if isempty(ibot)
          ibot = nz;
      end
      Zout(itop:ibot) = zone;
 end
%
 zone = 6;
 sigmin = 26.70;
 sigmax = 26.90;
 itop = find(CTD.sig_theta > sigmin, 1);
 ibot = find(CTD.sig_theta <= sigmax, 1, 'last');
 if ~isempty(itop) 
      if isempty(ibot)
          ibot = nz;
      end
      Zout(itop:ibot) = zone;
 end
%
 zone = 7;
 sigmin = 26.90;
 sigmax = 27.10;
 itop = find(CTD.sig_theta > sigmin, 1);
 ibot = find(CTD.sig_theta <= sigmax, 1, 'last');
  if ~isempty(itop) 
      if isempty(ibot)
          ibot = nz;
      end
      if ibot >= itop
         Zout(itop:ibot) = zone;
      end
 end
%
 zone = 8;
 sigmin = 27.10;
 sigmax = 27.30;
 itop = find(CTD.sig_theta > sigmin, 1);
 ibot = find(CTD.sig_theta <= sigmax, 1, 'last');
  if ~isempty(itop)
     if isempty(ibot)
         ibot= nz;
     end
     if ibot >= itop
        Zout(itop:ibot) = zone;
     end
 end
%
 zone = 9;
 sigmin = 27.30;
 sigmax = 27.50;
 itop = find(CTD.sig_theta > sigmin, 1);
 ibot = find(CTD.sig_theta <= sigmax, 1, 'last');
 if ~isempty(itop)
   if isempty(ibot)
       ibot= nz;
   end
   if ibot >= itop
       Zout(itop:ibot) = zone;
   end
end
%
 zone = 10;
 indx = find(CTD.sig_theta > sigmax);
  if ~isempty(indx) 
     Zout(indx) = zone;
 end

end % function
