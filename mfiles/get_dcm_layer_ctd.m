function [DCM] = get_dcm_layer_ctd(CHL, DE, MLD, percent, izmax)
% function [DCM] = get_dcm_layer(CHL, DE, MLD, percent, izmax)
% Evaluates the layer associated with the deep chlorophyll max (DCM)
% as the portion of the profile falling within the specified % of the max
% value
%
%  Input:  
%   CHL:  array [nz, nprof]
%   DE:  array [nz, nprof]
%   MLD :  row vector [1,nprof]  MLD for each profile
%   percent : of chlor max value to define layer (scalar)
%   izmax :  index of max zlevel to evaluate (scalar)
%
%  Output:
%   DCM with fields (all row vectors [1,nprof]
%       .depth : depth of DCM
%       .chlor_val : value of chlor at DCM
%       .itop  :  index of top of layer [1, nprof]  : NaN if DCM is in ML
%       .ibot  :  index of bottom of layer   == NaN
%       .de_top :  depth at top == NaN if DCM is within the ML
%       .de_bot : depth at bottom == MLD if DCM is within the ML
%       .DCMinML : 1 = yes; 0 = no (true/false, is the DCM in the ML)
%
%
% Ruth Curry, BIOS / ASU
% Uploaded for BIOS-SCOPE project 19 October 2023
% Krista Longnecker 2 July 2024 (saved Ruth's version as get_dcm_layer_ctd_Ruth_v0.m)
% adding a flag if the top or bottom of the DCM is in the ML, but change to 
% keeping the actual DCM value; also do not set DCM to MLD depth...use
% logic later to track this

[~,nprof] = size(CHL);
DCM = struct();
XX = ones(1,nprof) .* NaN;
DCM.depth = XX;
DCM.chlor_val = XX;
DCM.itop = XX;
DCM.ibot = XX;
DCM.de_top = XX;
DCM.de_bot = XX;
DCM.DCMinML = zeros(1,nprof); %KL adding 7/2/2024
clear XX

[dcm_indx, chl_pct] = chlor_percent(CHL, izmax);
Ipct = chl_pct >= percent;

for ii = 1:nprof

    idcm = dcm_indx(ii);
    profdepth = max(DE(:,ii));
    if ~isnan(idcm)
        if profdepth < 100 || profdepth == DE(idcm,ii)
            disp('Profile too shallow to evaluate DCM');
            continue
        end
        
        DCM.depth(ii) = DE(idcm,ii);
        DCM.chlor_val(ii) = CHL(idcm,ii); 
 
        % DCM is within ML , layer is not defined
        if DCM.depth(ii) < MLD(ii)  
            DCM.ibot(ii) = NaN;   
            DCM.itop(ii) = NaN;  
            DCM.de_top(ii) = MLD(ii);  
            DCM.de_bot(ii) = MLD(ii);  
            DCM.DCMinML(ii) = 1; %KL adding 7/2/2024
        else
            %get top of layer
            indx = find(Ipct(:,ii) > 0,1,'first');
            if ~isempty(indx)
                DCM.itop(ii) = indx;
                DCM.de_top(ii) = DE(indx,ii);
              % top of DCM within ML, set to MLD and index to NaN %KL comment 7/2/2024 - testing
                if DCM.de_top(ii) < MLD(ii)   
                    %DCM.itop(ii) = NaN; %KL comment 7/2/2024 - testing
                    %DCM.de_top(ii) = MLD(ii);%KL comment 7/2/2024 - testing
                    DCM.DCMinML(ii) = 1; %KL adding 7/2/2024
                end
            end
            %
            %get bottom of layer
            indx = find(Ipct(:,ii) > 0,1,'last');
            if ~isempty(indx)
                DCM.ibot(ii) = indx;
                DCM.de_bot(ii) = DE(indx,ii);
                % bottom of DCM within ML; set to MLD and index to NaN %KL comment 7/2/2024 - testing
                if DCM.de_bot(ii) <= MLD(ii) 
                    %DCM.ibot(ii) = NaN;     %KL comment 7/2/2024 - testing
                    %DCM.de_bot(ii) = MLD(ii);%KL comment 7/2/2024 - testing
                    DCM.DCMinML(ii)= 1;
                end
            end
        end
    end
end %for ii

end %function
