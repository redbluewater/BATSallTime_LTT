%
% Ruth Curry, BIOS / ASU
% Uploaded for BIOS-SCOPE project 19 October 2023
%
%  edit bad segments of fluorometer profiles and overwrite individual *.csv
%  files and CRU*.csv files; save corrected *.mat files

rootdir = '/Users/rcurry/GliderData/BIOSSCOPE/CTD_BOTTLE';
workdir = fullfile(rootdir,'00_CTD');

% read in MAT file
cd(workdir)
MATfile = fullfile(workdir,'20190817_10362_CTD.mat');
load(MATfile);

%   read in CRU*.csv file
fmt = '/Users/rcurry/GliderData/BIOSSCOPE/CTD_BOTTLE/00_CTD/CRU_%1d%04d_ctd.csv';
CRU_fname = sprintf(fmt,CTD.type(1),CTD.cruise(1));
CRUtab = readtable(CRU_fname);
CRU = table2struct(CRUtab,'ToScalar',true);
clear CRUtab

%  define cast directory
fmt = [workdir,'/b%1d%04d/'];
castdir = sprintf(fmt,CTD.type(1),CTD.cruise(1));
%
figure; 
   hold on; axis ij;
    for iprof = 1:23
        plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k');
        plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-')
    end
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    plot([0. 0.],ylim(),'--g');
ibad = [];   % 4     9    19    20    21    22    23
for iprof = 1:23
    if any(find(CTD.fluor_filt(:,iprof) > 0.5 | CTD.fluor_filt(:,iprof) < -0.01))
        ibad = [ibad; iprof];
    end
end
%%%%%%%%%%%%%%%%%%%%%%

% iprof = ibad(1);  % 4 is OK
%     figure; hold on; axis ij;
%     plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
%     plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-r','Linewidth',1.5);
%     plot([-0.01 -0.01],ylim(),'--b');
%     plot([0.5 0.5],ylim(),'--b');
%     title(['Cast #: ' num2str(iprof)]);
%     plot([0. 0.],ylim(),'--g','Linewidth',1.5);
%%%%%%%%%%%%%%%%%%%%%%

iprof = ibad(2);  % 9
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.01 -0.01],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);

    dmin = 425;
    dmax = 451;
    indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
    xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof));
    CTD.fluor(indx,iprof) = xnew;

    subplot(1,2,2); hold on; axis ij;
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    XX.Pressure = CTD.pr(:,iprof);
    XX.Fluor = CTD.fluor(:,iprof);
        filt_width = 3;
    [Xfilt,bias] = smooth_fluor(XX,filt_width);
    clear XX
        Xfilt = Xfilt - bias;
        CTD.fluor_offset(iprof) = bias;
    CTD.fluor_filt(:,iprof) = Xfilt;
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-c','Linewidth',1.5);
    
    percent = 0.33;
    izmax = find(CTD.de(:,iprof) > 400,1);
    if isnan(izmax)
        izmax = find(~isnan(CTD.de(:,prof)),1,'last');
    end
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);
    if isnan(DCM.depth)
        DCM.depth = -999.0;
    end
    CTD.DCM(iprof) = DCM.depth;
    
%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
   
   nz = length(XX.Cast);
   XX.DCM(:) = DCM.depth;
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);
       
   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   

   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;
   
   % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;
   
   clear XX TTcast DCM
%%%%%%%%%%%%%%%%%%%%%
% iprof = ibad(3);  % 19  is OK -- just at top
%     figure; 
%     subplot(1,2,1); hold on; axis ij;
%     xlim([-10 10]);
%     plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
%     plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
%     plot([-0.02 -0.02],ylim(),'--b');
%     plot([0.5 0.5],ylim(),'--b');
%     title(['Cast #: ' num2str(iprof)]);
%     plot([0. 0.],ylim(),'--g','Linewidth',1.5);
%     xlim([-0.2 0.5]);

iprof = ibad(4);  % 20
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);
    xlim([-5 0.5]);
    
    dmin = 455;
    dmax = 470;
    indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
    xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof));
    CTD.fluor(indx,iprof) = xnew;

    dmin = 230;
    dmax = 240;
    indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
    xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof))
    CTD.fluor(indx,iprof) = xnew;
    
    subplot(1,2,2); hold on; axis ij;
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    XX.Pressure = CTD.pr(:,iprof);
    XX.Fluor = CTD.fluor(:,iprof);
        
        filt_width = 3;
    [Xfilt,bias] = smooth_fluor(XX,filt_width);
     
    if abs(bias) > 0.005 
        Xfilt = Xfilt - bias;
        CTD.fluor_offset(iprof) = bias;
    end
    CTD.fluor_filt(:,iprof) = Xfilt;
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-c','Linewidth',1.5);
    
    percent = 0.33;
    izmax = find(CTD.de(:,iprof) > 400,1);
    if isnan(izmax)
        izmax = find(~isnan(CTD.de(:,prof)),1,'last');
    end
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);
    CTD.DCM(iprof) = DCM.depth;
%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
         
   nz = length(XX.Cast);
   XX.DCM(:) = DCM.depth;
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);
   
   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   

   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;

   % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;

   clear XX TTcast DCM
%%%%%%%%%%%%%%%%%%%%%

iprof = ibad(5);  % 21
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);
    xlim([-5 0.5]);

    dmin = 394;
    dmax = 421;
        indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
        xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof));
        CTD.fluor(indx,iprof) = xnew;
    
    subplot(1,2,2); hold on; axis ij;
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    XX.Pressure = CTD.pr(:,iprof);
    XX.Fluor = CTD.fluor(:,iprof);
        
        filt_width = 3;
    [Xfilt,bias] = smooth_fluor(XX,filt_width);
     
    if abs(bias) > 0.005 
        Xfilt = Xfilt - bias;
        CTD.fluor_offset(iprof) = bias;
    end
    CTD.fluor_filt(:,iprof) = Xfilt;
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-c','Linewidth',1.5);
    
    percent = 0.33;
    izmax = find(CTD.de(:,iprof) > 400,1);
    if isnan(izmax)
        izmax = find(~isnan(CTD.de(:,prof)),1,'last');
    end
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);
    CTD.DCM(iprof) = DCM.depth;
%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
   
   nz = length(XX.Cast);
   XX.DCM(:) = DCM.depth;
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);

   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   

   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;

   % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;
   
   clear XX TTcast DCM
    
%%%%%%%%%%%%%%%%%%%%%

iprof = ibad(6);  % 22
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);
    xlim([-5 0.5]);
    
    dmin = 315;
    dmax = 328;
    indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
    xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof))
    CTD.fluor(indx,iprof) = xnew;
    
    dmin = 800;
    dmax = 940;
        indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
        xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof))
        CTD.fluor(indx,iprof) = xnew;
    
    subplot(1,2,2); hold on; axis ij;
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    XX.Pressure = CTD.pr(:,iprof);
    XX.Fluor = CTD.fluor(:,iprof);
        
        filt_width = 3;
    [Xfilt,bias] = smooth_fluor(XX,filt_width);
     
    if abs(bias) > 0.005 
        Xfilt = Xfilt - bias;
        CTD.fluor_offset(iprof) = bias;
    end
    CTD.fluor_filt(:,iprof) = Xfilt;
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-c','Linewidth',1.5);
    
    percent = 0.33;
    izmax = find(CTD.de(:,iprof) > 400,1);
    if isnan(izmax)
        izmax = find(~isnan(CTD.de(:,prof)),1,'last');
    end
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);
    CTD.DCM(iprof) = DCM.depth;
%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
   
   nz = length(XX.Cast);
   XX.DCM(:) = DCM.depth;
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);
      
   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   
   
   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;

   % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;
   
   clear XX TTcast DCM
%%%%%%%%%%%%%%%%%%%%%
% iprof = ibad(7);  % 23  This had a duplicate cast in the original
% *ctd.txt file which I removed
%%%%%%%%%%%%%%%%%%

% Overwrite CRU*.csv 
disp(['Overwriting ',CRU_fname]);
TTcast = struct2table(CRU);
writetable(TTcast,CRU_fname);

% Save Matfile
disp(['Saving ',MATfile]);
save(MATfile,'CTD');

%%
% read in MAT file
cd(workdir)
MATfile = fullfile(workdir,'20190908_10363_CTD.mat');
load(MATfile);

%   read in CRU*.csv file
fmt = '/Users/rcurry/GliderData/BIOSSCOPE/CTD_BOTTLE/00_CTD/CRU_%1d%04d_ctd.csv';
CRU_fname = sprintf(fmt,CTD.type(1),CTD.cruise(1));
CRUtab = readtable(CRU_fname);
CRU = table2struct(CRUtab,'ToScalar',true);
clear CRUtab

%  define cast directory
fmt = [workdir,'/b%1d%04d/'];
castdir = sprintf(fmt,CTD.type(1),CTD.cruise(1));
%
nprof = length(CTD.cast);
figure; 
   hold on; axis ij;
    for iprof = 1:nprof
        plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k');
        plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-')
    end
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    plot([0. 0.],ylim(),'--g');
ibad = [];   % 2    13    15    18    22
for iprof = 1:nprof
    if any(find(CTD.fluor_filt(:,iprof) > 0.5 | CTD.fluor_filt(:,iprof) < -0.01))
        ibad = [ibad; iprof];
    end
end

%%%%%%%%%%%%%
iprof = ibad(1);  % 2
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);
    xlim([-5 5]);
    
    %  Entire profile is bad
    CTD.fluor(:,iprof) = NaN;
    CTD.fluor_filt(:,iprof) = NaN;
    CTD.fluor_offset(iprof) = NaN;
    CTD.DCM(iprof) = NaN;
    
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);
%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
   
   nz = length(XX.Cast);
   XX.DCM(:) = CTD.DCM(iprof);
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);
   
   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   

   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;

   % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;
  
   clear XX TTcast DCM
%%%%%%%%%%%%%
iprof = ibad(2);  % 13
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);
    xlim([-5 5]);
    
    %  Upper 100m of profile is bad -- NaN entire cast
    CTD.fluor(:,iprof) = NaN;
    CTD.fluor_filt(:,iprof) = NaN;
    CTD.fluor_offset(iprof) = NaN;
    CTD.DCM(iprof) = NaN;
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);
%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
   
   nz = length(XX.Cast);
   XX.DCM(:) = CTD.DCM(iprof);
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);

   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   
   
   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;

   % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;
   
   clear XX TTcast DCM

   %%%%%%%%%%%%%
iprof = ibad(3);  %15 
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);
    xlim([-5 5]);

    %  Entire profile is bad
    CTD.fluor(:,iprof) = NaN;
    CTD.fluor_filt(:,iprof) = NaN;
    CTD.fluor_offset(iprof) = NaN;
    CTD.DCM(iprof) = NaN;
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);

%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
   
   nz = length(XX.Cast);
   XX.DCM(:) = CTD.DCM(iprof);
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);
   
   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   

   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;

   % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;
   
   clear XX TTcast DCM
%%%%%%%%%%%%%
% iprof = ibad(4);  %18 Seems okay
%     figure; 
%     subplot(1,2,1); hold on; axis ij;
%     xlim([-10 10]);
%     plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
%     plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
%     plot([-0.02 -0.02],ylim(),'--b');
%     plot([0.5 0.5],ylim(),'--b');
%     title(['Cast #: ' num2str(iprof)]);
%     plot([0. 0.],ylim(),'--g','Linewidth',1.5);
%     xlim([-5 5]);
%%%%%%%%%%%%%

iprof = ibad(5);  % 
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);
    xlim([-.05 .1]);

    dmin = 240;
    dmax = 310;
    indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
    xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof))
    CTD.fluor(indx,iprof) = xnew;

    dmin = 365;
    dmax = 385;
    indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
    xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof))
    CTD.fluor(indx,iprof) = xnew;

    subplot(1,2,2); hold on; axis ij;
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    XX.Pressure = CTD.pr(:,iprof);
    XX.Fluor = CTD.fluor(:,iprof);
        
    filt_width = 3;
    [Xfilt,bias] = smooth_fluor(XX,filt_width);
     
    if abs(bias) > 0.005 
        Xfilt = Xfilt - bias;
        CTD.fluor_offset(iprof) = bias;
    end
    CTD.fluor_filt(:,iprof) = Xfilt;
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-c','Linewidth',1.5);

    percent = 0.33;
    izmax = find(CTD.de(:,iprof) > 400,1);
    if isnan(izmax)
        izmax = find(~isnan(CTD.de(:,prof)),1,'last');
    end
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);
    CTD.DCM(iprof) = DCM.depth;
%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
   
   nz = length(XX.Cast);
   XX.DCM(:) = DCM.depth;
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);
   
   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   

   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;
   
 % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;
   
   clear XX TTcast DCM
%%%%%%%%%%%%%%%%%%%%%%%%%   
iprof = 14;  %  has a spike not detected
    figure; 
    subplot(1,2,1); hold on; axis ij;
    xlim([-10 10]);
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'r','Linewidth',1);
    plot([-0.02 -0.02],ylim(),'--b');
    plot([0.5 0.5],ylim(),'--b');
    title(['Cast #: ' num2str(iprof)]);
    plot([0. 0.],ylim(),'--g','Linewidth',1.5);
    xlim([-.05 .1]);

    dmin = 20;
    dmax = 35;
    indx = find(CTD.de(:,iprof) >= dmin & CTD.de(:,iprof) <= dmax);
    xnew = 0.5 * ( CTD.fluor(indx(1),iprof) + CTD.fluor(indx(end),iprof))
    CTD.fluor(indx,iprof) = xnew;


    subplot(1,2,2); hold on; axis ij;
    plot(CTD.fluor(:,iprof),CTD.de(:,iprof),'-k','Linewidth',1.5);
    XX.Pressure = CTD.pr(:,iprof);
    XX.Fluor = CTD.fluor(:,iprof);
        
    filt_width = 3;
    [Xfilt,bias] = smooth_fluor(XX,filt_width);
     
    if abs(bias) > 0.005 
        Xfilt = Xfilt - bias;
        CTD.fluor_offset(iprof) = bias;
    end
    CTD.fluor_filt(:,iprof) = Xfilt;
    plot(CTD.fluor_filt(:,iprof),CTD.de(:,iprof),'-c','Linewidth',1.5);

    percent = 0.33;
    izmax = find(CTD.de(:,iprof) > 400,1);
    if isnan(izmax)
        izmax = find(~isnan(CTD.de(:,prof)),1,'last');
    end
    DCM = get_dcm_layer_ctd(CTD.fluor_filt(:,iprof),CTD.de(:,iprof), CTD.MLD_dens125(iprof),percent,izmax);
    CTD.DCM(iprof) = DCM.depth;
%%%
%  Read in cast file and update fields
   fmt = '%4d%02d%02d_%1d%04d_%03d_ctd.csv';
   CAST_fname = sprintf(fmt,CTD.year(iprof), CTD.month(iprof), CTD.day(iprof),CTD.type(iprof),CTD.cruise(iprof),CTD.cast(iprof));
   TTin = readtable([castdir,CAST_fname]);
   XX = table2struct(TTin,'ToScalar',true);
   clear TTin
   
   nz = length(XX.Cast);
   XX.DCM(:) = DCM.depth;
   XX.Fluor_offset(:) =  CTD.fluor_offset(iprof);
   XX.Fluor_filt(:) = CTD.fluor_filt(1:nz,iprof);
   XX.Fluor(:) = CTD.fluor(1:nz,iprof);
   
   XX.Fluor_filt(isnan(XX.Fluor_filt)) = -999.0;   % substitute missing values
   XX.Fluor(isnan(XX.Fluor)) = -999.0;   
   XX.Fluor_offset(isnan(XX.Fluor_offset)) = -999.0;   
   XX.DCM(isnan(XX.DCM)) = -999.0;   

   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,XX.MLD_dens125);
   CTD.vertZone(1:length(XX.VertZone),iprof) = XX.VertZone(:);
   CTD.vertZone(CTD.vertZone(:,iprof)< -990,iprof) = NaN;

   % Output the cast
   TTcast = struct2table(XX);
   writetable(TTcast,[castdir,CAST_fname]);
   
   % Update CRU fields
   icru = find(CRU.Cast == CTD.cast(iprof));
   if isempty(icru)
       error(['Unable to find cast ',num2str(CTD.cast(iprof)),' in ',CRU_fname]);
   end
   
   CRU.DCM(icru) = XX.DCM;
   CRU.Fluor_offset(icru) = XX.Fluor_offset;
   CRU.Fluor_filt(icru) = XX.Fluor_filt;
   CRU.Fluor(icru) = XX.Fluor;
   CRU.VertZone(icru) = XX.VertZone;
   
   clear XX TTcast DCM
%%%%%%%%%%%%%%%%%%%%%%   
% Overwrite CRU*.csv 
disp(['Overwriting ',CRU_fname]);
TTcast = struct2table(CRU);
writetable(TTcast,CRU_fname);

% Save Matfile
disp(['Saving ',MATfile]);
save(MATfile,'CTD');










