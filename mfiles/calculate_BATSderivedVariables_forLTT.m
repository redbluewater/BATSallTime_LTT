function Xout = calculate_BATSderivedVariables_forLTT(infile,trans_dates,do_plots,outdir,giveNotice)
% function Xout = calculate_BATSderivedVariables_forLTT(infile,trans_dates,do_plots,outdir,giveNotice)
% Reads a *mat file (from BATS group), creates a structure with added 
% fields, and writes a .csv file for each individual cruise. 
% INPUT:
%  infile is the name of the file to read.
%  do_plots is 1/0 depending if you do (1) or do not (0) want figures
%  outdir is where to put the final CSV files
%
% notes from Krista, this calls a script to generate the seasons based on
% pre-set date windows and does NOT rely on analysis of the water column
% data to change the transition dates each year
%
% OUTPUT:
%  Writes individual files for each cruise in csv format
%  Xout  :  a matlab structure with fields storing info for entire
%           cruise in row vectors and rectangular matrices (KL note this is
%           not currently being exported)

% NOTE:  Vertical zones are computed using ML_dens125
%
% Original script from Ruth Curry BIOS/ASU (2023); Krista Longnecker
% editing to add more MLD definitions 24 June 2024
% KL using this to get BV values for LTT paper; 4 April 2025
%KL added filter to BV data; 10 April 2025

%%  Read file into rectangular array, and store each column as a field in structure CTD
%KL version, read in the MATLAB file and then convert to match the legacy format
%first, get the name of the variables in the file
S = whos('-file',infile); 
Sn = {S.name}; 
%then pull the data from the structure, will be gf#####_data (where ##### is
%the BATS cruise id)
r = contains(Sn,'_data');
kr = find(r==1); %which part of the structure has the data?
riMAT = load(infile,Sn{kr});
%ugly MATLAB hack to pull out the data into a matrix and put into Ruth's
%format; TTin is a 1 x 12 cell, one variable per cell
v = fieldnames(riMAT);
TTin = riMAT.(v{1}); %this will read in matrix, only do this to get the number of rows
% # of rows is the same variable that Ruth was calling MAXZ
% MAXZ is the deepest water sample in this cruise (Ruth was setting this to
% 2500 m, but KL changed this to dynamically determine for each cast)
MAXZ = size(TTin,1);
%now use the number of rows to make an array matching Ruth's format
TTin = mat2cell(riMAT.(v{1}),[MAXZ],[ones(12,1)]); %12 is fixed - number of variables
%tidy up this mess
clear S Sn r kr riMAT v 

%%  Assign columns 
%KL note - no conductivity in the CTD data BATS data files, change the indexing
icol.cast_id = 1;
icol.dec_yr = 2;
icol.lat = 3;
icol.lon =4;
icol.pr = 5;
icol.de = 6;
icol.te = 7;
% icol.co = 8; %no conductivity in these data KL 2/8/2024
icol.sa = 8;
icol.o2 = 9;
icol.beam = 10;
icol.fluor = 11;
icol.par = 12;

% create structure to store all the casts     
CTD = struct();   % for conversion to csv files

CTD.BATS_id = TTin{icol.cast_id};
    yy = floor(CTD.BATS_id .* 1e-7);
    xx =  CTD.BATS_id - yy .* 1e7;
CTD.Cruise = floor(xx .* 1e-3);
CTD.Cast = xx - CTD.Cruise .* 1e3;
    clear xx yy

ZZ = ones(size(CTD.Cruise)) .* -999; %blank array

CTD.decy = TTin{icol.dec_yr};
    dvec = datevec(decyear2dnum(CTD.decy));
CTD.yyyymmdd = dvec(:,1).* 1e4 + dvec(:,2) .* 1e2 + dvec(:,3);
CTD.hhmm = dvec(:,4).*1e2 + dvec(:,5);
CTD.latN = TTin{icol.lat};
CTD.lonW = TTin{icol.lon};
CTD.Pressure = TTin{icol.pr};
CTD.Depth  = TTin{icol.de};
CTD.Temp = TTin{icol.te};
CTD.Salt = TTin{icol.sa};
% CTD.Conductivity = TTin{icol.co}; %no conductivity in these data KL 2/8/2024
CTD.O2 = TTin{icol.o2};
CTD.Beam = TTin{icol.beam};
CTD.PAR = TTin{icol.par};
CTD.Fluor = TTin{icol.fluor};

% Check for missing salts, o2 
indx = find(CTD.Salt < 0);
if ~isempty(indx)
    CTD.Salt(indx) = NaN;
end
indx = find(CTD.O2 < 0);
if ~isempty(indx)
    CTD.O2(indx) = NaN;
end
indx = find(CTD.Fluor < -990);
if ~isempty(indx)
    CTD.Fluor(indx) = NaN;
end

% Add derived variables
CTD.Fluor_filt = ZZ;
CTD.Fluor_offset = zeros(size(ZZ));
CTD.sig_theta = sw_pden(CTD.Salt,CTD.Temp,CTD.Pressure,0) - 1000;
CTD.AOU = aou(CTD.Salt,CTD.Temp,CTD.O2);
CTD.bvfrq = ZZ;
CTD.bvfilt = ZZ; %KL added 4/10/2025
CTD.vort = ZZ; %KL added 4/11/2025
CTD.Sunrise = ZZ;
CTD.Sunset = ZZ;
CTD.Season = ZZ;
CTD.VertZone = ZZ;
CTD.MLD_dens125 = ZZ;
CTD.MLD_bvfrq = ZZ;
CTD.MLD_densT2 = ZZ;
CTD.MLD_densT3 = ZZ; %KL adding MLD definitions, 6/24/2024
CTD.MLD_densT1 = ZZ;
CTD.MLD_dens2 = ZZ;
CTD.MLD_densGR = ZZ;
CTD.MLD_te2 = ZZ;
CTD.DCM = ZZ;
CTD.DCMde_top = ZZ;
CTD.DCMde_bot = ZZ;
CTD.DCMinML = ZZ;
CTD.par_est = ZZ;
CTD.par0 = ZZ;
CTD.kpar = ZZ;
CTD.z_par_1pcnt = ZZ;
CTD.z_par_halfpcnt = ZZ;
CTD.z_par_tenthpcnt = ZZ;
clear ZZ TTin dvec mtime

%%  create structure Xout to be saved as .mat file 
%leave here because Ruth later pulls values from here that get
%exported in the CSV file. The Xout structure was exported as *mat file in
%earlier versions of this script.
Xout = struct();  

castlist = unique(CTD.BATS_id);
ncast = length(castlist);
% row vectors
XX = ones(1,ncast) .* NaN;   
Xout.BATS_id = XX;
Xout.type = XX;
Xout.cruise = XX;
Xout.cast = XX;
Xout.mtime = XX;
Xout.year = XX;
Xout.month = XX;
Xout.day = XX;
Xout.hour = XX;
Xout.doy = XX;
Xout.lat = XX;
Xout.lon = XX;
Xout.MLD_dens125 = XX;
Xout.MLD_bvfrq = XX;
Xout.MLD_densT2 = XX;
Xout.MLD_densT3 = XX; %KL adding MLD 6/24/2024
Xout.MLD_densT1 = XX;
Xout.MLD_dens2 = XX;
Xout.MLD_densGR = XX;
Xout.MLD_te2 = XX;
Xout.DCM = XX;
Xout.DCMde_top = XX;
Xout.DCMde_bot = XX;
Xout.DCMinML = XX;
Xout.Season = XX;
Xout.Sunrise = XX;
Xout.Sunset = XX;
Xout.par0 = XX;
Xout.kpar = XX;
Xout.z_par_1pcnt = XX;
Xout.z_par_halfpcnt = XX;
Xout.z_par_tenthpcnt = XX;
Xout.fluor_offset = XX;
clear XX

%  rectangular arrays
XX= ones(MAXZ,ncast) .* NaN;
Xout.pr = XX;
Xout.de = XX;
Xout.te = XX;
% Xout.co = XX; %KL 2/8/2024 no conductivity data
Xout.sa = XX;
Xout.o2 = XX;
Xout.bac = XX;
Xout.fluor = XX;
Xout.fluor_filt = XX;
Xout.par = XX;
Xout.par_est = XX;
Xout.th = XX;
Xout.sig0 = XX;
Xout.rho = XX;
Xout.bvfrq = XX;
Xout.bvfilt = XX; %KL added 4/10/2025
Xout.vertZone = XX;
clear XX


% Copy data from CTD struct into Xout fields
for ii = 1:ncast
   theCast = castlist(ii);
   AA = theCast;
   BB = mod(AA, 1e7);
   CC = mod(BB, 1e3);

   Xout.BATS_id(ii) = AA;
   Xout.type(ii) = floor(AA/1e7);
   Xout.cruise(ii) = (BB - CC)/1e3;
   Xout.cast(ii) = CC;
   clear AA BB CC
  
   indx = find(CTD.BATS_id == theCast);
   itop = indx(1);
   nz = length(indx);
   if nz > MAXZ
       %KL updated error 1/21/2024 - this issue also occurs if you are
       %running the script on data that has already been processed
       error('foo:bar',['Note from Krista: You may have tried to run this script on data that has already been processed, \n', ...
           'or data arrays not long enough. In the later case, increase MAXZ to at least ',num2str(nz)])
   end
   Xout.mtime(ii) = decyear2dnum(CTD.decy(itop));
      dvec = datevec(Xout.mtime(ii));
   Xout.year(ii) = dvec(:,1);
   Xout.month(ii) = dvec(:,2);
   Xout.day(ii) = dvec(:,3);
   Xout.hour(ii) = dvec(:,4) + dvec (:,5) / 60;
   Xout.doy(ii) = datevec2doy(dvec)';
   Xout.lat(ii) = CTD.latN(itop);
   Xout.lon(ii) = CTD.lonW(itop) * -1;
%  Sunrise() uses lon, W is positive
     [rhr,rmin,shr,smin]=sunrise(dvec(2),dvec(3),dvec(1),Xout.lat(ii),-Xout.lon(ii));
   Xout.Sunrise(ii) = rhr * 100 + rmin;
   Xout.Sunset(ii) = shr * 100 + smin;
    clear dvec rhr shr rmin smin
   %  get profile data
   Xout.pr(1:nz,ii) = CTD.Pressure(indx);
   Xout.de(1:nz,ii) = CTD.Depth(indx);
   Xout.te(1:nz,ii) = CTD.Temp(indx);
%    Xout.co(1:nz,ii) = CTD.Conductivity(indx); %KL 2/8/2024 no conductivity data
   Xout.sa(1:nz,ii) = CTD.Salt(indx);
   Xout.o2(1:nz,ii) = CTD.O2(indx);
   Xout.bac(1:nz,ii) = CTD.Beam(indx);
   Xout.fluor(1:nz,ii) = CTD.Fluor(indx);
   Xout.par(1:nz,ii) = CTD.PAR(indx);
   
   % Check for missing salts, replace if they are near surface 
   ibad = find(isnan(Xout.sa(1:nz,ii)));
   if ~isempty(ibad)
       if max(ibad) < 5
           Xout.sa(ibad,ii) = Xout.sa(max(ibad)+1,ii);
           CTD.Salt(indx)= Xout.sa(1:nz,ii);
           %disp('Replaced missing salts at top of cast') %turn off warnings
       end
   end
    % compute derived variables
   Xout.th(:,ii) = sw_ptmp(Xout.sa(:,ii),Xout.te(:,ii),Xout.pr(:,ii),0);
   Xout.rho(:,ii) = sw_dens(Xout.sa(:,ii),Xout.te(:,ii),Xout.pr(:,ii));
   Xout.sig0(:,ii) = sw_pden(Xout.sa(:,ii),Xout.te(:,ii),Xout.pr(:,ii),0) - 1000;
   
   if ~isempty(ibad)
       CTD.sig_theta(indx) = Xout.sig0(1:nz,ii);
   end
   
   %change to also export potential vorticity (VORT)
   [BVFRQ,VORT,~] = sw_bfrq(Xout.sa(:,ii),Xout.te(:,ii),Xout.pr(:,ii),Xout.lat(ii));
   Xout.bvfrq(:,ii) = [NaN; BVFRQ];
   CTD.bvfrq(indx) = Xout.bvfrq(1:nz,ii);   % plug into CTD struct

   %KL 4/10/2025 addition: filter the BV data using Ruth's get_bvfilt.m
   Xout.bvfilt(:,ii) = get_bvfilt([NaN; BVFRQ],5);
   CTD.bvfilt(indx) = Xout.bvfilt(1:nz,ii);    

   %KL 4/11/2025 addition: calculation potential vorticity
   Xout.vort(:,ii) = [NaN; VORT];
   CTD.vort(indx) = Xout.vort(1:nz,ii);   % plug into CTD struct

   clear itop indx nz BVFRQ VORT
end %for ii

%%  loop through each cast again, add some more info, and write to individual files

flist = fieldnames(CTD);
nfields = length(flist);
for ii = 1:ncast
    theCast = castlist(ii);
    indx = find(CTD.BATS_id == theCast);
    nz = length(indx);
 %    
    XX = struct();
    for jj= 1:nfields
        fname = flist{jj};
        XX.(fname) = CTD.(fname)(indx);
    end
    clear jj
%    
    MLD = get_mld_ctd(XX);
    ML_ToUse = MLD.dens125;  % ensure same ML is used by various functions
    %KL adding options to XX,CTD,Xout; 6/24/2024 
    XX.MLD_dens125(:) = MLD.dens125;
    XX.MLD_bvfrq(:) = MLD.bvfrq;
    XX.MLD_densT2(:) = MLD.densT2;    
    XX.MLD_densT3(:) = MLD.densT3;    
    XX.MLD_densT1(:) = MLD.densT1;
    XX.MLD_dens2(:) = MLD.dens2;    
    XX.MLD_densGR(:) = MLD.densGR;    
    XX.MLD_te2(:) = MLD.te2;    
    
    CTD.MLD_dens125(indx) = MLD.dens125;
    CTD.MLD_bvfrq(indx) = MLD.bvfrq;
    CTD.MLD_densT2(indx) = MLD.densT2;
    CTD.MLD_densT3(indx) = MLD.densT3; 
    CTD.MLD_densT1(indx) = MLD.densT1; 
    CTD.MLD_dens2(:) = MLD.dens2;    
    CTD.MLD_densGR(:) = MLD.densGR;    
    CTD.MLD_te2(:) = MLD.te2;    
    
    Xout.MLD_dens125(ii) = MLD.dens125;
    Xout.MLD_bvfrq(ii) = MLD.bvfrq;
    Xout.MLD_densT2(ii) = MLD.densT2;    
    Xout.MLD_densT3(ii) = MLD.densT3;
    Xout.MLD_densT1(ii) = MLD.densT1;
    Xout.MLD_dens2(ii) = MLD.dens2;    
    Xout.MLD_densGR(ii) = MLD.densGR;    
    Xout.MLD_te2(ii) = MLD.te2;      
   
%   Sunrise was computed for Xout above
    XX.Sunrise(:) = Xout.Sunrise(ii);
    XX.Sunset(:) = Xout.Sunset(ii);
    CTD.Sunrise(indx) = XX.Sunrise(:);
    CTD.Sunset(indx) = XX.Sunset(:);
    
 % smooth the pr-binned fluorometer profile and identify CM layer.
 %  The bias is the average difference from zero between 400-600 meters
     filt_width = 3;
     [Xfilt,bias] = smooth_fluor(XX, filt_width);
if do_plots 
         figure; 
         subplot(1,2,1); hold on; axis ij; ylim([0 900]);
         plot(XX.sig_theta,XX.Pressure,'-','Linewidth',1);
         plot(xlim(),[MLD.densT2,MLD.densT2],'--y','Linewidth',1)
         plot(xlim(),[MLD.dens125,MLD.dens125],'--m','Linewidth',1)
         xlabel('Sigma-theta');
         ylabel('Pressure')
         subplot(1,2,2); hold on; axis ij; ylim([0 900]); xlim([-0.05 0.3]);
         title([num2str(XX.Cruise(1)),' - ',num2str(XX.Cast(1))]);
         plot(XX.Fluor,XX.Pressure,'-c','Linewidth',1.5);
         plot(Xfilt,XX.Pressure,'-k','Linewidth',1.5);
end
     if abs(bias) > 0.005 & giveNotice
       disp('Applying bias to fluor profile')
       Xfilt = Xfilt - bias;
if do_plots
    plot(Xfilt,XX.Pressure,'-r','Linewidth',1.5);
end
       XX.Fluor_offset(:) = bias;
       CTD.Fluor_offset(indx) = bias;
       Xout.fluor_offset(ii) = bias;
     end
    XX.Fluor_filt(:) = Xfilt;
    CTD.Fluor_filt(indx) = Xfilt;
    Xout.fluor_filt(1:nz,ii) = Xfilt;

%  
%  dcm layer
    percent = 0.33;
    izmax = find(XX.Depth > 400,1);
    if isempty(izmax)
        izmax = length(XX.Depth);
    end
    
    DCM = get_dcm_layer_ctd(Xfilt,XX.Depth, ML_ToUse,percent,izmax);   
    XX.DCM(:) = DCM.depth;
    CTD.DCM(indx) = DCM.depth;
    Xout.DCM(ii) = DCM.depth;
    %KL adding option to store DCM.itop and DCMinML because will need that to define
    %season 7/1/2024
    XX.DCMde_top(:) = DCM.de_top;
    CTD.DCMde_top(indx) = DCM.de_top;
    Xout.DCMde_top(ii) = DCM.de_top;

    XX.DCMde_bot(:) = DCM.de_bot;
    CTD.DCMde_bot(indx) = DCM.de_bot;
    Xout.DCMde_bot(ii) = DCM.de_bot;

    XX.DCMinML(:) = DCM.DCMinML;
    CTD.DCMinML(indx) = DCM.DCMinML;
    Xout.DCMinML(ii) = DCM.DCMinML;
       
    if ML_ToUse < -990;   % no ML defined -- likely a surface cast
        XX.DCM(:) = -999;
        CTD.DCM(indx) = -999;
        Xout.DCM(ii) = NaN;
        %KL adding these variables so they also end up not defined
        XX.DCMde_top(:) = -999;
        CTD.DCMde_top(indx) = -999;
        Xout.DCMde_top(ii) = NaN;

        XX.DCMde_bot(:) = -999;
        CTD.DCMde_bot(indx) = -999;
        Xout.DCMde_bot(ii) = NaN;

        XX.DCMinML(:) = -999;
        CTD.DCMinML(indx) = -999;
        Xout.DCMinML(ii) = NaN;
    end
    
if do_plots
    plot(xlim(),[DCM.depth,DCM.depth],'--g','Linewidth',1)    
    plot(xlim(),[DCM.de_top,DCM.de_top],'--k','Linewidth',1)    
    plot(xlim(),[DCM.de_bot,DCM.de_bot],'--k','Linewidth',1)  
    plot(xlim(),[ML_ToUse,ML_ToUse],'--m','Linewidth',1)
end
if isnan(DCM.itop) & giveNotice
       disp('top of DCM within ML')
end
    
    clear filt_width  Xfilt izmax 
%
%  fit PAR profile, get 1% 0.5% and 0.1% light levels
%  
   fitrange = [0 200];
   zmax = 300;
   PAR = get_BATS_par_vars(XX.PAR,XX.Depth,fitrange,zmax);
   
%    figure; hold on; axis ij;
%    parin = XX.PAR;
%    parin(parin < -900) = NaN;
%    din = XX.Depth;
%    din(din < -900) = NaN;
%    plot(parin,din,'-k','Linewidth',1.5)
%    plot(PAR.par_est,din,'-r','Linewidth',1.5);
%    ylim([0 300]);
%    plot(1,PAR.z_par_1pcnt,'*');
%    plot(1,PAR.z_par_halfpcnt,'*');
%    plot(1,PAR.z_par_tenthpcnt,'*');
%

clear din parin

%  save values     
   Pfields = fieldnames(PAR);
  
   for kk=1:length(Pfields)
       Pname = Pfields{kk};
       AA = PAR.(Pname);
       if strcmp(Pname,'par_est')
           Xout.(Pname)(1:nz,ii) = AA;
           AA(isnan(AA)) = -999;
           XX.(Pname)(:) = AA;
           CTD.(Pname)(indx) = AA;
       else
            Xout.(Pname)(ii) = AA;
           if isnan(AA)
               AA = -999;
           end
           BB = ones(size(XX.(Pname))) .* AA;
           XX.(Pname)(:) = BB;
           CTD.(Pname)(indx) = BB;
       end
   end
   
   clear AA BB zmax fitrange PAR
   
% Label Vertical Zones
   XX.VertZone(:) = label_vertical_zones_ctd(XX,DCM,ML_ToUse);
%    if any(XX.VertZone == -999)
%        blah = label_vertical_zones_ctd(XX,DCM,ML_ToUse);
%    end
   CTD.VertZone(indx) = XX.VertZone(:);
   Xout.vertZone(1:length(indx),ii) = XX.VertZone(:);
      ibad = find(Xout.vertZone(:,ii) < -990); %nan any missing values
   Xout.vertZone(ibad,ii) = NaN;  
   
% Label Seasons
%KL April 2025: use Excel file from July 2024 with season information
%send it into this function so I am not opening file for each CTD cast
   %trans_dates = get_season_dates(Xout.year); %new function from KL
   theCode = label_seasons_ctd(XX,trans_dates); %send in trans_dates
   disp([num2str(Xout.year(ii)),' ',num2str(Xout.month(ii)),' ',num2str(Xout.day(ii)),'  Season: ', num2str(theCode)]);
   XX.Season(:) = theCode;
   CTD.Season(indx) = theCode;
   Xout.Season(ii) = theCode;
   clear theCode;
% Change NaNs to -999 for output as csv
    for jj= 1:nfields
        fname = flist{jj};
        BB = XX.(fname);
        BB(isnan(BB)) = -999.;
        XX.(fname) = BB;       
    end

% Output the cast %KL note: for the LTT project, need the individual casts
%KL note - these are square matrices, to whatever depth is set in MAXZ  
    %%slight change, this is what Ruth had:
    % fmt = '%8d_%1d%04d_%03d_ctd.csv';
    % outfile = sprintf(fmt,XX.yyyymmdd(1),Xout.type(ii),XX.Cruise(1),XX.Cast(1));
    % TTcast = struct2table(XX);
    % writetable(TTcast,strcat(outdir,filesep,outfile));
    %I think it will be easier to save the table itself
    fmt = '%8d_%1d%04d_%03d_ctd.mat';
    outfile = sprintf(fmt,XX.yyyymmdd(1),Xout.type(ii),XX.Cruise(1),XX.Cast(1));
    TTcast = struct2table(XX);
    save(strcat(outdir,outfile),'TTcast');
    

    clear XX DCM MLD TTcast
  
end % for ii

   if 0 
       %don't think I need the CSV files for the LTT project
       fmt = 'CRU_%1d%04d_ctd.csv';
       outfile = sprintf(fmt,Xout.type(1),CTD.Cruise(1));
       %stick in the directory I want
       outfile = strcat(outdir,outfile);
       
       % replace any NaNs in CTD struct with -999.
       flist = fieldnames(CTD);
       nfields = length(flist);
       for jj=1:nfields
           CC=CTD.(fname);
           CC(isnan(CC)) = -999;
           CTD.(fname) = CC;
       end
       TTcruise = struct2table(CTD);
       writetable(TTcruise,outfile);
   end
   
    if do_plots
        tilefigs
        reply = input(' Press ENTER to close figures....');
        close all
    end

end %function
