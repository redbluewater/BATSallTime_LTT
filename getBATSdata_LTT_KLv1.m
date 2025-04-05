%Set up to calculate Ruth's derived variables from all BATS cruises (e.g.,
%MLD, season, and VZ); starting with BATS cruise #1
%based on Ruth's prior code : create_biosscope_files_2022_2023.m
% Original code from Ruth Curry, BIOS / ASU
% Krista Longnecker; 8 February 2024
% Krista Longnecker; 21 June 2024
%Modify existing MATLAB code to get the BV data from all time at BATS;
%using MATLAB because Ruth made these calculations in MATLAB
%Krista Longnecker; 27 June 2024
%Krista Longnecker ; 4 April 2025 (turn off Dropbox syncing while running
%this : will massively speed up the code)
clear all 
close all

%%add options depending on computer, KL is jumping between computers
if isequal(getenv('COMPUTERNAME'),'ESPRESSO')
    %% add ./BIOSSCOPE/CTD_BOTTLE/mfiles into matlab path
    addpath(genpath('C:\Users\klongnecker\Documents\GitHub_espresso\BATSallTime_LTT\mfiles'));    
    
    %% update the folder information before getting started
    rootdir = 'C:\Users\klongnecker\Documents\Dropbox\Current projects\Kuj_BIOSSCOPE\RawData\';
    %Krista has put the next two folders outside the space accessible by GitHub
    %These files are too large to put into GitHub
    workdir = fullfile(rootdir,'RCcalcBATS\data_temporary\');
    % workdir = fullfile(rootdir,'RCcalcBATS\data_copySmall_testing\');
    outdir = fullfile(rootdir,'RCcalcBATS\data_holdingZone\');
elseif isequal(getenv('COMPUTERNAME'),'LONGNECKER-1650')
    %% add ./BIOSSCOPE/CTD_BOTTLE/mfiles into matlab path    
    addpath(genpath('C:\Users\klongnecker\Documents\Dropbox\GitHub\BATSallTime_LTT\mfiles'));

    %% update the folder information before getting started
    rootdir = 'C:\Users\klongnecker\Documents\Dropbox\Current projects\Kuj_BIOSSCOPE\RawData\';
    %Krista has put the next two folders outside the space accessible by GitHub
    %These files are too large to put into GitHub
    workdir = fullfile(rootdir,'RCcalcBATS\data_temporary\');
    % workdir = fullfile(rootdir,'RCcalcBATS\data_copySmall_testing\'); up to 1992
    outdir = fullfile(rootdir,'RCcalcBATS\data_holdingZone\');
end

gitdir = pwd;
%set this to one if you want to see a plot for each cast (that will clearly
%be a ton of plots, so this is best used if you to look at a preset number
%of casts within the full set of casts)
do_plots = 0;
   
NameOfFile = 'BATSdataForBVplots.2025.04.04.mat';

%use the Excel file from July 2024 to set the seasons - do here and then
%send result (trans_dates) into calculate_BATSderivedVariables_forLTT
%Use this function to make a MATLAB structure with transition dates
seasonsFile = fullfile(gitdir,'BATS_seasons_wKLedits.2024.07.05.xlsx');
%use this function to reformat the dates, set fName in calcDerivedVariables
trans_dates = reformat_season_dates(seasonsFile) ; 

%now do the calculations

%Ruth set this up for txt files, but the BATS txt files are a pain, use
%the BATS *mat files
cd(workdir)
dirlist = dir('*.mat');
%delete some names...not the best way to do this, but will work
s = contains({dirlist.name},'YR');
ks = find(s==1);
dirlist(ks) = []; clear s ks
s = contains({dirlist.name},'bats_ctd.mat');
ks = find(s==1);
dirlist(ks) = []; clear s ks
s = contains({dirlist.name},'bval_ctd.mat');
ks = find(s==1);
dirlist(ks) = []; clear s ks
s = contains({dirlist.name},'working.mat');
ks = find(s==1);
dirlist(ks) = []; clear s ks

nfiles = length(dirlist);
stepOne = table();

% doFiles = 1; %use for testing, smaller number of files
% doFiles = 3; %use for testing, more files in case you need to test the loop
doFiles = nfiles; %do  everything
for ii = 1:doFiles;
   fname = dirlist(ii).name;      
   infile = fullfile(workdir,fname);
   
   %use modified function from KL, 4/4/2025
   CTD = calculate_BATSderivedVariables_forLTT(infile,trans_dates,do_plots,outdir,0);
   %make a table, easier to manipulate; have to write a custom script to
   %make a table given the complexity of these structures
   trim = 1; %set this to one to only keep the values that are one per cruise/cast
   T = convert_RCstructure2table(CTD,trim); %this is a new KL function 6/21/2024
   stepOne = [stepOne;T];
   
   clear idx T trim CTD infile fname
end
clear ii nfiles doFiles dirlist do_plots

%Can stop here first....next up will be to concatenate all the CSV files
%sitting in the outdir.
%go to the folder with the MAT files, read them in and concatenate into one
%matrix...will be pretty big
cd(outdir)
dirlist = dir('*.mat');

%use the first file to set the stage
load(dirlist(1).name)
%Only keeping the upper water column bc that is the focus of this paper.
%Delete depths more than 220m.
k = find(TTcast.Depth<225);
allData = TTcast(k,:);
clear TTcast
for a = 2:length(dirlist)
    load(dirlist(a).name);
    k = find(TTcast.Depth<225); %only keep upper water column
    allData = [allData ; TTcast(k,:)];
    clear TTcast
end
clear a dirlist


%now have many rows of data, save this (after some housecleaning)
cd(gitdir)

%export this as a CSV file...that will be better to share then the MATLAB
%file; file is big, do not put in gitdir (or setup .gitignore file)
writetable(allData,strcat(outdir,'BATS_CTDdata_2db.2025.04.04.csv'))

clear outdir rootdir workdir gitdir stepOne 

save(NameOfFile)


