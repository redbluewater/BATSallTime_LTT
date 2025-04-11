%plot up the BV data...down to 200 m only
%starting in MATLAB bc that is where the data are
%KL 4 April 2025
clear all
close all
%KL compilation from the BATS 2db CTD casts
load BATSdataForBVplots.2025.04.04.mat

%Craig's calculation from the synoptic cast
T = readtable('Brunt-V_is_l__Frequency_Data_2.xlsx');
%set rows from Craig's data where BV ==0 to NaN (calculation issue)
k = find(T.N_squared==0);
T.N_squared(k) = NaN;
clear k

%KL 4/9/2025
%following Ruth - change negative values to a small number and then plot as
%log10
small = 10^-100;
k = find(T.N_squared<0);
T.N_squared(k) = small;
k = find(allData.bvfrq<0);
allData.bvfrq(k) = small;

k = find(allData.bvfilt<0);
allData.bvfilt(k) = small;


figure
%%% plot the filtered BV data (add 4/10/2025
subplot(211)

x = allData.decy;
y = allData.Depth; 
% z = allData.bvfrq;
z = power(log10(allData.bvfilt),10);

k = find(y==-999);
y(k) = [];
x(k) =[];
z(k)= [];

k = find(z==-999);
y(k) = [];
x(k) =[];
z(k)= [];

% forMax = max(z);
forMax = 1e7;
useC = [min(z) forMax]; 
doPlot(x,y,z,useC)

xlabel('time')
ylabel('depth (m)')
title('BV frequency calculated from all BATS data using 2db resolution CTD cast data, FILTERED')
h = colorbar();


%%plot up Craig's data
%%plot up Craig's data
%%plot up Craig's data
%%plot up Craig's data
subplot(212)

x = T.decy;
y = T.Depth_m_; 
% z = T.N_squared;
z = power(log10(T.N_squared),10);

k = find(y==-999);
y(k) = [];
x(k) =[];
z(k)= [];

k = find(z==-999);
y(k) = [];
x(k) =[];
z(k)= [];

forMax = max(z);
useC = [min(z) forMax];     
doPlot(x,y,z,useC)

xlabel('time')
ylabel('depth (m)')
title('Carlson BV freq')
h = colorbar()

title_up('Colors are BV as ''power(log10(BV),10)''')


if 0
    %saveas(gcf,char(strcat(nDir,filesep,mtabNames(a),'.jpg')),'jpg')
    %if I don't do the first bit, I don't get vectors for Illustrator
    set(gcf,'paperpositionmode','auto')
    set(gcf,'renderer','Painters')
    print(gcf,'BATS_bv_allYears.svg','-dpdf')   
elseif 1 
    saveas(gcf,'BATS_bvfrq_allYears.jpg')
end

% figure
% subplot(311)
% 
% x = allData.decy;
% y = allData.Depth; 
% % z = allData.bvfrq;
% z = power(log10(allData.bvfrq),10);
% 
% k = find(y==-999);
% y(k) = [];
% x(k) =[];
% z(k)= [];
% 
% k = find(z==-999);
% y(k) = [];
% x(k) =[];
% z(k)= [];
% 
% forMax = max(z);
% useC = [min(z) forMax]; 
% doPlot(x,y,z,useC)
% 
% xlabel('time')
% ylabel('depth (m)')
% title('BV frequency calculated from all BATS data using 2db resolution CTD cast data')
% h = colorbar();
% 
% 
% %%% plot the filtered BV data (add 4/10/2025
% subplot(312)
% 
% x = allData.decy;
% y = allData.Depth; 
% % z = allData.bvfrq;
% z = power(log10(allData.bvfilt),10);
% 
% k = find(y==-999);
% y(k) = [];
% x(k) =[];
% z(k)= [];
% 
% k = find(z==-999);
% y(k) = [];
% x(k) =[];
% z(k)= [];
% 
% forMax = max(z);
% useC = [min(z) forMax]; 
% doPlot(x,y,z,useC)
% 
% xlabel('time')
% ylabel('depth (m)')
% title('BV frequency calculated from all BATS data using 2db resolution CTD cast data, FILTERED')
% h = colorbar();
% 
% 
% %%plot up Craig's data
% %%plot up Craig's data
% %%plot up Craig's data
% %%plot up Craig's data
% subplot(313)
% 
% x = T.decy;
% y = T.Depth_m_; 
% % z = T.N_squared;
% z = power(log10(T.N_squared),10);
% 
% k = find(y==-999);
% y(k) = [];
% x(k) =[];
% z(k)= [];
% 
% k = find(z==-999);
% y(k) = [];
% x(k) =[];
% z(k)= [];
% 
% forMax = max(z);
% useC = [min(z) forMax];     
% doPlot(x,y,z,useC)
% 
% xlabel('time')
% ylabel('depth (m)')
% title('Carlson BV freq')
% h = colorbar()
% 
% title_up('Colors are BV as ''power(log10(BV),10)''')


if 0
    %saveas(gcf,char(strcat(nDir,filesep,mtabNames(a),'.jpg')),'jpg')
    %if I don't do the first bit, I don't get vectors for Illustrator
    set(gcf,'paperpositionmode','auto')
    set(gcf,'renderer','Painters')
    print(gcf,'BATS_bv_allYears.svg','-dpdf')   
elseif 1 
    saveas(gcf,'BATS_bvfrq_allYears.jpg')
end



function doPlot(x,y,z,useC)
    sp = 100;
    x1=linspace(min(x),max(x),sp);
    y1=linspace(min(y),max(y),sp);
    [Xi,Yi]=meshgrid(x1,y1);
    %Zi=griddata(x,y,z,Xi,Yi,'linear');
    SN=100;
%     LX = 2;
%     LY = 2;
    % xq(1)=0.5
    % yq(1)=0.5
    %function[zq,eq]=divagrid(x,y,z,xq,yq,SN,LX,LY,masking)
    %using defaults for now
    %Zi=divagrid(x,y,z,Xi,Yi); %DIVA, replaces griddata, slow, but OK
    %let divagrid chose the length scales (LX and LY), also what ODV does 
    Zi=divagrid(x,y,z,Xi,Yi,SN); %DIVA, replaces griddata
%     Zi=divagrid(x,y,z,Xi,Yi,SN,LX,LY); %DIVA, replaces griddata

    pcolor(x1,y1,Zi);
    colormap(parula)

    shading flat

    caxis([useC])
    set(gca,'ydir','reverse')
%     datetick('x','mmmyy','keepticks')

%     colorbar
    %hold on
    %plot(x,y,'.','color',0.5*ones(1,3)) %turn off the dots
    set(gcf,'position',[-1070 33 844 1081])
    
end

function doContour(x,y,z,useC, useV)
    x1=linspace(min(x),max(x),50);
    y1=linspace(min(y),max(y),50);
    [Xi,Yi]=meshgrid(x1,y1);
    Zi=griddata(x,y,z,Xi,Yi);
    if nargin==5
        [c h]=contour(x1,y1,Zi,useV,'k');
    elseif nargin ==4
        [c h]=contour(x1,y1,Zi,'k');
    elseif nargin==3
        [c h]=contour(x1,y1,Zi,'k','linewidth',2);
        ra
        %pcolor(x1,y1,Zi)
        %shading flat
    else       
        error('Wrong number of arguments')
    end

    clabel(c,h)
%     caxis(useC) %otherwise the contour map will reset the colormap
end

