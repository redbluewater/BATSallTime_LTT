%Compare KL calculation for BV to Craig's calculation
%KL 8 April 2025
%run this using MATLAB 2023a (different versions are reading in the Excel
%file differently
clear all
close all

load('BATSdataForBVplots.2025.04.04.mat');
%need cruise to compare with Craig's calculation
for a = 1:size(allData,1)
    if rem(a,50000)==1 
        fprintf('iteration %d for %d peaks\n',a,size(allData,1)) 
    end

    one = char(string(allData.BATS_id(a)));
    cruise5(a,1) = {one(1:5)}; %do this, rewriting the table is slow!
    clear one
end
clear a
allData.depth_rounded = round(allData.Depth);

k = find(allData.bvfrq==-999);
allData.bvfrq(k) = NaN;
clear k

k = find(allData.bvfilt==-999);
allData.bvfilt(k) = NaN;
clear k

T = readtable('Brunt-V_is_l__Frequency_Data_2.xlsx'); %this is Craig's file

small = 10^-100;
k = find(T.N_squared<0);
T.N_squared(k) = small;
clear k

k = find(allData.bvfrq<0);
allData.bvfrq(k) = small;
clear k
k = find(allData.bvfilt<0);
allData.bvfilt(k) = small;
clear k

for a = 1:size(T,1)
    if rem(a,1000)==1 
        fprintf('iteration %d for %d peaks\n',a,size(T,1)) 
    end
   one = T.cruise_depth{a};
   r = regexp(one,'_');
   T.cruise{a} = one(1 : r-1);

   %dtm for Craig as well (follow Ruth's code)
   T.BV_CC(a) = power(nansum(log10(T.N_squared(a)))./(1),10);
   clear one r

   s = strcmp(T.cruise(a),cruise5);
   ks = find(s==1 & T.Depth_m_(a)==allData.depth_rounded);

   if ~isempty(ks)
        %From Ruth (4/9/2025, change to this);
        T.BV_KL_count(a) = length(ks);
        %T.BV_KL(a) = power(nansum(log10(allData.bvfrq(ks)))./length(ks),10);
        T.BV_KL(a) = power(nansum(log10(allData.bvfilt(ks)))./length(ks),10);
   end
   clear s ks
end
clear a

% plot(T.BV_KL,T.N_squared,'.')
% gscatter(T.BV_KL, T.BV_CC, T.BV_KL_count,[],[],20)
plot(T.BV_KL, T.BV_CC,'.')
hold on
% YL = ylim;
% line(xlim,xlim,'color','k')
% ylim(YL)
xlabel('BV from KL, filtered with Ruth''s code')
ylabel('BV from Craig')
title('Values at 10^2^0 are negative BV (unstable, converted to fixed value)')
set(gca,'yscale','log')
set(gca,'xscale','log')

saveas(gcf,'BATS_bvfrq_cfKLandCC.jpg')


   % 
   % if ~isempty(ks)
   %     if isequal(length(ks),1); 
   %      T.BV_KL(a) = allData.bvfrq(ks); 
   %      T.BV_KL_count(a) = 1;
   %     elseif length(ks)>1
   %         %take the mean? From Ruth (4/9/2025, change to this);
   %          T.BV_KL_count(a) = length(ks);
   %          T.BV_KL(a) = nansum(log10(allData.bvfrq(ks)))
   %          T.BV_KL(a) = power(nansum(log10(allData.bvfrq(ks)))./length(ks),10);
   % 
   % 
   %         % ./length(ks),20;
   % 
   %     end
