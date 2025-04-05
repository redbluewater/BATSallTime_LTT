function season_dates = get_season_dates(year)
%function season_dates = get_season_dates(year)
% Make a separate function to set up the seasons - I need to call this each
% time I open a BATS data file - the year parameter is set dynamically from
% there
% Krista Longnecker, 8 February 2024

% If glider not available can use general dates:  15-Dec: 01-Apr : 20-Apr : 01-Nov
%but still need things as MATLAB datenum format
season_dates.mixed = [datenum(year,12,15), datenum(year,4, 1)];
season_dates.spring = [datenum(year,4,1), datenum(year,4, 20)];
season_dates.strat = [datenum(year,4,20), datenum(year,11, 1)];
season_dates.fall = [datenum(year,11,1), datenum(year,12, 15)];


% %this might be useful later so I don't have to use datenum
% for a = 1:size(season_dates.mixed,1)
%     season_dates.year(a,1) = year(datetime(datestr(season_dates.mixed(a,1))));
%     dt.mixed(a,1) = datetime(datestr(season_dates.mixed(a,1)));
%     dt.mixed(a,2) = datetime(datestr(season_dates.mixed(a,2)));
%     
%     dt.spring(a,1) = datetime(datestr(season_dates.spring(a,1)));
%     dt.spring(a,2) = datetime(datestr(season_dates.spring(a,2)));
% 
%     dt.strat(a,1) = datetime(datestr(season_dates.strat(a,1)));
%     dt.strat(a,2) = datetime(datestr(season_dates.strat(a,2)));
%     
%     dt.fall(a,1) = datetime(datestr(season_dates.fall(a,1)));
%     dt.fall(a,2) = datetime(datestr(season_dates.fall(a,2)));
% end
    

