function season_dates = reformat_season_dates(fName)
%reformat the Excel file with dates into a structure matching Ruth's format
%Krista Longnecker 3 July 2024

t = readtable(fName);

%make into the format expected by Ruth's code (datenums in a structure)
season_dates = struct();
season_dates.mixed = datenum(t{:,{'mixed_begin','mixed_end'}});
season_dates.spring = datenum(t{:,{'spring_begin','spring_end'}});
season_dates.strat = datenum(t{:,{'strat_begin','strat_end'}});
season_dates.fall = datenum(t{:,{'fall_begin','fall_end'}});
season_dates.year = t.year;

%these are datenum - that is what Ruth uses - convert to datetime:
%datetime(season_dates.mixed(1),'ConvertFrom','datenum');