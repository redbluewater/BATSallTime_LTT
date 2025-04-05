function [theCode] = label_seasons_ctd(CTD,trans_dates)
% function [theCode] = add_seasons(CTD,trans_dates)
%   Assigns a season code based dates listed in trans_dates
%INPUT
% CTD : struct containing a single BIOSSCOPE CTD cast 
% trans_dates can be [] or a struct with fields:
%    .mixed  (n x 2 array of start/end dates)  
%    .spring 
%    .strat
%    .fall
%OUTPUT
% theCode:
%   1 :  Mixed  begins when top of CM layer no longer defined
%   2 :  Spring begins first day MLD shoals above top of CM layer
%   3 :  Strat  begins when MLD is consistently above top of CM
%   4 :  Fall   begins when MLD first goes below top of CM 
%
% Ruth Curry, BIOS / ASU
% Uploaded for BIOS-SCOPE project 19 October 2023
% KL edited, now entirely relying on dates in pre-determined Excel file
% Krista Longnecker; 3 July 2024

mtime = decyear2dnum(CTD.decy(1));

%only use pre-defined dates to set the season
if isstruct(trans_dates)  % check the date of cast against the dates provided
    
    TD=trans_dates.mixed;
    [ndates,~] = size(TD);
    for ii=1:ndates
        if mtime >= TD(ii,1) && mtime < TD(ii,2)
            theCode = 1;
            return
        end
    end
    
    TD=trans_dates.spring;
    [ndates,~] = size(TD);
    for ii=1:ndates
        if mtime >= TD(ii,1) && mtime < TD(ii,2)
            theCode = 2;
            return
        end
    end
    
    TD=trans_dates.strat;
    [ndates,~] = size(TD);
    for ii=1:ndates
        if mtime >= TD(ii,1) && mtime < TD(ii,2)
            theCode = 3;
            return
        end
    end
    TD=trans_dates.fall;
    [ndates,~] = size(TD);
    for ii=1:ndates
        if mtime >= TD(ii,1) && mtime < TD(ii,2)
            theCode = 4;
            return
        end
    end
end %if (trans_dates)

%  if reach here, stop the code as I have no date information and I don't
%  want to add nominal dates without stopping to think about what happened.
% This is likely an error in the Excel file used to set the boundaries.
keyboard
error('There is no pre-defined date, consider what happened\n')

end

