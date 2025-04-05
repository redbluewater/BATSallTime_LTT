function [ mtime ] = decyear2dnum( dec_year )
% function [ mtime ] = decyear2dnum( dec_year )
%  convert decimal year to matlab time

year=floor(dec_year);
  partialYear = mod(dec_year,1);
date0 = datenum(num2str(year),'yyyy');
date1 = datenum(num2str(year+1),'yyyy');
daysInYear = date1 - date0;
mtime=date0 + partialYear .* daysInYear;

end
