function out = ce(in,type);
%function out = ce(in,type);
% Check if a value is empty. 
% If yes, enter NaN or NaT; 
% If no, enter that value; with the default as NaN
% Krista Longnecker, 3 July 2024

if nargin == 1
    type = 'n';
end

if ~isempty(in)
    %just use as is
    out = in;
elseif isempty(in) && isequal(type,'n')
    out = NaN;
elseif isempty(in) && isequal(type,'t')
    out = NaT;
end
    