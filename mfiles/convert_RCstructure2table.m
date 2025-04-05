function tableOut = convert_RCstructure2table(structureIn,trimToShort)
%function tableOut = convert_RCstructure2table(structureIn,trimToShort)
%The structures holding the CTD have a mix of sizes, so I need a custom
%function to make them into tables (which will be easier to navigate).
%MATLAB will not transpose items within the structure 
%Krista Longnecker, 21 June 2024
%
% structureIn is one of the CTD structures
% trimToShort is 1/0 if I only want to keep the items with one value per
% cruise/cast (easier to put that in here because of the way MATLAB handles
% tables)
% tableOut is the resulting table, will always be one row
tableOut = table();

fn = fieldnames(structureIn);

for a = 1:length(fn)
    tableOut(1,fn(a)) = {structureIn.(fn{a})'};    
end
clear a

%will always have two values - one small number (equal to number of casts)
%and one bigger number  that is the length of the concatenatedata
ni = structfun(@length,structureIn);
uni = unique(ni);
k = find(ni==min(uni));
tableOut = tableOut(:,k);

