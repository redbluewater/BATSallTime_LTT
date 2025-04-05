function cruise5 = id2cruise(dataIn);
%getting tired of pulling out the five digit cruise ID from a full BATS id
% input:
%if dataIn is multiple rows, return multipe rows
%if dataIn is one value, return one value
% 1 July 2024, Krista Longnecker

if size(dataIn,1) > 1
    for a = 1:size(dataIn,1)
        one = char(string(dataIn(a)));
        one = one(1:5);
        cruise5(a,1) = str2double(one);
        clear one
    end
    clear a
    
else 
    one = char(string(dataIn));
    one = one(1:5);
    cruise5 = str2double(one)
end
	