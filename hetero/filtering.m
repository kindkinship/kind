function [outped, outmap, outfrq, validIdxesPed, validIdxesMap, validIdxesFrq] = filtering(ped, map, frq)

outped = {};
outmap = {};
outfrq = [];
validIdxesPed = [];
validIdxesMap = [];
validIdxesFrq = [];

validIdxesFrq = find(frq > 0);
outfrq = frq(validIdxesFrq);
map = map(validIdxesFrq,:);
idxes = 2*validIdxesFrq;
idxes = 6+sort([idxes;idxes-1],'ascend');
ped = [ped(:,1:6) ped(:,idxes)];

validIdxesPed = [1:6];
validIdxesMap = [];
for i=7:2:size(ped,2)
    missingIdxes1 = find(strcmp('0', ped(:,i)));
    missingIdxes2 = find(strcmp('0', ped(:,i+1)));
        
    if isempty(missingIdxes1) && isempty(missingIdxes2)
        validIdxesPed = [validIdxesPed i i+1];
        validIdxesMap = [validIdxesMap;(i-5)/2];
    end;    
end;
outped = ped(:,validIdxesPed);
outmap = map(validIdxesMap,:);
