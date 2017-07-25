function bps = updateBPa(map)

bps = [];

uniqueChrnos = unique(map(:,1));
for chrno = uniqueChrnos'
    idxes = find(strcmpi(chrno, map(:,1)));
    if isempty(bps)
        bps = str2double(map(idxes,4));
    else
        tmp = str2double(map(idxes,4));
        bps = [bps;tmp-tmp(1)+1+bps(end)];
    end;
end;
