function out = calculateAFclusters(cmbrship, genotypes)

out(1).af = [];

for c=min(cmbrship):max(cmbrship)
    cnt_ma = [];
    idxes = find(c == cmbrship);
    for k=idxes'
        if isempty(cnt_ma)
            cnt_ma = sum(reshape(genotypes(k,:),2,size(genotypes,2)/2));
        else
            cnt_ma = cnt_ma + sum(reshape(genotypes(k,:),2,size(genotypes,2)/2));
        end;
    end;    
    
    out(c).af = cnt_ma'/(2*length(idxes));
end;    