function [sampleIds, genotypes] = getGenotypes(d)

sampleIds = d(:,1);
genotypes = repmat(0,size(d,1), size(d,2)-6);

for i=1:size(d,1)
    idxes = find(strcmp('1',d(i,7:end)));
    genotypes(i,idxes) = 1;
end;    