function cmbrship = clustering(genotypes, minK, maxK)

alleleCount = [];
for i=1:size(genotypes,1)
    alleleCount = [alleleCount;sum(reshape(genotypes(i,:),2,size(genotypes,2)/2))];
end;

clusterQuality = repmat(NaN,maxK,1);
cmbr = repmat(NaN,size(genotypes,1),maxK);
for k = minK:maxK
    [IDX, centroid, sumd] = kmeans(alleleCount, k, 'Distance', 'correlation', 'Replicates', 20, 'EmptyAction', 'drop');     
    clusterQuality(k,1) = sum(sumd);
    cmbr(1:numel(IDX),k) = IDX;
end;

[C, numClusters] = max(clusterQuality);
[cmbrship, centroid, sumd] = kmeans(alleleCount, numClusters, 'Distance', 'correlation', 'Replicates', 10, 'EmptyAction', 'drop');
