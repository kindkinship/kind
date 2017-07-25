function estimateParamsClusters(sampleIds, afClusters, genotypes, updatedBPs, cmbrship, samplePairsKnown, coeffKnown)

out = repmat(NaN,max(cmbrship),1);

for c=unique(cmbrship)'
    idxes = find(c == cmbrship);
    expectedDist = calculateExpectedDistance(bps, afClusters(c).af);
    p = estimateParam(sampleIds, genotypes, bps, samplePairsKnown, coeffKnown, expectedDist);
end;    