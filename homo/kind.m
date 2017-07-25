% kind.m

clear;

% individual pairs of known relationship and their kinship coefficients
knownrel = 'known.txt';
mapfile = 'example.map';
pedfile = 'example.ped';
frqfile = 'example.frq';

if 0
    [map,ped,frq,samplePairsKnown,coeffKnown] = readInputFiles(mapfile, pedfile, frqfile, knownrel);
    save example.mat map ped frq samplePairsKnown coeffKnown;
else
    load example.mat;
end;

[filteredped, filteredmap, filteredfrq, validIdxesPed, validIdxesMap, validIdxesFrq] = filtering(ped, map, frq);
[sampleIds, genotypes] = getGenotypes(filteredped);
updatedBPs = updateBPs(filteredmap);
expDist = calculateExpectedDistance(updatedBPs, filteredfrq);

p = estimateParam(sampleIds, genotypes, updatedBPs, samplePairsKnown, coeffKnown, expDist);
coeffest = estimateKinshipCoefficients(sampleIds, genotypes, updatedBPs, expDist, p);
outputToFile(coeffest, 'estimatedcoeff.txt');