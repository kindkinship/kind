% kind_hetero.m

clear;

% individual pairs of known relationship and their kinship coefficients
knownrel = 'known.txt';
mapfile = 'example.map';
pedfile = 'example.ped';
frqfile = 'example.frq';

if 0 
disp('readInputFiles ...');
[map,ped,frq,samplePairsKnown,coeffKnown] = readInputFiles(mapfile, pedfile, frqfile, knownrel);
%[map,ped,frq] = readInputFiles(mapfile, pedfile, frqfile);
save example.mat map ped frq samplePairsKnown coeffKnown;
%save example.mat map ped frq;
else
    load example.mat;
end;

disp('filtering ...');
[filteredped, filteredmap, filteredfrq, validIdxesPed, validIdxesMap, validIdxesFrq] = filtering(ped, map, frq);
disp('getGenotypes ...');
[sampleIds, genotypes] = getGenotypes(filteredped);
disp('updateBPs ...');
updatedBPs = updateBPs(filteredmap);
disp('clustering ...');
cmbrship = clustering(genotypes, 2, 5);

disp('calculateAFclusters ...');
afClusters = calculateAFclusters(cmbrship, genotypes);

disp('estimateKinshipCoefficientsWBC ...');
[coeffest, clustersWOknownrel] = estimateKinshipCoefficientsWBC(sampleIds, genotypes, updatedBPs, cmbrship, samplePairsKnown, coeffKnown, afClusters);
outputToFile(coeffest, 'estimatedcoeff.txt');
