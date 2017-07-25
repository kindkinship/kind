function [coeffest, clustersWOknownrel] = estimateKinshipCoefficientsWBC(sampleIds, genotypes, bps, cmbrship, samplePairsKnown, coeffKnown, afClusters)

coeffest = {};

% within-cluster estimation
for c=unique(cmbrship)'
    clustersWOknownrel(c).samples = {}; % because there is no sample pairs with known relationship for the samples in the clusters
    
    idxes = find(cmbrship == c);
    sIds = sampleIds(idxes);
    gts = genotypes(idxes,:);
    
    trainSamplePairs = {};
    trainCoeffs = [];    
    for i=1:size(samplePairsKnown,1)
        if isempty(find(strcmpi(samplePairsKnown{i,1},sIds)))
            continue; 
        end;
        
        if isempty(find(strcmpi(samplePairsKnown{i,2},sIds)))
            continue;
        end;
        
        trainSamplePairs = [trainSamplePairs;samplePairsKnown(i,1:2)];
        trainCoeffs = [trainCoeffs;coeffKnown(i)];
    end;    
    
    if isempty(trainSamplePairs)
        clustersWOknownrel(c).samples = [clustersWOknownrel(c).samples;sIds];
        continue;
    end;
    
    indices = find(afClusters(c).af > 0);
    expectedDist = calculateExpectedDistance(bps(indices), afClusters(c).af(indices));
    
    genotypeIndices = sort([2*indices-1 2*indices],'ascend');
    p = estimateParam(sampleIds, genotypes(:,genotypeIndices), bps(indices), trainSamplePairs, trainCoeffs, expectedDist);    
    
    cnt = 0;
    for i=1:size(sIds,1)-1
        for j=i+1:size(sIds,1)
            idx1 = find(strcmpi(sIds{i}, sampleIds));
            gt1 = reshape(genotypes(idx1,:),2,size(genotypes,2)/2)';

            idx2 = find(strcmpi(sIds{j}, sampleIds));
            gt2 = reshape(genotypes(idx2,:),2,size(genotypes,2)/2)';
   
            cnt = cnt + 1;
            coeffest{cnt,1} = sIds{i};
            coeffest{cnt,2} = sIds{j};
            coeffest{cnt,3} = num2str(calc_kinship_coef_distance(bps(indices), gt1(indices,:), gt2(indices,:), p, expectedDist));
        end;
    end;
end;

% between-cluster estimation
unique_cmbrship = unique(cmbrship);
for i=1:(length(unique_cmbrship)-1)
    c1 = unique_cmbrship(i);
    idxes1 = find(cmbrship == c1);
    sIds1 = sampleIds(idxes1);
    gts1 = genotypes(idxes1,:);
    
    trainSamplePairs = {};
    trainCoeffs = [];
    for j=i+1:length(unique_cmbrship)
        c2 = unique_cmbrship(j);
        idxes2 = find(cmbrship == c2);
        sIds2 = sampleIds(idxes2);
        gts2 = genotypes(idxes2,:);
        
        for k=1:size(samplePairsKnown,1)
            if isempty(find(strcmpi(samplePairsKnown{k,1},sIds1)))
                continue; 
            end;
        
            if isempty(find(strcmpi(samplePairsKnown{k,2},sIds2)))
                continue;
            end;
        
            trainSamplePairs = [trainSamplePairs;samplePairsKnown(k,1:2)];
            trainCoeffs = [trainCoeffs;coeffKnown(k)];
        end;
                
        if isempty(trainSamplePairs)
            clustersWOknownrel(c1).samples = [clustersWOknownrel(c1).samples;sIds];
            clustersWOknownrel(c2).samples = [clustersWOknownrel(c2).samples;sIds];
            continue;
        end;
        
        indices = sort(unique([find(afClusters(c1).af > 0);find(afClusters(c2).af > 0)]), 'ascend');        
        expectedDist = calculateExpectedDistance(bps(indices), afClusters(c1).af(indices), afClusters(c2).af(indices));
    
        genotypeIndices = sort([2*indices-1 2*indices],'ascend');
        p = estimateParam(sampleIds, genotypes(:,genotypeIndices), bps(indices), trainSamplePairs, trainCoeffs, expectedDist);    
            
        for k1=1:size(sIds1,1)
            for k2=1:size(sIds2,1)
                idx1 = find(strcmpi(sIds1{k1}, sampleIds));
                gt1 = reshape(genotypes(idx1,:),2,size(genotypes,2)/2)';

                idx2 = find(strcmpi(sIds2{k2}, sampleIds));
                gt2 = reshape(genotypes(idx2,:),2,size(genotypes,2)/2)';

                cnt = cnt + 1;
                coeffest{cnt,1} = sIds1{k1};
                coeffest{cnt,2} = sIds2{k2};
                coeffest{cnt,3} = num2str(calc_kinship_coef_distance(bps(indices), gt1(indices,:), gt2(indices,:), p, expectedDist));
            end;
        end;
    end;    
end;    
