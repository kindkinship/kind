function coeffest = estimateKinshipCoefficients(sampleIds, genotypes, bps, expectedDist, p)

coeffest = {};

cnt = 0;
for i=1:size(sampleIds,1)-1
    for j=i+1:size(sampleIds,1)
        idx1 = find(strcmpi(sampleIds{i}, sampleIds));
        gt1 = reshape(genotypes(idx1,:),2,size(genotypes,2)/2)';

        idx2 = find(strcmpi(sampleIds{j}, sampleIds));
        gt2 = reshape(genotypes(idx2,:),2,size(genotypes,2)/2)';
   
        cnt = cnt + 1;
        coeffest{cnt,1} = sampleIds{i};
        coeffest{cnt,2} = sampleIds{j};
        coeffest{cnt,3} = num2str(calc_kinship_coef_distance(bps, gt1, gt2, p, expectedDist));
    end;
end;