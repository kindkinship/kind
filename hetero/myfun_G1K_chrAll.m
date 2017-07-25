function fit_fun_val = myfun_G1K_chrAll(param, sampleIds, genotypes, samplePairsKnown, coeffKnown, expectedDist, bps)

fit_fun_val = [];

p = param;
%p

kinship_coef = [];
for i=1:size(samplePairsKnown,1)
    sid1 = samplePairsKnown{i,1};
    idx1 = find(strcmpi(sid1, sampleIds));
    gt1 = reshape(genotypes(idx1,:),2,size(genotypes,2)/2)';
        
    sid2 = samplePairsKnown{i,2};
    idx2 = find(strcmpi(sid2, sampleIds));
    gt2 = reshape(genotypes(idx2,:),2,size(genotypes,2)/2)';
        
    kc = calc_kinship_coef_distance(bps, gt1, gt2, p, expectedDist);

    kinship_coef = [kinship_coef;kc];
end;

fit_fun_val = [fit_fun_val;sum((kinship_coef - coeffKnown).^2)];
%fit_fun_val = sum(fit_fun_val)