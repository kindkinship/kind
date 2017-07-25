function kinship_coef = calculate_kinshp_coef_chrAll_WP_v2(D, p, expected_dist, af)

kinship_coef = [];
for j=1:size(D.pair_indices,1)
    myidx = D.pair_indices(j,1);
    pidx = D.pair_indices(j,2);
    
    snp_positions = [];    
    genotype1 = [];
    genotype2 = [];
    
    for k = 1:size(D.chr,2)
        if k == 1
            snp_positions = D.chr(k).snp_positions;
        else
            snp_positions = [snp_positions;D.chr(k).snp_positions-D.chr(k).snp_positions(1)+1+snp_positions(end)];
        end;
                    
        genotype1 = [genotype1;D.chr(k).D_genotype(myidx).genotype];
        genotype2 = [genotype2;D.chr(k).D_genotype(pidx).genotype];
    end;    
    
    kinship_coef = [kinship_coef;calc_kinship_coef_distance_v2(snp_positions, genotype1, genotype2, p, expected_dist, af)];
end;