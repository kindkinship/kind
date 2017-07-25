% calculateExpectedDistanceCluster.m

maf = [];


maf = [];
if strcmp(prefix_data,'HapMap')
    maf = calculate_MAF_cluster_HapMap(prefix_data, mbr_ids, num_chromosomes, target_pops);
    chromosomes_of_interest = [1:num_chromosomes];
elseif strcmp(prefix_data,'G1K')   
    maf = calculate_MAF_cluster_G1K(prefix_data, mbr_ids, chromosomes_of_interest, target_pops, blk_size);
end;

sps = [];
for k = 1:numel(chromosomes_of_interest)
    chrno = chromosomes_of_interest(k);
    sp_fname = [prefix_data '_chr' num2str(chrno) '_' target_pops '_snp_positions.mat'];
    load(sp_fname);

    if k == 1
        sps = snp_positions;
    else
        sps = [sps;snp_positions-snp_positions(1)+1+sps(end)];
    end;
end;

indices = find(maf > 0);
expected_dist(c).dist = calculateExpectedDistance(sps, maf, indices);    
save(expected_dist_fname, 'expected_dist');