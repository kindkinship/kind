function out = calc_kinship_coef_distance(bps, genotype1, genotype2, p, expectedDist)

out = [];

[C,ia,ic] = unique(bps);
multiple_minor_alleles_indices = find(abs(diff(bps)) == 0);

if ~isempty(multiple_minor_alleles_indices)    
    genotype1(multiple_minor_alleles_indices,1) = genotype1(multiple_minor_alleles_indices,1) + genotype1(multiple_minor_alleles_indices+1,1); 
    genotype1(multiple_minor_alleles_indices,2) = genotype1(multiple_minor_alleles_indices,2) + genotype1(multiple_minor_alleles_indices+1,2);
    genotype1 = genotype1(ia,:);
    
    genotype2(multiple_minor_alleles_indices,1) = genotype2(multiple_minor_alleles_indices,1) + genotype2(multiple_minor_alleles_indices+1,1); 
    genotype2(multiple_minor_alleles_indices,2) = genotype2(multiple_minor_alleles_indices,2) + genotype2(multiple_minor_alleles_indices+1,2);
    genotype2 = genotype2(ia,:);    
    
    bps = bps(ia);
end;

num_ma_s1 = sum(genotype1');
num_ma_s2 = sum(genotype2');

s1 = zeros(1,size(genotype1,1));
s2 = s1;

% use the same bp positions used for ds to dt
% when a bp has two minor alleles and the previous bp has only one minor
% allele, calculate the distance for one of the two minor alleles but not
% both
nonzero_indices_s1 = find(num_ma_s1 > 0);
nonzero_indices_s2 = find(num_ma_s2 > 0);
nonzero_indices = sort(unique([nonzero_indices_s1 nonzero_indices_s2]));
bps_s1s2 = bps(nonzero_indices);

num_ma_s1 = num_ma_s1(nonzero_indices);
num_ma_s2 = num_ma_s2(nonzero_indices);

% positions where the number of minor allele is at least 1
indices_s1 = find(num_ma_s1 > 0);
s1(indices_s1) = abs(diff([(bps_s1s2(1) - (1e-16));bps_s1s2(indices_s1)]))';

indices_s2 = find(num_ma_s2 > 0);
s2(indices_s2) = abs(diff([(bps_s1s2(1) - (1e-16));bps_s1s2(indices_s2)]))';

% for S1
indices_s1_num_ma1 = find(num_ma_s1 == 1);
indices_s1_num_ma2 = find(num_ma_s1 == 2);
indices_s1 = indices_s1_num_ma2 - 1;
if indices_s1(1) == 0
    indices_s1 = indices_s1(2:end);
end;
indices_s1 = indices_s1(find(ismember(indices_s1, indices_s1_num_ma1) == 0));
indices_s1 = indices_s1 + 1;
if num_ma_s1(1) == 2
    indices_s1 = [indices_s1_num_ma2(1) indices_s1];
end;

% for S2
indices_s2_num_ma1 = find(num_ma_s2 == 1);
indices_s2_num_ma2 = find(num_ma_s2 == 2);
indices_s2 = indices_s2_num_ma2 - 1;
if indices_s2(1) == 0
    indices_s2 = indices_s2(2:end);
end;
indices_s2 = indices_s2(find(ismember(indices_s2, indices_s2_num_ma1) == 0));
indices_s2 = indices_s2 + 1;
if num_ma_s2(1) == 2
    indices_s2 = [indices_s2_num_ma2(1) indices_s2];
end;

% calculate the similarity score between S1 and S2
if ~isempty(indices_s1)
    s1(indices_s1) = s1(indices_s1) + abs(diff([(bps_s1s2(1) - 1);bps_s1s2(indices_s1)]))';
end;

if ~isempty(indices_s2)
    s2(indices_s2) = s2(indices_s2) + abs(diff([(bps_s1s2(1) - 1);bps_s1s2(indices_s2)]))';
end;

nominator = power(sum(abs(s1-s2).^p),1/p);
denominator = power(sum(abs(expectedDist).^p),1/p);

out = (1-(nominator/denominator))/2;