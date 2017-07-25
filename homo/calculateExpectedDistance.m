function out = calculateExpectedDistance(bps, maf, varargin)

out = [];

if nargin > 3   % between-population
    maf2 = varargin{1};
    out = expectedDistanceBP(bps, maf, maf2);
else     % within-population
    out = expectedDistanceWP(bps, maf);
end;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% between-population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = expectedDistanceBP(bps, maf1, maf2)

out = [];

delta_d = abs(diff([(bps(1) - (1e-16));bps]));

% for population1
af_A_global = [1;maf1];
af_A_global = 1-af_A_global;
p1_A = af_A_global(1:end-1);
p2_A = af_A_global(2:end);
len_p = length(p1_A);

g0_A = p1_A.^2+(repmat(1,len_p,1) - p1_A.^2).*p2_A.^2;
gd_A = 2*p1_A.*(repmat(1,len_p,1) - p1_A).*(repmat(1,len_p,1) - p2_A.^2)+(repmat(1,len_p,1) - p1_A).^2*2.*p2_A.*(repmat(1,len_p,1) - p2_A);
g2d_A = (repmat(1,len_p,1) - p1_A).^2.*(repmat(1,len_p,1) - p2_A).^2;

% for population2
af_B_global = [1;maf2];
af_B_global = 1-af_B_global;
p1_B = af_B_global(1:end-1);
p2_B = af_B_global(2:end);

g0_B = p1_B.^2 + (repmat(1,len_p,1)-p1_B.^2).*p2_B.^2;
gd_B = 2*p1_B.*(repmat(1,len_p,1) - p1_B).*(repmat(1,len_p,1) - p2_B.^2) + (repmat(1,len_p,1) - p1_B).^2*2.*p2_B.*(repmat(1,len_p,1) - p2_B);
g2d_B = (repmat(1,len_p,1) - p1_B).^2.*(repmat(1,len_p,1) - p2_B).^2;

denominator = (repmat(1,len_p,1) - p1_A.^2.*p1_B.^2 - p2_A.^2.*p2_B.^2 + p1_A.^2.*p1_B.^2.*p2_A.^2.*p2_B.^2);
t0 = (g0_A.*g0_B - p1_A.^2.*p1_B.^2 - p2_A.^2.*p2_B.^2 + p1_A.^2.*p1_B.^2.*p2_A.^2.*p2_B.^2 + gd_A.*gd_B + g2d_A.*g2d_B)./denominator;
td = (g0_A.*gd_B + gd_A.*g0_B + gd_A.*g2d_B + g2d_A.*gd_B)./denominator;
t2d = (g0_A.*g2d_B + g2d_A.*g0_B)./denominator;

out = td.*delta_d + t2d.*delta_d*2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for within-population
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = expectedDistanceWP(bps, maf)

out = [];

delta_d = abs(diff([(bps(1) - (1e-16));bps]));
af_global = [1;maf];
af_global = 1-af_global;
p1 = af_global(1:end-1);
p2 = af_global(2:end);
len_p = length(p1);

g0 = p2.^2;
gd = 2*p2.*(repmat(1,len_p,1) - p2);
g2d = (repmat(1,len_p,1)-p2).^2;

denominator = repmat(1,len_p,1) - (p2.^4) + (p1.^4).*(p2.^4);
t0 = (g0.*g0+gd.*gd+g2d.*g2d-p2.^4+p1.^4.*p2.^4)./denominator; 
td = (2*g0.*gd+2*gd.*g2d)./denominator;
t2d = (2*g0.*g2d)./denominator;

out = td.*delta_d+t2d.*delta_d*2;

length(find(isnan(out)))