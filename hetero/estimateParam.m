function p = estimateParam(sampleIds, genotypes, bps, samplePairsKnown, coeffKnown, expectedDist)

p = NaN;

% optimization
LB = 1e-1;
UB = 1;

% initial guess
Bound = repmat([LB;UB],1,1);
x0 = unifrnd(LB,UB,1,1);
    
% with user-specifed gradient of the objective function        
options = optimset('UseParallel','always','Algorithm','active-set','FunValCheck','on');
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)myfun_G1K_chrAll(x, sampleIds, genotypes, samplePairsKnown, coeffKnown, expectedDist, bps), x0, [], [], [], [], ...
    Bound(1,:), Bound(2,:), [], options);   
p = x;