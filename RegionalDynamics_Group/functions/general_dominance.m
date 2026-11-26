function [weights,fullR2,relative,conditional,subsetR2,subsetMasks] = general_dominance(X,y)
% General dominance: average incremental ordinary R^2 across subset sizes.
% The intercept-only model has R^2=0. Five predictors give 32 subsets.
X = zscore(double(X),0,1);
y = zscore(double(y(:)),0,1);
[n,p] = size(X);
assert(n==numel(y) && all(isfinite(X(:))) && all(isfinite(y)), ...
    'Complete, nonconstant maps are required.');
subsetR2 = zeros(2^p,1);
subsetMasks = false(2^p,p);
for s = 0:2^p-1
    keep = logical(bitget(repmat(uint32(s),1,p),1:p));
    subsetMasks(s+1,:) = keep;
    if any(keep)
        mdl = fitlm(X(:,keep),y);
        subsetR2(s+1) = mdl.Rsquared.Ordinary;
    end
end
conditional = nan(p,p);
for j = 1:p
    for k = 0:p-1
        s = find(~subsetMasks(:,j) & sum(subsetMasks,2)==k)-1;
        conditional(k+1,j) = mean(subsetR2(s+2^(j-1)+1)-subsetR2(s+1));
    end
end
weights = mean(conditional,1)';
fullR2 = subsetR2(end);
assert(abs(sum(weights)-fullR2)<1e-8,'Dominance weights must sum to full-model R^2.');
relative = nan(p,1);
if fullR2>0, relative = 100*weights/fullR2; end
end
