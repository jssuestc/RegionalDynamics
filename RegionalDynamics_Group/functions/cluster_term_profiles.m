function C = cluster_term_profiles(W,TermNames,fixedK)
% W: 124 terms x N participants, original bootstrap Z weights.
% Row z-score across participants; Euclidean k-means with 100 restarts
% for every candidate K. Choose the best mean silhouette (squared Euclidean).
% Retain the winning candidate itself rather than randomly refitting it.
% Optional fixedK bypasses K selection; silhouette is descriptive.
if nargin < 3, fixedK = []; end
W = double(W); TermNames = string(TermNames(:));
assert(size(W,1)==124 && size(W,2)>1 && numel(TermNames)==124 && ...
    numel(unique(TermNames))==124 && all(isfinite(W(:))) && ...
    all(std(W,0,2)>0),'Expected 124 complete nonconstant term profiles.');
M = normalize(W,2);
C.ValidTermMask = true(124,1); C.OriginalWeights = W;
C.NormalizedWeights = M; C.TermNames = TermNames;
if isempty(fixedK)
    C.KList = (2:8)';
else
    validateattributes(fixedK,{'numeric'},{'scalar','integer','>=',2,'<=',8});
    C.KList = double(fixedK);
end
nCandidates = numel(C.KList);
C.Silhouette = nan(nCandidates,1); C.TotalWithinSS = nan(nCandidates,1);
C.Replicates = 100; C.MaxIter = 1000; C.Seed = 2026;
C.CandidateLabels = zeros(124,nCandidates); C.CandidateSilhouettes = nan(124,nCandidates);
centroids = cell(nCandidates,1); sumd = cell(nCandidates,1); distances = cell(nCandidates,1);
oldRng = rng; cleanup = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(C.Seed,'twister');
for j = 1:numel(C.KList)
    k = C.KList(j);
    [idx,centroids{j},sumd{j},distances{j}] = kmeans(M,k, ...
        'Distance','sqeuclidean','Start','plus','Replicates',C.Replicates, ...
        'MaxIter',C.MaxIter,'EmptyAction','singleton');
    sv = silhouette(M,idx,'sqeuclidean');
    assert(all(isfinite(sv)), 'Undefined silhouette for K=%d.',k);
    C.CandidateLabels(:,j) = idx; C.CandidateSilhouettes(:,j) = sv;
    C.Silhouette(j) = mean(sv); C.TotalWithinSS(j) = sum(sumd{j});

end
if isempty(fixedK), [~,j] = max(C.Silhouette); else, j = 1; end
C.K = C.KList(j);
C.ClusterID = C.CandidateLabels(:,j); C.Centroids = centroids{j};
C.SumDistance = sumd{j}; C.Distance = distances{j};
C.CenterTermIndex = nan(C.K,1); C.CenterTerm = strings(C.K,1);
C.RToCenter = nan(124,1); C.RToCentroid = C.RToCenter; C.PToCenter = C.RToCenter;
C.OriginalClusterMean = nan(C.K,size(W,2));
for k = 1:C.K
    members = find(C.ClusterID==k); profiles = M(members,:)';
    r = corr(profiles,C.Centroids(k,:)','Type','Pearson');
    C.RToCentroid(members) = r;
    assert(any(isfinite(r)),'Undefined centroid correlation for cluster %d.',k);
    v = find(isfinite(r)); [~,i] = max(r(v)); i = v(i);
    C.CenterTermIndex(k) = members(i); C.CenterTerm(k) = TermNames(members(i));
    [C.RToCenter(members),C.PToCenter(members)] = corr(profiles,profiles(:,i),'Type','Pearson');
    C.OriginalClusterMean(k,:) = mean(W(members,:),1);
end
end
