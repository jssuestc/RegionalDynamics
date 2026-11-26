function [counts,nValid] = compute_windowed_FCD(X,threshold,chunkSize)
% X: Nmask x window length. All masked voxels are both seeds and targets.
% Center and unit-normalize each time series; dot product = Pearson r.
% Zero-variance time series have undefined r; their connections count as 0.
% Use double precision and blocks to avoid an Nmask x Nmask allocation.
X=double(X);
nonconstant=any(X~=X(:,1),2);
X=X-mean(X,2);
norms=sqrt(sum(X.^2,2));
valid=nonconstant & isfinite(norms) & norms>0;
Z=zeros(size(X));
Z(valid,:)=X(valid,:)./norms(valid);
n=size(X,1); nValid=sum(valid);
counts=zeros(n,1,'uint32');
for first=1:chunkSize:n
    last=min(n,first+chunkSize-1); rows=first:last;
    C=Z(rows,:)*Z';
    d=sub2ind(size(C),(1:numel(rows))',rows(:));
    % Set diagonal explicitly; legacy includes self for nonconstant voxels.
    C(d)=double(valid(rows));
    counts(rows)=uint32(sum(C>threshold,2));
end
end
