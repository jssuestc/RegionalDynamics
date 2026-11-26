function M=calculate_RAD_iRTD_metrics(Patterns,TargetROIIDs,ROI_IDs)
% Patterns: cortical mask voxels x windows. Fixed voxel support for iRTD.
% RAD: sample SD of window-wise regional means; exclude zero FCD values, as in tvFCD.m.
% iRTD: 1 - mean Pearson correlation over all distinct valid window pairs.
nR=numel(ROI_IDs); nW=size(Patterns,2);
M.RAD=nan(nR,1); M.iRTD=nan(nR,1); M.MeanFCD=nan(nR,1);
M.ROI_VoxelCount=zeros(nR,1); M.WindowMeanFCD=nan(nR,nW);
M.RAD_ValidWindows=zeros(nR,1); M.iRTD_ValidWindows=zeros(nR,1);
M.iRTD_ValidPairs=zeros(nR,1); M.RAD_MinWindowVoxels=zeros(nR,1);
for r=1:nR
    Q=double(Patterns(TargetROIIDs==ROI_IDs(r),:)); nV=size(Q,1);
    M.ROI_VoxelCount(r)=nV;
    if nV<3, continue; end
    % Regional means: exclude zero FCD values as in the original RAD code.
    n=sum(Q~=0,1); A=sum(Q,1)./n;
    M.RAD_MinWindowVoxels(r)=min(n); M.WindowMeanFCD(r,:)=A;
    M.RAD_ValidWindows(r)=sum(isfinite(A));
    % Require every window for RAD; never silently drop empty windows.
    if all(isfinite(A))
        M.MeanFCD(r)=mean(A); M.RAD(r)=std(A,0,2);
    end
    % iRTD: fixed voxel set, including zero FCD. Undefined spatial
    % correlations (constant patterns) are omitted and counted explicitly.
    Q=Q-mean(Q,1); d=sqrt(sum(Q.^2,1)); valid=isfinite(d) & d>0;
    k=sum(valid); M.iRTD_ValidWindows(r)=k; M.iRTD_ValidPairs(r)=k*(k-1)/2;
    if k<2, continue; end
    Z=Q(:,valid)./d(valid); R=Z'*Z;
    vals=R(triu(true(k),1)); vals=max(-1,min(1,vals));
    M.iRTD(r)=1-mean(vals);
end
end
