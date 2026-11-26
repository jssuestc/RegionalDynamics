function Step01_FCD(cfg)
% Step 01: real 4D fMRI -> raw voxel-wise global FCD for all windows.
P=cfg.Paths;
FCDParameters=cfg.FCDParameters;
V=spm_vol(cfg.FMRIFile); VM=spm_vol(cfg.MaskFile); VA=spm_vol(cfg.AtlasFile);
assert(numel(VM)==1 && numel(VA)==1,'Mask and atlas must be 3D.');
assert(numel(V)>=FCDParameters.WindowLength+FCDParameters.WindowStep,'At least two complete windows are required.');
check_grid(V(1),V,'fMRI volumes'); check_grid(V(1),VM,'Grey-matter mask');
check_grid(V(1),VA,'BN atlas');
voxelSize=sqrt(sum(V(1).mat(1:3,1:3).^2,1));
assert(all(abs(voxelSize-5)<0.05),'Expected 5-mm isotropic input; check the atlas/fMRI grid.');
mask=spm_read_vols(VM);
MaskLinearIndices=find(isfinite(mask) & mask>FCDParameters.MaskThreshold);
NMaskVoxels=numel(MaskLinearIndices); assert(NMaskVoxels>=2,'Mask has fewer than two voxels.');
SpatialSize=V(1).dim(1:3); Affine=V(1).mat; NVolumes=numel(V);
X=zeros(NMaskVoxels,NVolumes);
for t=1:NVolumes
    vol=spm_read_vols(V(t)); X(:,t)=vol(MaskLinearIndices);
end
assert(all(isfinite(X(:))),'Nonfinite fMRI values inside the mask. Correct the input before analysis.');
NConstantScanVoxels=sum(all(X==X(:,1),2));
if NConstantScanVoxels>0
    warning('%d masked voxels have constant full-scan signals; retained in mask, undefined FC counts as zero.',NConstantScanVoxels);
end
WindowStarts=1:FCDParameters.WindowStep:(NVolumes-FCDParameters.WindowLength+1);
WindowEnds=WindowStarts+FCDParameters.WindowLength-1; nW=numel(WindowStarts);
RawFCD=zeros(NMaskVoxels,nW,'uint32');
WindowGlobalMeanFCD=zeros(1,nW); WindowValidSignalVoxels=zeros(1,nW);

for w=1:nW
    [RawFCD(:,w),WindowValidSignalVoxels(w)]=compute_windowed_FCD( ...
        X(:,WindowStarts(w):WindowEnds(w)),FCDParameters.FCDThreshold,FCDParameters.ChunkSize);
    WindowGlobalMeanFCD(w)=mean(double(RawFCD(:,w)));

end
SubjectID=cfg.SubjectID;
SourceFiles=struct('FMRIFile',cfg.FMRIFile,'MaskFile',cfg.MaskFile,'AtlasFile',cfg.AtlasFile);
if ~isfolder(P.FCD), mkdir(P.FCD); end
save(P.FCDMat,'RawFCD','MaskLinearIndices','WindowGlobalMeanFCD','WindowValidSignalVoxels', ...
    'WindowStarts','WindowEnds','NMaskVoxels','NVolumes','NConstantScanVoxels','SpatialSize','Affine', ...
    'SubjectID','FCDParameters','SourceFiles','-v7.3');
Window=(1:nW)';
T=table(Window,WindowStarts(:),WindowEnds(:),WindowGlobalMeanFCD(:),WindowValidSignalVoxels(:), ...
    'VariableNames',{'Window','StartVolume','EndVolume','GlobalMeanRawFCD','ValidSignalVoxels'});
writetable(T,fullfile(P.FCD,'FCD_window_summary.csv'));

end
