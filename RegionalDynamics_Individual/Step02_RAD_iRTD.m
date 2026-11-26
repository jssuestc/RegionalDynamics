function Step02_RAD_iRTD(cfg)
% Step 02: saved raw FCD -> RAD and iRTD across BN cortical IDs 1:210.
P=cfg.Paths;
S=load(P.FCDMat);
VA=spm_vol(cfg.AtlasFile); assert(numel(VA)==1,'Atlas must be 3D.');
ref=struct('dim',S.SpatialSize,'mat',S.Affine); check_grid(ref,VA,'BN atlas');
atlas=spm_read_vols(VA);
assert(all(isfinite(atlas(:)) & atlas(:)>=0 & atlas(:)==fix(atlas(:))), 'Atlas labels must be finite nonnegative integers.');
ROI_IDs=(1:210)'; labels=unique(atlas(:));
assert(all(ismember(ROI_IDs,labels)) && all(ismember(labels,0:246)), ...
    'Expected BN cortical IDs 1:210 and optional subcortical IDs 211:246.');
maskROI=atlas(S.MaskLinearIndices); target=ismember(maskROI,ROI_IDs);
TargetVoxelLinearIndices=S.MaskLinearIndices(target); TargetVoxelROIIDs=maskROI(target);
assert(any(target),'No cortical voxels overlap the FCD mask.');
% Use the original integer FCD counts directly, without any scaling.
FCDWindowPatterns=S.RawFCD(target,:);
M=calculate_RAD_iRTD_metrics(FCDWindowPatterns,TargetVoxelROIIDs,ROI_IDs);
if ~isfolder(P.Metrics), mkdir(P.Metrics); end
R=M; R.ROI_IDs=ROI_IDs; R.SubjectID=cfg.SubjectID; R.Parameters=cfg;
R.FCDParameters=S.FCDParameters; R.WindowGlobalMeanFCD=S.WindowGlobalMeanFCD;
R.MetricDefinition='raw FCD; RAD excludes zero FCD; iRTD uses fixed mask support';
R.WindowStarts=S.WindowStarts; R.WindowEnds=S.WindowEnds;
R.TargetVoxelLinearIndices=TargetVoxelLinearIndices; R.TargetVoxelROIIDs=TargetVoxelROIIDs;
R.FCDWindowPatterns=FCDWindowPatterns;
R.MeanRAD=mean(M.RAD); R.MeaniRTD=mean(M.iRTD); R.RawFCDFile=P.FCDMat;
R.CompleteCorticalCoverage=all(isfinite(M.RAD) & isfinite(M.iRTD));
save(P.MetricsMat,'-struct','R','-v7.3');
T=table(ROI_IDs,M.ROI_VoxelCount,M.MeanFCD,M.RAD,M.iRTD,M.RAD_ValidWindows, ...
    M.iRTD_ValidWindows,M.iRTD_ValidPairs,M.RAD_MinWindowVoxels, ...
    'VariableNames',{'ROI_ID','ValidVoxelCount','MeanFCD','RAD','iRTD','RAD_ValidWindows', ...
    'iRTD_ValidWindows','iRTD_ValidPairs','RAD_MinWindowVoxels'});
writetable(T,fullfile(P.Metrics,'RAD_iRTD_ROI.csv'));
SubjectID=string(cfg.SubjectID);
MeanRAD=R.MeanRAD; MeaniRTD=R.MeaniRTD; NValidROI=sum(isfinite(M.RAD)&isfinite(M.iRTD));
writetable(table(SubjectID,MeanRAD,MeaniRTD,NValidROI),fullfile(P.Metrics,'RAD_iRTD_summary.csv'));
write_atlas_map(VA,atlas,ROI_IDs,M.RAD,fullfile(P.Metrics,'RAD_map.nii'),'Regional aggregate dynamics');
write_atlas_map(VA,atlas,ROI_IDs,M.iRTD,fullfile(P.Metrics,'iRTD_map.nii'),'Intraregional topographic dynamics');
write_atlas_map(VA,atlas,ROI_IDs,M.MeanFCD,fullfile(P.Metrics,'MeanFCD_map.nii'),'Mean regional FCD');

if NValidROI<210, warning('Incomplete cortical coverage. Inspect RAD_iRTD_ROI.csv; coupling and balance require all 210 ROIs.'); end
if any(M.iRTD_ValidWindows<numel(S.WindowStarts)), warning('Some iRTD patterns are undefined. See iRTD_ValidWindows in ROI CSV.'); end

end
