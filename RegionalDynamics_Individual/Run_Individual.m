%% Single-participant regional dynamics demo
clear; clc;

%% Input paths (relative to this script, independent of the working directory)
cfg.CodeRoot = fileparts(mfilename('fullpath'));
cfg.DataRoot = cfg.CodeRoot;
cfg.SubjectID = 'sub01';
cfg.FMRIFile = fullfile(cfg.DataRoot,'individual_data',[cfg.SubjectID '_RS.nii']);
cfg.AtlasFile = fullfile(cfg.DataRoot,'atlas_mask','BN_Atlas_246_5mm.nii');
cfg.MaskFile = fullfile(cfg.DataRoot,'atlas_mask','GreyMatter_mask.nii');
cfg.NeurosynthFile = fullfile(cfg.DataRoot,'atlas_mask','NeuroSynth124_BN_210.mat');
cfg.SPMDir = ''; % SPM12 installation folder, or leave empty if already on path.
cfg.AnnotationDir = fullfile(cfg.CodeRoot,'atlas');
cfg.SphereDir = fullfile(cfg.CodeRoot,'atlas_mask','fsaverage','surf');

%% Output paths (edit folder names here)
cfg.OutputRoot = fullfile(cfg.DataRoot,'results',cfg.SubjectID);
P.FCD = fullfile(cfg.OutputRoot,'01_FCD');
P.FCDMat = fullfile(P.FCD,'FCD_windows.mat');
P.Metrics = fullfile(cfg.OutputRoot,'02_RAD_iRTD');
P.MetricsMat = fullfile(P.Metrics,'RAD_iRTD_results.mat');
P.Coupling = fullfile(cfg.OutputRoot,'03_SpatialCoupling');
P.SpinFile = fullfile(P.Coupling,'SpinPermutations_used.mat');
P.Balance = fullfile(cfg.OutputRoot,'04_Balance');
P.PLS = fullfile(cfg.OutputRoot,'05_Pls_Neurosynth');
cfg.Paths = P;

%% FCD settings (unchanged)
cfg.FCDParameters = struct('WindowLength',50,'WindowStep',1,'FCDThreshold',0.6, ...
    'MaskThreshold',0,'IncludeSelfConnections',true,'ChunkSize',256);

%% Setup and input check
addpath(cfg.CodeRoot,fullfile(cfg.CodeRoot,'functions'),'-begin');
% If the configured directory is absent, try an installed FreeSurfer template.
if ~all(isfile(fullfile(cfg.SphereDir,{'lh.sphere.reg','rh.sphere.reg'})))
    fsCandidates = {};
    if ~isempty(getenv('SUBJECTS_DIR'))
        fsCandidates{end+1} = fullfile(getenv('SUBJECTS_DIR'),'fsaverage','surf');
    end
    if ~isempty(getenv('FREESURFER_HOME'))
        fsCandidates{end+1} = fullfile(getenv('FREESURFER_HOME'),'subjects','fsaverage','surf');
    end
    for k = 1:numel(fsCandidates)
        if all(isfile(fullfile(fsCandidates{k},{'lh.sphere.reg','rh.sphere.reg'})))
            cfg.SphereDir = fsCandidates{k};
            break;
        end
    end
end
requiredFiles = {cfg.FMRIFile,cfg.AtlasFile,cfg.MaskFile,cfg.NeurosynthFile, ...
    fullfile(cfg.AnnotationDir,'lh.BN_Atlas.annot'),fullfile(cfg.AnnotationDir,'rh.BN_Atlas.annot'), ...
    fullfile(cfg.AnnotationDir,'lh.cortex.label'),fullfile(cfg.AnnotationDir,'rh.cortex.label'), ...
    fullfile(cfg.SphereDir,'lh.sphere.reg'),fullfile(cfg.SphereDir,'rh.sphere.reg')};
missingFiles = requiredFiles(~isfile(requiredFiles));
assert(isempty(missingFiles),'Missing required input file(s):\n%s',strjoin(missingFiles,newline));
check_spm(cfg);

%% Run
Step01_FCD(cfg);
Step02_RAD_iRTD(cfg);
Step03_SpatialCoupling(cfg);
Step04_Balance(cfg);
Step05_Pls_Neurosynth(cfg,cfg.NeurosynthFile);
fprintf('Results: %s\n',cfg.OutputRoot);
