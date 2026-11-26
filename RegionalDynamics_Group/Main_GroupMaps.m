%% Group-map analyses
clear; clc;

%% Input paths (relative to this script, independent of the working directory)
codeRoot = fileparts(mfilename('fullpath'));
repoRoot = fileparts(codeRoot);
dataRoot = codeRoot;
radMeanFile = fullfile(dataRoot,'AgeBin_RAD','Age_00-90_Mean.mat');
irtMeanFile = fullfile(dataRoot,'AgeBin_iRTD','Age_00-90_Mean.mat');
neuromapsFile = fullfile(dataRoot,'Neuromaps','Neuromaps_BN_210.mat');
neurosynthFile = fullfile(dataRoot,'NeursynthMaps','NeuroSynth124_BN_210.mat');
radWeightFile = fullfile(dataRoot,'Individual_pls_RAD','PLS1w_Z.mat');
irtWeightFile = fullfile(dataRoot,'Individual_pls_iRTD','PLS1w_Z.mat');
% Each supplied PLS weight MAT also contains its aligned Age vector.
% Set either age file to '' to skip the corresponding age fits.
radAgeFile = radWeightFile;
irtAgeFile = irtWeightFile;
symptomFiles = {
    fullfile(dataRoot,'SymptomMap_RAD_MDD','Symptom_Map.mat');
    fullfile(dataRoot,'SymptomMap_iRTD_MDD','Symptom_Map.mat');
    fullfile(dataRoot,'SymptomMap_RAD_SZ','Symptom_Map.mat');
    fullfile(dataRoot,'SymptomMap_iRTD_SZ','Symptom_Map.mat')};
analysisNames = {'MDD_RAD';'MDD_iRTD';'SZ_RAD';'SZ_iRTD'};
atlasDir = fullfile(codeRoot,'atlas');
sphereDir = fullfile(repoRoot,'RegionalDynamics_Individual','atlas_mask','fsaverage','surf');

%% Output paths (edit folder names here)
outputRoot = fullfile(dataRoot,'results');
dominanceDir = fullfile(outputRoot,'01_Dominance');
clusteringDir = fullfile(outputRoot,'02_Individual_Clustering');
radClusterDir = fullfile(clusteringDir,'RAD');
irtClusterDir = fullfile(clusteringDir,'iRTD');
symptomDir = fullfile(outputRoot,'03_Symptom_PLS');
symptomOutputDirs = {
    fullfile(symptomDir,'MDD_RAD');
    fullfile(symptomDir,'MDD_iRTD');
    fullfile(symptomDir,'SZ_RAD');
    fullfile(symptomDir,'SZ_iRTD')};
spinFile = fullfile(symptomDir,'SpinPermutations_used.mat'); % Or point to a validated existing spin file.

%% Clustering plot settings
plotByCluster = true; % false: skip Cluster_age_trajectories_by_cluster.png

%% Setup (run this section before individual analysis sections)
addpath(fullfile(codeRoot,'functions'),'-begin');
clear run_individual_clustering cluster_term_profiles
assert(strcmpi(which('run_individual_clustering'),fullfile(codeRoot,'functions','run_individual_clustering.m')) && ...
    strcmpi(which('cluster_term_profiles'),fullfile(codeRoot,'functions','cluster_term_profiles.m')) && ...
    nargin('run_individual_clustering')==6 && nargin('cluster_term_profiles')==3, ...
    'Clustering function version/path mismatch. Use the functions folder supplied with these main scripts.');
requiredFiles = [{radMeanFile,irtMeanFile,neuromapsFile,neurosynthFile,radWeightFile,irtWeightFile},symptomFiles(:)'];
if ~isempty(radAgeFile), requiredFiles{end+1} = radAgeFile; end
if ~isempty(irtAgeFile), requiredFiles{end+1} = irtAgeFile; end
missingFiles = requiredFiles(~isfile(requiredFiles));
assert(isempty(missingFiles),'Missing required input file(s):\n%s',strjoin(missingFiles,newline));
% A valid existing spin cache does not require the external sphere files.
if ~isfile(spinFile)
    % If the configured directory is absent, try an installed FreeSurfer template.
    if ~all(isfile(fullfile(sphereDir,{'lh.sphere.reg','rh.sphere.reg'})))
        fsCandidates = {};
        if ~isempty(getenv('SUBJECTS_DIR'))
            fsCandidates{end+1} = fullfile(getenv('SUBJECTS_DIR'),'fsaverage','surf');
        end
        if ~isempty(getenv('FREESURFER_HOME'))
            fsCandidates{end+1} = fullfile(getenv('FREESURFER_HOME'),'subjects','fsaverage','surf');
        end
        for k = 1:numel(fsCandidates)
            if all(isfile(fullfile(fsCandidates{k},{'lh.sphere.reg','rh.sphere.reg'})))
                sphereDir = fsCandidates{k};
                break;
            end
        end
    end
    requiredFiles = {fullfile(atlasDir,'lh.BN_Atlas.annot'),fullfile(atlasDir,'rh.BN_Atlas.annot'), ...
        fullfile(atlasDir,'lh.cortex.label'),fullfile(atlasDir,'rh.cortex.label'), ...
        fullfile(sphereDir,'lh.sphere.reg'),fullfile(sphereDir,'rh.sphere.reg')};
    missingFiles = requiredFiles(~isfile(requiredFiles));
    assert(isempty(missingFiles),'Missing required input file(s):\n%s',strjoin(missingFiles,newline));
end

%% Dominance
run_group_dominance(radMeanFile,irtMeanFile,neuromapsFile,dominanceDir);

%% Individual clustering
[~,TermNames] = load_reference_maps(neurosynthFile,124);
A = load(radWeightFile,'PLS1w_Z_individual');
B = load(irtWeightFile,'PLS1w_Z_individual');
assert(isequal(size(A.PLS1w_Z_individual),size(B.PLS1w_Z_individual)) && ...
    ~isequal(A.PLS1w_Z_individual,B.PLS1w_Z_individual), 'Check RAD/iRTD weight inputs.');
run_individual_clustering(radWeightFile,TermNames,radAgeFile,radClusterDir,[],plotByCluster);
run_individual_clustering(irtWeightFile,TermNames,irtAgeFile,irtClusterDir,[],plotByCluster);

%% Symptom PLS
S = get_bn210_spins(spinFile,atlasDir,sphereDir);
if ~isfolder(symptomDir), mkdir(symptomDir); end
run_symptom_map_pls(symptomFiles,analysisNames,neurosynthFile,symptomDir,S.SpinIndices,symptomOutputDirs);
fprintf('Results: %s\n',outputRoot);
