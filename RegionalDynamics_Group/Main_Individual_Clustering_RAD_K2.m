%% RAD Neurosynth term clustering with K fixed at 2
clear; clc;

%% Input and output paths (relative to this script)
codeRoot = fileparts(mfilename('fullpath'));
dataRoot = codeRoot;
radWeightFile = fullfile(dataRoot,'Individual_pls_RAD','PLS1w_Z.mat');
neurosynthFile = fullfile(dataRoot,'NeursynthMaps','NeuroSynth124_BN_210.mat');
radAgeFile = radWeightFile;
outputDir = fullfile(dataRoot,'results','02_Individual_Clustering','RAD_K2');

%% Clustering plot settings
plotByCluster = true; % false: skip Cluster_age_trajectories_by_cluster.png

%% Run fixed-K clustering
addpath(fullfile(codeRoot,'functions'),'-begin');
clear run_individual_clustering cluster_term_profiles
assert(strcmpi(which('run_individual_clustering'),fullfile(codeRoot,'functions','run_individual_clustering.m')) && ...
    strcmpi(which('cluster_term_profiles'),fullfile(codeRoot,'functions','cluster_term_profiles.m')) && ...
    nargin('run_individual_clustering')==6 && nargin('cluster_term_profiles')==3, ...
    'Clustering function version/path mismatch. Use the functions folder supplied with these main scripts.');
[~,TermNames] = load_reference_maps(neurosynthFile,124);
run_individual_clustering(radWeightFile,TermNames,radAgeFile,outputDir,2,plotByCluster);
