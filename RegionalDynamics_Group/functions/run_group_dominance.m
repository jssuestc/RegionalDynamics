function run_group_dominance(radFile,irtFile,neuromapsFile,outputDir)
% General dominance of the supplied all-age mean RAD and iRTD maps.
% Each response is the pooled cortical map, not a mean of age-bin means.
[X,MapNames] = load_reference_maps(neuromapsFile,5);
expectedNames = ["T1T2ratio";"G1_fMRI";"G1_dMRI";"PET_AG";"Cortical_Thickness"];
assert(isequal(sort(MapNames),sort(expectedNames)), ...
    'Expected T1T2ratio, G1_fMRI, G1_dMRI, PET_AG and Cortical_Thickness. Use the corrected Neuromaps file.');
% Preserve the column order and names of the supplied Neuromaps matrix.
assert(~strcmpi(radFile,irtFile),'Specify separate RAD and iRTD input files.');
RAD = load_mean_map(radFile);
iRTD = load_mean_map(irtFile);
assert(~isequal(RAD,iRTD), ...
    'RAD and iRTD maps are identical. Check that the two pooled maps are not duplicate files.');
if ~isfolder(outputDir), mkdir(outputDir); end
metricNames = {'RAD','iRTD'};
sourceFiles = {radFile,irtFile};
D = {RAD,iRTD};
allTables = cell(2,1);
for m = 1:2
    metric = metricNames{m};
    Y = D{m};
    [weights,fullR2,relative,conditional,subsetR2,subsetMasks] = general_dominance(X,Y);
    T = table(repmat(string(metric),5,1),repmat("Age_00-90",5,1), ...
        MapNames,weights,100*weights,relative,repmat(fullR2,5,1), ...
        'VariableNames',{'Metric','Model','Map','Dominance_R2','ExplainedVariance_pct', ...
        'ShareOfExplainedVariance_pct','FullModel_R2'});
    allTables{m} = T;

    R = struct('Metric',metric,'Model','Age_00-90', ...
        'X_raw',X,'Y_raw',Y,'MapNames',MapNames,'weights',weights, ...
        'fullR2',fullR2,'relative',relative,'conditional',conditional, ...
        'subsetR2',subsetR2,'subsetMasks',subsetMasks,'ROI_IDs',(1:210)');
    R.sourceFile = sourceFiles{m}; R.neuromapsFile = neuromapsFile;
    R.MeanDefinition = 'Supplied all-age mean map (Age_00-90_Mean.mat); no age-bin averaging is performed.';
    prefix = [metric,'_Age_00-90_Dominance'];
    save(fullfile(outputDir,[prefix,'.mat']),'-struct','R');
    writetable(T,fullfile(outputDir,[prefix,'.csv']));
    f = figure('Color','w','Visible','off');
    bar(100*weights);
    xticks(1:5); xticklabels(MapNames); xtickangle(30);
    ylabel('Explained spatial variance (%)');
    title(sprintf('%s: all-age mean (R^2 = %.3f)',metric,fullR2)); box off;
    exportgraphics(f,fullfile(outputDir,[prefix,'.png']),'Resolution',300);
    close(f);
end
% Overwrite the previous summary with only the two all-age models.
writetable(vertcat(allTables{:}),fullfile(outputDir,'Dominance_all_models.csv'));
end
