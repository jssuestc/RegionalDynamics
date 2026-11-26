function run_symptom_map_pls(inputFiles,analysisNames,neurosynthFile,outputDir,SpinIndices,outputDirs)
% Use existing symptom-associated maps directly; do not recompute correlations.
% Fit all three symptom dimensions jointly in each disorder/metric model.
[X,TermNames] = load_reference_maps(neurosynthFile,124);
assert(numel(inputFiles)==4 && numel(analysisNames)==4,'Expected four symptom-map analyses.');
% Validate all inputs before the expensive resampling loops.
Y = cell(4,1); names = cell(4,1);
for a = 1:4
    S = load(inputFiles{a},'Y_all','symNames');
    assert(isfield(S,'Y_all') && isfield(S,'symNames'), ...
        'Expected Y_all and symNames in %s.',inputFiles{a});
    assert(isnumeric(S.Y_all) && isreal(S.Y_all) && isequal(size(S.Y_all),[210,3]), ...
        'Expected Y_all of size 210 x 3 in %s.',inputFiles{a});
    Y{a} = double(S.Y_all); names{a} = string(S.symNames(:));
    assert(all(isfinite(Y{a}(:))) && all(std(Y{a},0,1)>0) && ...
        numel(names{a})==3 && all(~ismissing(names{a})) && ...
        all(strlength(strtrim(names{a}))>0) && numel(unique(names{a}))==3, ...
        'Invalid symptom maps or labels: %s',inputFiles{a});
end
allSummary = cell(4,1);
for a = 1:4

    R = fit_map_pls(X,Y{a},2,20,SpinIndices);
    R.TermNames = TermNames; R.ResponseNames = names{a};
    R.symNames = cellstr(names{a}); R.VarNames = cellstr(TermNames);
    R.Analysis = analysisNames{a}; R.SourceFile = inputFiles{a}; R.NeurosynthFile = neurosynthFile;
    out = outputDirs{a};
    write_pls_outputs(R,out);
    T = readtable(fullfile(out,'PLS_summary.csv'));
    allSummary{a} = addvars(T,repmat(string(analysisNames{a}),height(T),1), ...
        'Before',1,'NewVariableNames','Analysis');
end
writetable(vertcat(allSummary{:}),fullfile(outputDir,'Symptom_PLS_summary.csv'));
end
