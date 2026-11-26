function write_pls_outputs(R,outputDir)
% Save full PLS results, term weights, response loadings and null distributions.
if ~isfolder(outputDir), mkdir(outputDir); end
save(fullfile(outputDir,'PLS_results.mat'),'-struct','R');
for c = 1:R.nCompResult
    T = table(R.TermNames,R.W_original(:,c),R.W_bootstrap_SE(:,c), ...
        R.W_bootstrap_CI95(:,c,1),R.W_bootstrap_CI95(:,c,2),R.W_z(:,c),R.W_p(:,c), ...
        'VariableNames',{'Term','Weight','BootstrapSE','CI95_Lower','CI95_Upper','BootstrapZ','NormalApproxP'});
    writetable(T,fullfile(outputDir,sprintf('PLS%d_NeurosynthWeights.csv',c)));
    T = table(R.ResponseNames,R.YL_original(:,c),R.YL_bootstrap_SE(:,c), ...
        R.YL_bootstrap_CI95(:,c,1),R.YL_bootstrap_CI95(:,c,2),R.YL_z(:,c),R.YL_p(:,c), ...
        'VariableNames',{'Response','Loading','BootstrapSE','CI95_Lower','CI95_Upper','BootstrapZ','NormalApproxP'});
    writetable(T,fullfile(outputDir,sprintf('PLS%d_ResponseLoadings.csv',c)));
end
Component = (1:R.nCompResult)';
T = table(Component,reshape(R.variance_each(Component),[],1),reshape(R.variance_cumulative(Component),[],1), ...
    R.VarianceSpinP(:),R.CumulativeVarianceSpinP(:),R.ScoreR(:),R.ScoreSpinP(:), ...
    repmat(R.nBootValid,R.nCompResult,1),R.valid_n_incremental(:), ...
    R.valid_n_cumulative(:),R.ValidScoreSpinN(:),R.p_incremental(:),R.p_cumulative(:), ...
    'VariableNames',{'Component','VarianceExplained_pct','CumulativeVariance_pct', ...
    'VarianceSpinP','CumulativeVarianceSpinP','ScoreR','ScoreSpinP', ...
    'ValidBootstrapN','ValidSpinN','ValidCumulativeSpinN','ValidScoreSpinN', ...
    'SpinP','CumulativeSpinP'});
% SpinP and CumulativeSpinP are legacy aliases of the variance tests.
writetable(T,fullfile(outputDir,'PLS_summary.csv'));
T = table((1:R.nCompCurve)',R.variance_each(:),R.variance_cumulative(:), ...
    'VariableNames',{'Component','VarianceExplained_pct','CumulativeVariance_pct'});
writetable(T,fullfile(outputDir,'Variance_curve.csv'));
T = table((1:R.nPerm)','VariableNames',{'SpinID'});
for c = 1:R.nCompResult
    T.(sprintf('PLS%d_IncrementalVariance',c)) = R.variance_null_each(:,c);
    T.(sprintf('PLS%d_CumulativeVariance',c)) = R.variance_null_cumulative(:,c);
end
writetable(T,fullfile(outputDir,'Spin_null_variance.csv'));
T = table((1:R.nPerm)','VariableNames',{'SpinID'});
for c = 1:R.nCompResult
    T.(sprintf('PLS%d_ScoreR',c)) = R.score_r_null(:,c);
end
writetable(T,fullfile(outputDir,'Spin_null_scores.csv'));
end
