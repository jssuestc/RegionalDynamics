function Step05_Pls_Neurosynth(cfg, neurosynthFile)
% Participant-specific PLS1, paired ROI bootstrap and joint spatial null tests.
[RAD,iRTD,ROI_IDs] = load_demo_metrics(cfg);
assert(all(std([RAD,iRTD],0,1) > 0), ...
    'RAD and iRTD must both be nonconstant across the 210 cortical ROIs.');

%% 2. Load the original 124 Neurosynth maps
N = load(neurosynthFile);
if isfield(N,'roi_labels')
    labels = double(N.roi_labels(:));
    assert(numel(labels)>=210 && isequal(labels(1:210),(1:210)'), ...
        'Neurosynth cortical labels must follow BN atlas order 1:210.');
end
assert(isfield(N,'roi_matrix') && isfield(N,'file_names'), ...
    'Neurosynth MAT must contain roi_matrix and file_names.');
assert(isnumeric(N.roi_matrix) && isreal(N.roi_matrix) && ...
    ismatrix(N.roi_matrix) && size(N.roi_matrix,1) == 124 && ...
    ismember(size(N.roi_matrix,2),[210,246]), ...
    'Expected roi_matrix of size 124 x 210 or 124 x 246.');
% Match the original symptom-map script: select cortical columns and
% transpose to ROI x term, then standardize each term across ROIs.
X_raw = double(N.roi_matrix(:,1:210)');
assert(all(isfinite(X_raw(:))) && all(std(X_raw,0,1)>0), ...
    'Neurosynth maps must be finite and nonconstant.');
TermNames = string(N.file_names);
TermNames = TermNames(:);
assert(numel(TermNames)==124 && ~any(ismissing(TermNames)) && ...
    all(strlength(strtrim(TermNames))>0), 'Expected 124 valid term names.');
X = zscore(X_raw,0,1);

%% 3. Fixed analysis settings and output directory
nBoot = 5000;
nPerm = 5000;
baseSeed = 20260714;  % Same base seed as the symptom-map PLSR script.
outputDir = cfg.Paths.PLS;
if ~isfolder(outputDir), mkdir(outputDir); end
% Reuse the Step03 rotation indices.
spinFile = cfg.Paths.SpinFile;
S = get_bn210_spins(spinFile,cfg.AnnotationDir,cfg.SphereDir);
SpinIndices = S.SpinIndices;
assert(isequal(double(S.ROI_IDs(:)),double(ROI_IDs(:))), 'Spin and metric ROI order mismatch.');
save(fullfile(outputDir,'SpinPermutations_used.mat'),'-struct','S');
SubjectID = string(cfg.SubjectID);
Metric = ["RAD";"iRTD"];
Y_all = [RAD,iRTD];
VarianceExplained_pct = nan(2,1);
PermutationP = nan(2,1);
ScoreR = nan(2,1);
ScoreSpinP = nan(2,1);
ValidScoreSpinN = zeros(2,1);
ValidBootstrapN = zeros(2,1);
ValidPermutationN = zeros(2,1);

%% 4. Fit two independent PLS1 models and save all results
for m = 1:2
    name = char(Metric(m));

    R = fit_pls1(X,Y_all(:,m),nBoot,nPerm,baseSeed,SpinIndices);
    R.SubjectID = SubjectID;
    R.Metric = Metric(m);
    R.ROI_IDs = ROI_IDs;
    R.TermNames = TermNames;
    R.VarNames = TermNames;  % Retain the legacy term-name convention.
    R.X_raw = X_raw;
    R.X = X;
    R.NeurosynthFile = neurosynthFile;
    R.MetricsFile = cfg.Paths.MetricsMat;
    % Row vectors are compatible with the existing subject-by-term
    % clustering inputs when results are collected across participants.
    R.PLS1w_individual = R.W_original';
    R.PLS1w_Z_individual = R.W_z';
    R.PLS1_variance = R.VarianceExplained_pct/100; % Fraction, not percent.
    R.PLS1_perm_p = R.PermutationP; % Legacy field; now spatial p.
    R.PLS1_spin_p = R.SpinP;
    R.PLS1_score_r = R.ScoreR;
    R.PLS1_score_spin_p = R.ScoreSpinP;
    save(fullfile(outputDir,[name,'_PLS1_results.mat']),'-struct','R');

    T = table(TermNames,R.W_original,R.W_bootstrap_SE, ...
        R.W_bootstrap_CI95(:,1),R.W_bootstrap_CI95(:,2),R.W_z,R.W_p, ...
        abs(R.W_z)>=1.96,abs(R.W_z)>=2.58, ...
        'VariableNames',{'Term','Weight','BootstrapSE','BootstrapCI95_Lower', ...
        'BootstrapCI95_Upper','BootstrapZ','BootstrapNormalApproxP', ...
        'AbsZ_GE_1_96','AbsZ_GE_2_58'});
    writetable(T,fullfile(outputDir,[name,'_PLS1_NeurosynthWeights.csv']));
    Tnull = table((1:nPerm)',R.variance_null_each,R.score_r_null, ...
        'VariableNames',{'SpinID','VarianceExplained_pct','ScoreR'});
    writetable(Tnull,fullfile(outputDir,[name,'_PLS1_Spin_null.csv']));
    VarianceExplained_pct(m) = R.VarianceExplained_pct;
    PermutationP(m) = R.PermutationP;
    ScoreR(m) = R.ScoreR;
    ScoreSpinP(m) = R.ScoreSpinP;
    ValidScoreSpinN(m) = R.ValidScoreSpinN;
    ValidBootstrapN(m) = R.ValidBootstrapN;
    ValidPermutationN(m) = R.ValidPermutationN;

end
T = table(repmat(SubjectID,2,1),Metric,VarianceExplained_pct,PermutationP, ...
    ScoreR,ScoreSpinP,ValidBootstrapN,ValidPermutationN,ValidScoreSpinN, ...
    PermutationP,'VariableNames', ...
    {'SubjectID','Metric','VarianceExplained_pct','VarianceSpinP', ...
    'ScoreR','ScoreSpinP','ValidBootstrapN','ValidSpinN','ValidScoreSpinN','SpinP'});
% SpinP is retained as a legacy alias of VarianceSpinP.
writetable(T,fullfile(outputDir,'PLS1_summary.csv'));

end

function R = fit_pls1(X,Y_raw,nBoot,nPerm,baseSeed,SpinIndices)
% Single-response PLS1; use the same standardization, paired bootstrap and
% spatial-spin framework for RAD and iRTD.
nROI = size(X,1);
nTerm = size(X,2);
Y = zscore(double(Y_raw(:)),0,1);
assert(all(isfinite(Y)), 'Invalid standardized response.');

%% Observed PLS1
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,1);
W_original = stats.W(:,1);
% Fix component orientation against the unchanged response. Flipping both
% X and Y scores together would not change their mutual correlation.
score_r = corr(XS(:,1),Y,'Type','Pearson');
assert(isfinite(score_r) && all(isfinite(W_original)) && std(W_original)>0, ...
    'Observed PLS1 is undefined or has a constant term-weight vector.');
if score_r < 0
    W_original = -W_original;
    XL = -XL; YL = -YL; XS = -XS; YS = -YS;
    score_r = -score_r;
end
% Use the actual fitted X and Y scores for the observed statistic.
% The original orientation against Y above is unchanged.
ScoreR = corr(XS(:,1),YS(:,1),'Type','Pearson');
assert(isfinite(ScoreR),'Undefined observed X-score/Y-score correlation.');
VarianceExplained_pct = 100*PCTVAR(2,1);
assert(isfinite(VarianceExplained_pct),'Observed explained variance is undefined.');

%% Bootstrap: resample paired X/Y cortical rows with replacement
% Do not re-standardize within each resample, matching the symptom script.
rng(baseSeed+1,'twister');
W_bootstrap = nan(nTerm,nBoot);
bootstrap_indices = zeros(nROI,nBoot,'uint16');
bootstrap_alignment_r = nan(nBoot,1);
bootstrap_valid = false(nBoot,1);
for b = 1:nBoot
    idx = randsample(nROI,nROI,true);
    bootstrap_indices(:,b) = uint16(idx);
    Xb = X(idx,:); Yb = Y(idx,:);
    if std(Yb)>0
        [~,~,~,~,~,~,~,st] = plsregress(Xb,Yb,1);
        wb = st.W(:,1);
        rr = corr(W_original,wb,'Type','Pearson');
        bootstrap_alignment_r(b) = rr;
        if all(isfinite(wb)) && isfinite(rr)
            if rr<0, wb=-wb; end
            W_bootstrap(:,b) = wb;
            bootstrap_valid(b) = true;
        end
    end

end
ValidBootstrapN = sum(bootstrap_valid);
assert(ValidBootstrapN>=2,'Fewer than two valid bootstrap solutions.');
if ValidBootstrapN<nBoot
    warning('Only %d/%d bootstrap solutions were valid.',ValidBootstrapN,nBoot);
end
W_bootstrap_SE = std(W_bootstrap(:,bootstrap_valid),0,2);
W_bootstrap_SE(W_bootstrap_SE==0) = NaN;
W_bootstrap_CI95 = prctile(W_bootstrap(:,bootstrap_valid),[2.5 97.5],2);
W_z = W_original./W_bootstrap_SE;
% Descriptive normal approximation to the bootstrap ratio, not a spatially
% corrected or multiplicity-adjusted term-level significance test.
W_p = 2*normcdf(-abs(W_z));
if any(~isfinite(W_z)), warning('Some term bootstrap Z values are undefined.'); end

%% Spin test: rotate Y on the cortical sphere while keeping X fixed
check_spin_indices(SpinIndices,(1:210)',nPerm);
permutation_indices = uint16(SpinIndices); % Legacy name; spatial indices.
variance_null_each = nan(nPerm,1);
score_r_null = nan(nPerm,1);
for p = 1:nPerm
    idx = SpinIndices(:,p);
    % Refit both scores from the rotated original Y (not a post-fit spin).
    [~,~,xs_spin,ys_spin,~,pv] = plsregress(X,Y(idx,:),1);
    variance_null_each(p) = 100*pv(2,1);
    score_r_null(p) = corr(xs_spin(:,1),ys_spin(:,1),'Type','Pearson');

end
assert(all(isfinite(variance_null_each)), ...
    'Invalid spin fit: all 5000 rotations must produce finite variance.');
permutation_valid = isfinite(variance_null_each);
ValidPermutationN = sum(permutation_valid);
PermutationP = (1+sum(variance_null_each(permutation_valid)>=VarianceExplained_pct)) ...
    /(1+ValidPermutationN);

% Two-sided, plus-one spatial null test for the fitted score correlation.
score_spin_valid = isfinite(score_r_null);
ValidScoreSpinN = sum(score_spin_valid);
assert(ValidScoreSpinN==nPerm, ...
    'Invalid score correlation: all 5000 spins must be finite.');
ScoreSpinP = (1+sum(abs(score_r_null)>=abs(ScoreR)))/(1+ValidScoreSpinN);

%% Collect observed estimates and complete resampling distributions
R = struct;
R.Y_raw = Y_raw(:); R.Y = Y;
R.Xloading = XL; R.Yloading = YL; R.XScore = XS; R.YScore = YS;
R.Beta = BETA; R.PCTVAR = PCTVAR; R.MSE = MSE;
R.ScoreR = ScoreR;
R.ScoreSpinP = ScoreSpinP;
R.score_r_null = score_r_null;
R.score_spin_valid = score_spin_valid;
R.ValidScoreSpinN = ValidScoreSpinN;
R.ScoreSpinMethod = 'Two-sided absolute Pearson r; refit PLS for each Y spin; plus-one';
R.W_original = W_original;
R.W_bootstrap = W_bootstrap;
R.bootstrap_indices = bootstrap_indices;
R.bootstrap_alignment_r = bootstrap_alignment_r;
R.bootstrap_valid = bootstrap_valid;
R.ValidBootstrapN = ValidBootstrapN;
R.W_bootstrap_SE = W_bootstrap_SE;
R.W_bootstrap_CI95 = W_bootstrap_CI95;
R.W_z = W_z; R.W_p = W_p;
R.VarianceExplained_pct = VarianceExplained_pct;
R.variance_null_each = variance_null_each;
R.permutation_indices = permutation_indices;
R.permutation_valid = permutation_valid;
R.ValidPermutationN = ValidPermutationN;
R.PermutationP = PermutationP;
R.VarianceSpinP = PermutationP;
R.nBoot = nBoot; R.nPerm = nPerm; R.baseSeed = baseSeed;
R.nComponents = 1;
R.PermutationMethod = 'BN210 bilateral spin of Y; X fixed';
R.SpinIndices = SpinIndices; R.nSpin = nPerm; R.SpinP = PermutationP;
R.BootstrapMethod = 'Paired ROI-row resampling; not spatially corrected';
R.SignAlignment = 'Observed X-score versus Y; bootstrap weight correlation';
end
