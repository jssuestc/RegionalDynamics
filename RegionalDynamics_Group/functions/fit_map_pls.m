function R = fit_map_pls(X_raw,Y_raw,nCompResult,nCompCurve,SpinIndices)
% Joint symptom-map PLSR: paired ROI bootstrap and spatial spin inference.
% Two tested components and a 20-component variance curve.
% Both variance and fitted X-score/Y-score correlation are tested per spin.
% SpinIndices maps target rows to source rows in BN cortical ID order.
nBoot = 5000; nPerm = 5000; baseSeed = 20260714;
check_spin_indices(SpinIndices,(1:210)',nPerm);
assert(size(X_raw,1)==210 && size(X_raw,2)==124 && size(Y_raw,1)==210, ...
    'Expected 210 cortical rows and 124 predictors.');
assert(all(isfinite(Y_raw(:))) && all(std(Y_raw,0,1)>0), ...
    'Response maps must be finite and nonconstant.');
X = zscore(double(X_raw),0,1); Y = zscore(double(Y_raw),0,1);
[nROI,nTerm] = size(X); nY = size(Y,2);
[XL,YL,XS,YS,BETA,PCTVAR,MSE,st] = plsregress(X,Y,nCompCurve);
W_original = st.W(:,1:nCompResult);
score_r = nan(nCompResult,1);
for c = 1:nCompResult
    % For single Y, orient against the unchanged response (as in Demo).
    % For multivariate Y, retain the original symptom-script convention.
    if nY==1
        ref = Y;
    else
        ref = YS(:,c);
    end
    r = corr(XS(:,c),ref);
    assert(isfinite(r),'Undefined observed PLS score correlation.');
    if r<0
        W_original(:,c) = -W_original(:,c);
        XL(:,c) = -XL(:,c); YL(:,c) = -YL(:,c);
        XS(:,c) = -XS(:,c); YS(:,c) = -YS(:,c);
    end
    % Report actual fitted-score correlation, not the orientation reference.
    % Flipping both scores cannot change their mutual correlation.
    score_r(c) = corr(XS(:,c),YS(:,c),'Type','Pearson');
    assert(isfinite(score_r(c)),'Undefined observed X-score/Y-score correlation.');
end
YL_original = YL(:,1:nCompResult);
assert(all(isfinite(W_original(:))) && all(std(W_original,0,1)>0), ...
    'Invalid observed PLS weights.');
variance_each = 100*PCTVAR(2,:);
variance_cumulative = cumsum(variance_each);
assert(all(isfinite(variance_each)),'Invalid observed explained variance.');

%% Paired bootstrap (no re-standardization within each resample)
rng(baseSeed+1,'twister');
W_bootstrap = nan(nTerm,nCompResult,nBoot);
YL_bootstrap = nan(nY,nCompResult,nBoot);
bootstrap_indices = zeros(nROI,nBoot,'uint16');
bootstrap_alignment_r = nan(nCompResult,nBoot);
bootstrap_valid = false(nBoot,1);
for b = 1:nBoot
    idx = randsample(nROI,nROI,true);
    bootstrap_indices(:,b) = uint16(idx);
    try
        [~,yl,~,~,~,~,~,sb] = plsregress(X(idx,:),Y(idx,:),nCompResult);
        wb = sb.W(:,1:nCompResult);
        rr = nan(nCompResult,1);
        for c = 1:nCompResult
            rr(c) = corr(W_original(:,c),wb(:,c));
            if rr(c)<0, wb(:,c)=-wb(:,c); yl(:,c)=-yl(:,c); end
        end
        bootstrap_alignment_r(:,b) = rr;
        if all(isfinite(wb(:))) && all(isfinite(yl(:))) && all(isfinite(rr))
            W_bootstrap(:,:,b) = wb;
            YL_bootstrap(:,:,b) = yl;
            bootstrap_valid(b) = true;
        end
    catch ME
        if b==1, warning('Bootstrap fit failed: %s',ME.message); end
    end

end
nBootValid = sum(bootstrap_valid);
assert(nBootValid>=2,'Fewer than two valid bootstrap fits.');
if nBootValid<nBoot, warning('Valid bootstrap fits: %d/%d.',nBootValid,nBoot); end
W_bootstrap_SE = std(W_bootstrap(:,:,bootstrap_valid),0,3);
YL_bootstrap_SE = std(YL_bootstrap(:,:,bootstrap_valid),0,3);
W_bootstrap_SE(W_bootstrap_SE==0) = NaN;
YL_bootstrap_SE(YL_bootstrap_SE==0) = NaN;
W_z = W_original./W_bootstrap_SE; YL_z = YL_original./YL_bootstrap_SE;
W_p = 2*normcdf(-abs(W_z)); YL_p = 2*normcdf(-abs(YL_z));
W_bootstrap_CI95 = prctile(W_bootstrap(:,:,bootstrap_valid),[2.5 97.5],3);
YL_bootstrap_CI95 = prctile(YL_bootstrap(:,:,bootstrap_valid),[2.5 97.5],3);

%% Spin test: rotate all Y columns jointly, keeping X fixed
permutation_indices = uint16(SpinIndices); % Legacy field name; spatial indices.
variance_null_each = nan(nPerm,nCompResult);
score_r_null = nan(nPerm,nCompResult);
for p = 1:nPerm
    idx = SpinIndices(:,p);
    % Apply the same spatial rotation to every symptom column and refit.
    [~,~,xs_spin,ys_spin,~,pv] = plsregress(X,Y(idx,:),nCompResult);
    variance_null_each(p,:) = 100*pv(2,:);
    for c = 1:nCompResult
        score_r_null(p,c) = corr(xs_spin(:,c),ys_spin(:,c),'Type','Pearson');
    end

end
assert(all(isfinite(variance_null_each(:))), ...
    'Invalid spin fit: all 5000 rotations must produce finite variance.');
variance_null_cumulative = cumsum(variance_null_each,2);
assert(all(isfinite(variance_null_cumulative(:))), 'Invalid cumulative spin variance.');
valid_n_incremental = repmat(nPerm,1,nCompResult);
valid_n_cumulative = valid_n_incremental;
p_incremental = (1+sum(variance_null_each>=variance_each(1:nCompResult),1))/(nPerm+1);
p_cumulative = (1+sum(variance_null_cumulative>=variance_cumulative(1:nCompResult),1))/(nPerm+1);

% Two-sided, plus-one spatial null tests of the fitted score correlations.
score_spin_valid = isfinite(score_r_null);
valid_n_score = sum(score_spin_valid,1)';
assert(all(valid_n_score==nPerm), ...
    'Invalid score correlation: all 5000 spins must be finite for each component.');
score_spin_p = nan(nCompResult,1);
for c = 1:nCompResult
    score_spin_p(c) = (1+sum(abs(score_r_null(:,c))>=abs(score_r(c)))) ...
        /(1+valid_n_score(c));
end

%% Save complete estimates and resampling results through the caller
R = struct;
R.X_raw = X_raw; R.Y_raw = Y_raw; R.X = X; R.Y = Y;
R.Xloading = XL; R.Yloading = YL; R.XScore = XS; R.YScore = YS;
R.BETA = BETA; R.PCTVAR = PCTVAR; R.MSE = MSE;
R.W_original = W_original; R.YL_original = YL_original;
R.score_r = score_r; % Legacy observed correlation field.
R.ScoreR = score_r;
R.ScoreSpinP = score_spin_p;
R.score_r_null = score_r_null;
R.score_spin_valid = score_spin_valid;
R.ValidScoreSpinN = valid_n_score;
R.ScoreSpinMethod = 'Two-sided absolute Pearson r; refit PLS for each joint Y spin; plus-one';
R.W_bootstrap = W_bootstrap; R.YL_bootstrap = YL_bootstrap;
R.bootstrap_indices = bootstrap_indices; R.bootstrap_alignment_r = bootstrap_alignment_r;
R.bootstrap_valid = bootstrap_valid; R.nBootValid = nBootValid;
R.W_bootstrap_SE = W_bootstrap_SE; R.YL_bootstrap_SE = YL_bootstrap_SE;
R.W_bootstrap_CI95 = W_bootstrap_CI95; R.YL_bootstrap_CI95 = YL_bootstrap_CI95;
R.W_z = W_z; R.YL_z = YL_z; R.W_p = W_p; R.YL_p = YL_p;
R.variance_each = variance_each; R.variance_cumulative = variance_cumulative;
R.variance_null_each = variance_null_each; R.variance_null_cumulative = variance_null_cumulative;
R.permutation_indices = permutation_indices;
R.p_incremental = p_incremental; R.p_cumulative = p_cumulative;
R.valid_n_incremental = valid_n_incremental; R.valid_n_cumulative = valid_n_cumulative;
R.nBoot = nBoot; R.nPerm = nPerm; R.baseSeed = baseSeed;
R.nCompResult = nCompResult; R.nCompCurve = nCompCurve;
R.PermutationMethod = 'BN210 bilateral spin; joint Y-row rotation; X fixed';
R.SpinIndices = SpinIndices; R.nSpin = nPerm;
R.p_spin_incremental = p_incremental; R.p_spin_cumulative = p_cumulative;
R.VarianceSpinP = p_incremental(:);
R.CumulativeVarianceSpinP = p_cumulative(:);
R.BootstrapMethod = 'Paired ROI-row resampling; not spatially corrected';
R.ROI_IDs = (1:210)';
end
