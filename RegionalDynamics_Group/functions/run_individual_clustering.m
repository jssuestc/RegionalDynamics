function run_individual_clustering(plsFile,TermNames,ageFile,outputDir,fixedK,plotByCluster)
% Cluster existing N x 124 bootstrap Z weights. Neurosynth file_names must
% match the column order used to create the supplied PLS1w_Z_individual.
% Optional fixedK: [] selects K by silhouette; 2 fixes K at two.
% Optional plotByCluster: true saves the per-cluster age-trajectory panels.
if nargin < 5, fixedK = []; end
if nargin < 6, plotByCluster = true; end
S = load(plsFile);
assert(isfield(S,'PLS1w_Z_individual'),'Missing PLS1w_Z_individual: %s',plsFile);
W = double(S.PLS1w_Z_individual);
assert(size(W,2)==124,'Expected subjects x 124 terms.');
if isfield(S,'TermNames')
    assert(isequal(string(S.TermNames(:)),string(TermNames(:))),'Term order mismatch.');
end
Age = []; SubjectID = []; % No invented IDs or automatic subject deletion.
if ~isempty(ageFile) && ~isfile(ageFile)
    warning('Age file not found; age fits skipped: %s',ageFile);
end
if ~isempty(ageFile) && isfile(ageFile)
    A = load(ageFile);
    if isfield(A,'Age'), Age = double(A.Age(:));
    elseif isfield(A,'AgeN'), Age = double(A.AgeN(:));
    else, error('Age file must contain Age or AgeN.'); end
    assert(numel(Age)==size(W,1) && all(isfinite(Age)) && numel(unique(Age))>=3, ...
        'Age must be aligned to PLS rows. Do not automatically delete row 231.');
end
C = cluster_term_profiles(W',TermNames,fixedK);
age_fit = [];
if ~isempty(Age), age_fit = linspace(min(Age),max(Age),200)'; end
% Reproduce the source figure's quadratic fits to normalized term profiles.
term_fit_coeff=nan(numel(C.TermNames),3); term_fit_y=nan(numel(C.TermNames),numel(age_fit));
centroid_fit_coeff=nan(C.K,3); centroid_fit_y=nan(C.K,numel(age_fit));
original_cluster_fit_y=centroid_fit_y;
if ~isempty(Age)
for t=1:numel(C.TermNames)
    term_fit_coeff(t,:)=polyfit(Age,C.NormalizedWeights(t,:)',2);
    term_fit_y(t,:)=polyval(term_fit_coeff(t,:),age_fit)';
end
for k=1:C.K
    centroid_fit_coeff(k,:)=polyfit(Age,C.Centroids(k,:)',2);
    centroid_fit_y(k,:)=polyval(centroid_fit_coeff(k,:),age_fit)';
    original_cluster_fit_y(k,:)=polyval(polyfit(Age,C.OriginalClusterMean(k,:)',2),age_fit)';
end
end
if ~isfolder(outputDir), mkdir(outputDir); end
Term=C.TermNames; Cluster=C.ClusterID; CenterTerm=C.CenterTerm(Cluster);
MeanOriginalZ=mean(C.OriginalWeights,2); Weight=C.RToCenter; RToCentroid=C.RToCentroid;
writetable(table(Term,Cluster,CenterTerm,Weight,RToCentroid,MeanOriginalZ),fullfile(outputDir,'Cluster_members.csv'));
writetable(table(C.KList,C.Silhouette,C.TotalWithinSS,'VariableNames',{'K','MeanSilhouette','WithinClusterSS'}),fullfile(outputDir,'Silhouette.csv'));
% Retain the main field names used by the historical Python plotting script.
age=Age; bestK=C.K; normalized_term_data=C.NormalizedWeights; original_term_data=C.OriginalWeights;
term_names=cellstr(Term); cluster_labels=C.ClusterID; cluster_centroids=C.Centroids;
cluster_centroid_subject_values=C.Centroids'; silhouette_k=C.KList; silhouette_score=C.Silhouette;
cluster_member_indices=cell(bestK,1); cluster_member_names=cell(bestK,1); cluster_member_count=zeros(bestK,1);
cluster_centroid_profiles=cell(bestK,1); cluster_original_term_zweight=cell(bestK,1); cluster_normalized_term_zweight=cell(bestK,1);
cluster_mean_original_term_zweight=cell(bestK,1); cluster_mean_normalized_term_zweight=cell(bestK,1);
for k=1:bestK
    j=find(cluster_labels==k); cluster_member_indices{k}=j; cluster_member_names{k}=term_names(j); cluster_member_count(k)=numel(j);
    cluster_centroid_profiles{k}=cluster_centroids(k,:); cluster_original_term_zweight{k}=original_term_data(j,:);
    cluster_normalized_term_zweight{k}=normalized_term_data(j,:);
    cluster_mean_original_term_zweight{k}=mean(original_term_data(j,:),2);
    cluster_mean_normalized_term_zweight{k}=mean(normalized_term_data(j,:),2);
end
cluster_center_term_indices=C.CenterTermIndex; cluster_center_term_names=cellstr(C.CenterTerm);
term_centroid_similarity_r=C.RToCentroid; term_center_similarity_r=C.RToCenter; term_center_similarity_p=C.PToCenter;
center_term_name_by_term=cellstr(CenterTerm); validTermMask=C.ValidTermMask; removedSubjectIndices=[]; polynomialDegree=2;
sumd=C.SumDistance; D=C.Distance;
term_fit_r2=nan(numel(Term),1); centroid_fit_r2=nan(bestK,1);
if ~isempty(Age)
for j=1:numel(Term), term_fit_r2(j)=fit_r2(Age,normalized_term_data(j,:)',term_fit_coeff(j,:)); end
for k=1:bestK, centroid_fit_r2(k)=fit_r2(Age,cluster_centroids(k,:)',centroid_fit_coeff(k,:)); end
end
ClusterID=(1:bestK)'; NTerms=cluster_member_count;
RepresentativeTerm=C.CenterTerm; CentroidFitR2=centroid_fit_r2;
writetable(table(ClusterID,NTerms,RepresentativeTerm,CentroidFitR2),fullfile(outputDir,'Cluster_summary.csv'));
if ~isempty(Age)
    T=table(age_fit,'VariableNames',{'Age'});
    for k=1:bestK
        T.(sprintf('Cluster_%d',k))=centroid_fit_y(k,:)';
    end
    writetable(T,fullfile(outputDir,'Cluster_age_trajectories.csv'));
end
save(fullfile(outputDir,'Cluster_Data_ForPython.mat'),'age','age_fit','normalized_term_data','original_term_data', ...
    'term_names','bestK','silhouette_k','silhouette_score','cluster_labels','cluster_centroids','cluster_centroid_subject_values', ...
    'cluster_member_indices','cluster_member_names','cluster_member_count','cluster_centroid_profiles','cluster_original_term_zweight', ...
    'cluster_normalized_term_zweight','cluster_mean_original_term_zweight','cluster_mean_normalized_term_zweight', ...
    'cluster_center_term_indices','cluster_center_term_names','term_centroid_similarity_r','term_center_similarity_r', ...
    'term_center_similarity_p','center_term_name_by_term','term_fit_coeff','term_fit_y','term_fit_r2','centroid_fit_coeff', ...
    'centroid_fit_y','centroid_fit_r2','validTermMask','removedSubjectIndices','polynomialDegree','sumd','D','original_cluster_fit_y','-v7');
save(fullfile(outputDir,'Cluster_results.mat'),'C','SubjectID','Age','age_fit','term_fit_coeff','term_fit_y','centroid_fit_coeff','centroid_fit_y','original_cluster_fit_y');
% Basic QC and representative-term word clouds.
f = figure('Color','w','Visible','off'); plot(C.KList,C.Silhouette,'o-','LineWidth',1.5);
xlabel('K'); ylabel('Mean silhouette'); xticks(C.KList); box off;
if ~isempty(fixedK), xlim([C.K-0.5,C.K+0.5]); title(sprintf('Fixed K = %d (descriptive silhouette)',C.K)); end
exportgraphics(f,fullfile(outputDir,'Silhouette.png'),'Resolution',300); close(f);
for k = 1:C.K
    j = find(C.ClusterID==k); words = C.TermNames(j); weights = C.RToCenter(j);
    plot_weights = max(weights,0); plot_weights(~isfinite(plot_weights)) = 0;
    f = figure('Color','w','Visible','off');
    wordcloud(words,plot_weights); title(sprintf('Cluster %d: %s',k,C.CenterTerm(k)));
    exportgraphics(f,fullfile(outputDir,sprintf('Cluster_%d_wordcloud.png',k)),'Resolution',300); close(f);
    writetable(table(words,weights,plot_weights),fullfile(outputDir,sprintf('Cluster_%d_wordcloud.csv',k)));
end
if ~isempty(Age)
    plot_cluster_age_trajectories(age_fit,term_fit_y,centroid_fit_y, ...
        cluster_labels,cluster_member_count,outputDir,plotByCluster);
end

end

function r2=fit_r2(age,y,b)
sst=sum((y-mean(y)).^2); if sst>0, r2=1-sum((y-polyval(b,age)).^2)/sst; else, r2=NaN; end
end

function plot_cluster_age_trajectories(age_fit,term_fit_y,centroid_fit_y,labels,counts,outputDir,plotByCluster)
% Quadratic fits to within-term standardized bootstrap Z profiles.
% The centroid curve is a fit to the mean normalized profile of its members.
% Age is used only for visualization, not as a clustering input.
K=numel(counts); colors=lines(K);
labelsText=compose('Cluster %d (%d terms)',(1:K)',counts(:));
% Combined centroid curves, with counts in the legend.
f=figure('Color','w','Visible','off','Position',[100 100 850 560]);
ax=axes(f); hold(ax,'on');
yline(ax,0,'--','Color',[0.72 0.72 0.72],'HandleVisibility','off');
h=gobjects(K,1);
for k=1:K
    h(k)=plot(ax,age_fit,centroid_fit_y(k,:)','Color',colors(k,:),'LineWidth',2.8);
end
xlabel(ax,'Age (years)'); ylabel(ax,'Normalized PLS1 term weight');
xlim(ax,[min(age_fit),max(age_fit)]);
set(ax,'Box','off','FontName','Arial','FontSize',13,'LineWidth',1,'TickDir','out');
legend(ax,h,cellstr(labelsText),'Location','best','Box','off');
exportgraphics(f,fullfile(outputDir,'Cluster_age_trajectories.png'),'Resolution',600); close(f);
if ~plotByCluster, return; end
% Per-cluster panels: individual term fits in light colors and centroid in bold.
% Use a common y-axis range across panels within this metric.
nCols=min(K,2); nRows=ceil(K/nCols);
f=figure('Color','w','Visible','off','Position',[100 100 610*nCols 390*nRows]);
t=tiledlayout(f,nRows,nCols,'TileSpacing','compact','Padding','compact');
allY=[0;term_fit_y(:);centroid_fit_y(:)];
ymin=min(allY); ymax=max(allY); pad=0.08*max(ymax-ymin,0.01);
for k=1:K
    ax=nexttile(t); hold(ax,'on');
    lightColor=0.82*[1 1 1]+0.18*colors(k,:);
    members=find(labels==k);
    for j=members(:)'
        plot(ax,age_fit,term_fit_y(j,:)','Color',lightColor,'LineWidth',0.8,'HandleVisibility','off');
    end
    yline(ax,0,'--','Color',[0.72 0.72 0.72],'HandleVisibility','off');
    plot(ax,age_fit,centroid_fit_y(k,:)','Color',colors(k,:),'LineWidth',3);
    title(ax,labelsText(k),'Color',colors(k,:),'FontWeight','bold');
    xlabel(ax,'Age (years)'); ylabel(ax,'Normalized PLS1 term weight');
    xlim(ax,[min(age_fit),max(age_fit)]); ylim(ax,[ymin-pad,ymax+pad]);
    set(ax,'Box','off','FontName','Arial','FontSize',12,'LineWidth',1,'TickDir','out');
end
exportgraphics(f,fullfile(outputDir,'Cluster_age_trajectories_by_cluster.png'),'Resolution',600); close(f);
end
