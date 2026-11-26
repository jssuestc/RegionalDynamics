function Step03_SpatialCoupling(cfg)
% Step 03: single-participant spatial coupling and bilateral cortical spin.
P=cfg.Paths; [RAD,iRTD,ROI_IDs]=load_demo_metrics(cfg);
NSpin=5000; Seed=1;
assert(numel(unique(RAD))>1 && numel(unique(iRTD))>1,'Spatial correlation undefined for a constant map.');
if ~isfolder(P.Coupling), mkdir(P.Coupling); end

[SpinIndices,spinROI,Hemisphere,ROIName,Centroids]=BN210_spin_indices( ...
    cfg.AnnotationDir,cfg.SphereDir,NSpin,Seed);
assert(isequal(spinROI(:),ROI_IDs),'Spin ROI IDs do not match metrics.');
SpinParameters=struct('NSpin',NSpin,'Seed',Seed,'SphereDir',cfg.SphereDir, ...
    'AnnotationDir',cfg.AnnotationDir, ...
    'RotatedMap','iRTD','Tail','two-sided absolute rho','Algorithm','Vasa one-to-one');
save(P.SpinFile,'SpinIndices','ROI_IDs','Hemisphere','ROIName','Centroids','SpinParameters');
writetable(table(ROI_IDs,ROIName,Hemisphere,Centroids(:,1),Centroids(:,2),Centroids(:,3), ...
    'VariableNames',{'ROI_ID','ROI_Name','Hemisphere','SphereX','SphereY','SphereZ'}), ...
    fullfile(P.Coupling,'BN210_sphere_centroids.csv'));
[CouplingRho,CouplingPspin,NullR]=spin_spearman(RAD,iRTD,SpinIndices);
[~,P_raw]=corr(RAD,iRTD,'Type','Spearman');
SubjectID=string(cfg.SubjectID);
NROI=numel(ROI_IDs); MeanRAD=mean(RAD); MeaniRTD=mean(iRTD);
T=table(SubjectID,NROI,MeanRAD,MeaniRTD,CouplingRho,P_raw,CouplingPspin,NSpin);
writetable(T,fullfile(P.Coupling,'SpatialCoupling_subject.csv'));
Rotation=(1:NSpin)'; writetable(table(Rotation,NullR(:),'VariableNames',{'Rotation','NullRho'}), ...
    fullfile(P.Coupling,'Spin_null_distribution.csv'));
save(fullfile(P.Coupling,'RAD_iRTD_SpatialCoupling_results.mat'), ...
    'SubjectID','ROI_IDs','RAD','iRTD','CouplingRho','P_raw','CouplingPspin','NullR', ...
    'MeanRAD','MeaniRTD','NSpin','SpinParameters','cfg');
f=figure('Color','w','Visible','off'); scatter(iRTD,RAD,25,'filled'); box off;
xlabel('iRTD'); ylabel('RAD'); title(sprintf('Spearman rho = %.3f, P_{spin} = %.4g',CouplingRho,CouplingPspin));
exportgraphics(f,fullfile(P.Coupling,'SpatialCoupling_scatter.png'),'Resolution',300); close(f);
f=figure('Color','w','Visible','off'); histogram(NullR,40); hold on;
xline(CouplingRho,'r-','LineWidth',1.8); box off; xlabel('Spearman rho under spatial rotations'); ylabel('Count');
title(sprintf('%d spins, two-sided P = %.4g',NSpin,CouplingPspin));
exportgraphics(f,fullfile(P.Coupling,'Spin_null_distribution.png'),'Resolution',300); close(f);

end
