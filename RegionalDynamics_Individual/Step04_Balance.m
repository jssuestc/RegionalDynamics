function Step04_Balance(cfg)
% Step 04: rank-based balance. Independent of spatial coupling and spin.
P=cfg.Paths; [RAD,iRTD,ROI_IDs]=load_demo_metrics(cfg);
if ~isfolder(P.Balance), mkdir(P.Balance); end
RAD_Rank=tiedrank(RAD); iRTD_Rank=tiedrank(iRTD);
if numel(unique(RAD))==1 || numel(unique(iRTD))==1
    warning('Constant regional map: balance uses average tied ranks.');
end
BalanceAngle=atan2d(RAD_Rank,iRTD_Rank); BalanceDelta=BalanceAngle-45;
% Exact equality is defined by ranks, not floating-point angle rounding.
Predominance=repmat("Equal",numel(ROI_IDs),1);
Predominance(RAD_Rank>iRTD_Rank)="RAD";
Predominance(RAD_Rank<iRTD_Rank)="iRTD";
MeanBalance=mean(BalanceAngle); FractionRAD=mean(RAD_Rank>iRTD_Rank);
FractioniRTD=mean(RAD_Rank<iRTD_Rank); FractionEqual=mean(RAD_Rank==iRTD_Rank);
SubjectID=string(cfg.SubjectID);
T=table(ROI_IDs,RAD,iRTD,RAD_Rank,iRTD_Rank,BalanceAngle,BalanceDelta,Predominance, ...
    'VariableNames',{'ROI_ID','RAD','iRTD','RAD_Rank','iRTD_Rank','BalanceAngle_deg','BalanceDelta45_deg','Predominance'});
writetable(T,fullfile(P.Balance,'Balance_ROI.csv'));
writetable(table(SubjectID,MeanBalance,FractionRAD,FractioniRTD,FractionEqual), ...
    fullfile(P.Balance,'Balance_subject_summary.csv'));
save(fullfile(P.Balance,'RAD_iRTD_Balance_results.mat'),'SubjectID','ROI_IDs','RAD','iRTD', ...
    'RAD_Rank','iRTD_Rank','BalanceAngle','BalanceDelta','Predominance','MeanBalance', ...
    'FractionRAD','FractioniRTD','FractionEqual','cfg');
VA=spm_vol(cfg.AtlasFile); assert(numel(VA)==1,'Atlas must be 3D.'); atlas=spm_read_vols(VA);
write_atlas_map(VA,atlas,ROI_IDs,BalanceAngle,fullfile(P.Balance,'BalanceAngle_map.nii'),'Balance angle (degrees)');
write_atlas_map(VA,atlas,ROI_IDs,BalanceDelta,fullfile(P.Balance,'BalanceDelta45_map.nii'),'Balance angle minus 45 degrees');
f=figure('Color','w','Visible','off');
scatter(iRTD_Rank,RAD_Rank,28,BalanceDelta,'filled'); hold on;
plot([0 211],[0 211],'k--'); axis equal; xlim([0 211]); ylim([0 211]); box off;
xlabel('Ranked iRTD'); ylabel('Ranked RAD'); title('Within-participant rank-based balance');
% Blue: iRTD predominance; white: equal; orange: RAD predominance.
cmap=[interp1([0 1],[0.16 0.38 0.70;1 1 1],linspace(0,1,128)); ...
    interp1([0 1],[1 1 1;0.83 0.34 0.13],linspace(0,1,128))];
colormap(f,cmap); caxis([-45 45]); cb=colorbar; cb.Label.String='Balance angle - 45 (degrees)';
exportgraphics(f,fullfile(P.Balance,'Balance_rankspace.png'),'Resolution',300); close(f);

end
