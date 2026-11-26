function [RAD,iRTD,ROI_IDs] = load_demo_metrics(cfg)
S=load(cfg.Paths.MetricsMat,'RAD','iRTD','ROI_IDs');
ROI_IDs=double(S.ROI_IDs(:)); RAD=double(S.RAD(:)); iRTD=double(S.iRTD(:));
assert(isequal(ROI_IDs,(1:210)') && numel(RAD)==210 && numel(iRTD)==210, ...
    'Expected 210 cortical values in BN atlas order.');
assert(all(isfinite(RAD) & isfinite(iRTD)), 'Nonfinite cortical values; inspect RAD_iRTD_ROI.csv.');
end
