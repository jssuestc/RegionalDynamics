function y = load_mean_map(inputFile)
% Read a single BN210 mean cortical map, retaining the original ROI order.
S = load(inputFile);
assert(isfield(S,'dataROI') && isnumeric(S.dataROI) && ...
    isreal(S.dataROI) && isvector(S.dataROI) && numel(S.dataROI)==210, ...
    'Expected a 210-element dataROI vector in %s.',inputFile);
y = double(S.dataROI(:));
assert(all(isfinite(y)) && std(y)>0, ...
    'Mean map must be finite and nonconstant: %s',inputFile);
if isfield(S,'roi_labels')
    assert(isequal(double(S.roi_labels(:)),(1:210)'), ...
        'ROI labels must be in BN cortical order 1:210: %s',inputFile);
end
end
