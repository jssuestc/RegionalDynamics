function [X,names] = load_reference_maps(inputFile,nMaps)
% Reference maps are stored as map x ROI, or ROI x map.
% Never truncate an unexpected matrix to force it to match the labels.
S = load(inputFile);
assert(isfield(S,'roi_matrix') && isfield(S,'file_names'), ...
    'Expected roi_matrix and file_names in %s.',inputFile);
A = S.roi_matrix;
assert(isnumeric(A) && isreal(A) && ismatrix(A),'Invalid roi_matrix: %s',inputFile);
if isequal(size(A),[nMaps,210])
    X = double(A');
elseif isequal(size(A),[210,nMaps])
    X = double(A);
else
    error('Expected %d x 210 or 210 x %d roi_matrix; found %d x %d in %s. Check the saved data.', ...
        nMaps,nMaps,size(A,1),size(A,2),inputFile);
end
names = string(S.file_names(:));
assert(numel(names)==nMaps && all(~ismissing(names)) && ...
    all(strlength(strtrim(names))>0) && numel(unique(names))==nMaps, ...
    'Expected %d distinct map names in %s.',nMaps,inputFile);
assert(all(isfinite(X(:))) && all(std(X,0,1)>0), ...
    'Reference maps must be finite and nonconstant: %s',inputFile);
if isfield(S,'roi_labels')
    assert(isequal(double(S.roi_labels(:)),(1:210)'), ...
        'ROI labels must be in BN cortical order 1:210: %s',inputFile);
end
end
