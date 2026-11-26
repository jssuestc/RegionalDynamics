function write_atlas_map(VA,atlas,ids,values,outFile,description)
% Float32 3D NIfTI, original atlas affine; zero outside selected labels.
assert(numel(ids)==numel(values),'Map length does not match ROI IDs.');
Y=zeros(size(atlas));
for j=1:numel(ids), Y(atlas==ids(j))=values(j); end
write_volume(VA,Y,outFile,description);
end
