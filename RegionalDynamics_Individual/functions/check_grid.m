function check_grid(Vref,V,label)
% Check dimensions and affine, not merely the number of voxels.
for k=1:numel(V)
    assert(isequal(Vref.dim(1:3),V(k).dim(1:3)), ...
        '%s dimensions differ from reference. Reslice explicitly before running.',label);
    assert(all(isfinite(V(k).mat(:))) && max(abs(Vref.mat(:)-V(k).mat(:)))<1e-4, ...
        '%s affine differs from reference. Check alignment; no automatic reslicing.',label);
end
end
