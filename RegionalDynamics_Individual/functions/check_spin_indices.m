function check_spin_indices(I,ROI_IDs,nSpin)
assert(isequal(ROI_IDs(:),(1:210)') && isequal(size(I),[210,nSpin]), ...
    'Spin matrix must be 210 x NSpin in BN atlas ID order.');
assert(all(isfinite(I(:)) & I(:)==fix(I(:)) & I(:)>=1 & I(:)<=210), ...
    'Spin indices must be integer ROI IDs 1:210.');
assert(all(all(sort(I,1)==(1:210)')),'Each spin must be a one-to-one permutation.');
assert(all(all(mod(I,2)==mod((1:210)',2))),'Spin crosses hemispheres.');
assert(~any(all(I==(1:210)',1)),'Identity spin detected.');
end
