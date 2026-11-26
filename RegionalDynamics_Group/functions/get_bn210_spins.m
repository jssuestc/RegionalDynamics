function S = get_bn210_spins(spinFile,atlasDir,sphereDir)
% Reuse a validated BN210 rotation matrix or generate it from actual spheres.
% Same bilateral Vasa-style matching and seed as Step03_SpatialCoupling.
nSpin = 5000; seed = 1;
if isfile(spinFile)
    S = load(spinFile);

else
    [SpinIndices,ROI_IDs,Hemisphere,ROIName,Centroids] = ...
        BN210_spin_indices(atlasDir,sphereDir,nSpin,seed);
    SpinParameters = struct('NSpin',nSpin,'Seed',seed, ...
        'SphereDir',sphereDir,'AnnotationDir',atlasDir, ...
        'Algorithm','Vasa one-to-one bilateral cortical spins');
    S = struct('SpinIndices',SpinIndices,'ROI_IDs',ROI_IDs, ...
        'Hemisphere',Hemisphere,'ROIName',ROIName,'Centroids',Centroids, ...
        'SpinParameters',SpinParameters);
    folder = fileparts(spinFile);
    if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
    save(spinFile,'-struct','S');
end
S.SpinIndices = double(S.SpinIndices);
check_spin_indices(S.SpinIndices,double(S.ROI_IDs),nSpin);
end
