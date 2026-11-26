function [SpinIndices,ROI_IDs,Hemisphere,ROIName,Centroids] = BN210_spin_indices(atlasDir,sphereDir,nPerm,seed)
% Adapt CorMap_spin / Schaefer_Sphere_Centroids supplied by the author.
% Actual fsaverage BN210 annotation color table indices are atlas ID+1.
% rotate_parcellation returns LH-first/RH-second *row positions*.
% Translate both its target rows and source indices to ROI-ID order 1:210.
ROI_IDs=(1:210)'; LH=(1:2:209)'; RH=(2:2:210)';
Hemisphere=repmat("R",210,1); Hemisphere(LH)="L";
[coordL,nameL]=centroids('lh',LH,atlasDir,sphereDir);
[coordR,nameR]=centroids('rh',RH,atlasDir,sphereDir);
Centroids=nan(210,3); Centroids(LH,:)=coordL; Centroids(RH,:)=coordR;
ROIName=strings(210,1); ROIName(LH)=nameL; ROIName(RH)=nameR;
oldRng=rng; restoreRng=onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(seed,'twister');
P=rotate_parcellation(coordL,coordR,nPerm);
order=[LH;RH];
assert(isequal(size(P),[210,nPerm]) && all(isfinite(P(:)) & P(:)==fix(P(:)) & P(:)>=1 & P(:)<=210), ...
    'Unexpected rotate_parcellation output.');
SpinIndices=zeros(210,nPerm);
SpinIndices(order,:)=order(P);
check_spin_indices(SpinIndices,ROI_IDs,nPerm);
assert(all(all(Hemisphere(SpinIndices)==Hemisphere)),'Cross-hemisphere spin assignments detected.');
end

function [C,names]=centroids(hemi,ids,atlasDir,sphereDir)
A=fullfile(atlasDir,[hemi '.BN_Atlas.annot']);
L=fullfile(atlasDir,[hemi '.cortex.label']);
F=fullfile(sphereDir,[hemi '.sphere.reg']);
assert(isfile(F),'Missing actual registered sphere: %s. Set the sphere directory in the main script to fsaverage/surf.',F);
xyz=double(read_surf(F));
[v,lab,ct]=read_annotation(A);
assert(isequal(size(xyz),[163842,3]) && numel(v)==163842 && ...
    isequal(sort(v(:)),(0:163841)') && numel(lab)==numel(v), ...
    'Sphere and BN annotation must be matched full-resolution fsaverage meshes.');
r=sqrt(sum(xyz.^2,2));
assert(all(isfinite(r) & r>0) && max(abs(r-median(r)))/median(r)<0.05, ...
    'Input surface does not appear spherical: use sphere.reg, not white/inflated/pial or cortex.label XYZ.');
assert(ct.numEntries==211 && strcmpi(ct.struct_names{1},'Unknown'), ...
    'Unexpected BN210 annotation color table.');
% Use vertex IDs, not assumed annotation row order. Cortex.label is only
% used for membership (its anatomical XYZ are NOT spherical coordinates).
fid=fopen(L,'rt'); assert(fid>=0,'Cannot open %s',L);
cleanup=onCleanup(@() fclose(fid));
fgetl(fid); nv=str2double(strtrim(fgetl(fid)));
D=textscan(fid,'%f%f%f%f%f','CollectOutput',true); D=D{1};
assert(isfinite(nv) && size(D,1)==nv && size(D,2)==5,'Unexpected cortex.label format.');
li=D(:,1); assert(all(li>=0 & li<163842 & li==fix(li)) && numel(unique(li))==nv,'Invalid cortex vertex IDs.');
inCortex=false(163842,1); inCortex(li+1)=true;
vertexLabels=zeros(163842,1); vertexLabels(v+1)=lab;
% Preserve original CorMap_spin definition: arithmetic mean of sphere.reg
% coordinates for each parcel, with NO extra radial projection.
C=nan(numel(ids),3); names=strings(numel(ids),1);
for j=1:numel(ids)
    entry=ids(j)+1; names(j)=string(ct.struct_names{entry});
    suffix=upper(hemi(1));
    assert(endsWith(names(j),['_' suffix]),'ROI %d name/hemisphere mismatch.',ids(j));
    code=ct.table(entry,5); idx=vertexLabels==code;
    assert(code~=0 && any(idx) && all(inCortex(idx)), ...
        'ROI %d missing or outside cortex.label.',ids(j));
    C(j,:)=mean(xyz(idx,:),1);
end
assert(all(isfinite(C(:))) && all(sqrt(sum(C.^2,2))>0),'Invalid centroids.');
end
