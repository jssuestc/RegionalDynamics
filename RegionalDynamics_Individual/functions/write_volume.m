function write_volume(V,Y,outFile,description)
% Construct a fresh header to avoid carrying over integer scaling.
H=struct('fname',outFile,'dim',V.dim(1:3),'dt',[16 0], ...
    'mat',V.mat,'pinfo',[1;0;0],'n',[1 1],'descrip',description);
spm_write_vol(H,double(Y));
end
