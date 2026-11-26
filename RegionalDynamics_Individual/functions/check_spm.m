function check_spm(cfg)
if ~isempty(cfg.SPMDir)
    assert(isfolder(cfg.SPMDir),'SPMDir does not exist: %s',cfg.SPMDir);
    addpath(cfg.SPMDir);
end
assert(exist('spm_vol','file')==2 && exist('spm_read_vols','file')==2 && ...
    exist('spm_write_vol','file')==2,'SPM12 is required. Set cfg.SPMDir in Run_Individual.m.');
end
