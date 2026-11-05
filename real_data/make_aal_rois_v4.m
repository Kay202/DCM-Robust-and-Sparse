function make_aal_rois_v4(aalPath, outDir)
    if nargin < 1 || isempty(aalPath)
        aalPath = getenv('AAL_ATLAS');
    end
    if nargin < 2 || isempty(outDir)
        narpsRoot = getenv('NARPS_ROOT');
        if ~isempty(narpsRoot)
            outDir = fullfile(narpsRoot, 'firstlevel', 'ROIs');
        end
    end
    if isempty(aalPath) || isempty(outDir)
        error('Provide aalPath and outDir, or set AAL_ATLAS and NARPS_ROOT.');
    end

    if ~exist(outDir,'dir'); mkdir(outDir); end

    % Load AAL atlas
    V = spm_vol(aalPath);
    [Y, XYZ] = spm_read_vols(V);

    % ROI index definitions derived from ROI_MNI_V4.txt
    roiDefs = {
        'vmPFC',      [2611 2612 2601 2602]
        'vStriatum',  [7001 7002 7011 7012]
        'Amygdala',   [4201 4202]
        'aIns',       [3001 3002]
    };

    for i = 1:size(roiDefs,1)
        name = roiDefs{i,1};
        idx  = roiDefs{i,2};

        mask = ismember(Y, idx);
        Vm = V;
        Vm.fname   = fullfile(outDir, [name '_mask.nii']);
        Vm.descrip = ['AAL V4 ROI mask: ' name];

        spm_write_vol(Vm, mask);
        fprintf('Wrote ROI mask: %s\n', Vm.fname);
    end
end
