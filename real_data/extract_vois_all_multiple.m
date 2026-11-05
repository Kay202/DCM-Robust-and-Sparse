function extract_vois_all_multiple(subjects, rootDir)
% EXTRACT_VOIS_ALL_RUNS  Extract VOIs for all runs for NARPS subject(s)
%   subjects: Can be:
%             - Single subject string: 'sub-001'
%             - Cell array of subjects: {'sub-001', 'sub-002', 'sub-003'}
%             - Numeric range: [1 125] or 1:125 (will find sub-001 through sub-125)
%
%   Examples:
%     extract_vois_all_multiple('sub-001', rootDir)
%     extract_vois_all_multiple({'sub-001', 'sub-002'}, rootDir)

    if nargin < 2 || isempty(rootDir)
        rootDir = getenv('NARPS_ROOT');
    end
    if isempty(rootDir)
        error('Provide rootDir or set the NARPS_ROOT environment variable.');
    end
    firstLvlRoot = fullfile(rootDir, 'firstlevel');
    
    % Handle different input types
    if isnumeric(subjects)
        % Range of subject numbers - find which ones exist
        fprintf('Searching for subjects in range %d to %d...\n', min(subjects), max(subjects));
        subjectList = {};
        for subNum = subjects
            subName = sprintf('sub-%03d', subNum);
            subDir = fullfile(firstLvlRoot, subName);
            if exist(subDir, 'dir')
                spmFile = fullfile(subDir, 'SPM.mat');
                if exist(spmFile, 'file')
                    subjectList{end+1} = subName;
                end
            end
        end
        fprintf('Found %d subjects with SPM.mat in the specified range.\n', numel(subjectList));
        subjects = subjectList;
        
    elseif ischar(subjects) || isstring(subjects)
        % Single subject string
        subjects = {char(subjects)};
        
    elseif ~iscell(subjects)
        error('Input must be a subject string, cell array of subjects, or numeric range');
    end
    
    if isempty(subjects)
        fprintf('No subjects found to process.\n');
        return;
    end
    
    % SPM setup
    spm('Defaults','fMRI');
    spm_jobman('initcfg');
    
    % Process each subject
    successCount = 0;
    errorCount = 0;
    
    for subIdx = 1:numel(subjects)
        sub = char(subjects{subIdx});
        fprintf('\n========================================\n');
        fprintf('Extracting VOIs for subject: %s (%d/%d)\n', sub, subIdx, numel(subjects));
        fprintf('========================================\n');
        
        try
            process_single_subject_voi(sub, rootDir);
            fprintf('SUCCESS: Completed VOI extraction for %s\n', sub);
            successCount = successCount + 1;
        catch ME
            fprintf('ERROR extracting VOIs for %s: %s\n', sub, ME.message);
            fprintf('Skipping to next subject...\n');
            errorCount = errorCount + 1;
            continue;
        end
    end
    
    fprintf('\n========================================\n');
    fprintf('VOI extraction complete!\n');
    fprintf('Successfully processed: %d subjects\n', successCount);
    fprintf('Errors encountered: %d subjects\n', errorCount);
    fprintf('========================================\n');
end

function process_single_subject_voi(sub, rootDir)
% Extract VOIs for a single subject - original function logic

    firstLvlDir = fullfile(rootDir, 'firstlevel', sub);
    roiDir      = fullfile(rootDir, 'firstlevel', 'ROIs');
    
    SPMmat = fullfile(firstLvlDir,'SPM.mat');
    
    if ~exist(SPMmat, 'file')
        error('SPM.mat not found: %s', SPMmat);
    end
    
    load(SPMmat,'SPM');
    nSess = numel(SPM.Sess);   % number of runs in the GLM
    
    % Define ROI files
    roiFiles = {
        fullfile(roiDir,'vmPFC_mask.nii')
        fullfile(roiDir,'vStriatum_mask.nii')
        fullfile(roiDir,'Amygdala_mask.nii')
        fullfile(roiDir,'aIns_mask.nii')
    };
    
    roiBaseNames = {'vmPFC','vStr','Amygd','aIns'};
    
    % Check if ROI files exist
    for k = 1:numel(roiFiles)
        if ~exist(roiFiles{k}, 'file')
            error('ROI mask file not found: %s', roiFiles{k});
        end
    end
    
    fprintf('  Extracting VOIs for %d sessions...\n', nSess);
    
    % Extract VOIs for each session and ROI
    for s = 1:nSess
        fprintf('  Session %d/%d:\n', s, nSess);
        
        for k = 1:numel(roiFiles)
            clear matlabbatch
            voiName = sprintf('%s_run%d', roiBaseNames{k}, s);
            
            fprintf('    - %s\n', voiName);
            
            matlabbatch{1}.spm.util.voi.spmmat  = {SPMmat};
            matlabbatch{1}.spm.util.voi.adjust  = 0;      % no adjustment
            matlabbatch{1}.spm.util.voi.session = s;      % <- loop over sessions
            matlabbatch{1}.spm.util.voi.name    = voiName;
            matlabbatch{1}.spm.util.voi.roi{1}.mask.image    = {[roiFiles{k} ',1']};
            matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
            matlabbatch{1}.spm.util.voi.expression           = 'i1';
            
            spm_jobman('run', matlabbatch);
        end
    end
    
    fprintf('  All VOIs extracted successfully for %s\n', sub);
end
