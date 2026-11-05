function run_within_subject_peb_multiple(subjects, revisionRoot, spmDir)
% RUN_WITHIN_SUBJECT_PEB  Perform within-subject PEB across runs for NARPS subject(s)
%   subjects: Can be:
%             - Empty (default): processes the 48 analyzed subjects
%             - Single subject string: 'sub-001'
%             - Cell array of subjects: {'sub-001', 'sub-002', 'sub-003'}
%             - Numeric range: [1 125] or 1:125 (will find sub-001 through sub-125)
%
%   Performs *within-subject* PEB across runs (1–4) for each subject.
%
%   Saves corrected within-subject PEB results into:
%     revisionRoot/04_within_subject_PEB/PEB_corrected_sub-XXX_runwise.mat
%
%   Assumes corrected run-wise DCMs exist as:
%     revisionRoot/01_run_level_DCM_files/sub-XXX/DCM_corrected_run01.mat
%     revisionRoot/01_run_level_DCM_files/sub-XXX/DCM_corrected_run02.mat
%     revisionRoot/01_run_level_DCM_files/sub-XXX/DCM_corrected_run03.mat
%     revisionRoot/01_run_level_DCM_files/sub-XXX/DCM_corrected_run04.mat
%
%   Examples:
%     run_within_subject_peb_multiple()                       % All 48 analyzed subjects
%     run_within_subject_peb_multiple('sub-001')              % Single subject
%     run_within_subject_peb_multiple({'sub-001','sub-002'})  % Selected subjects
%     run_within_subject_peb_multiple(1:10)                   % Existing subjects in range

    clc;

    if nargin < 2 || isempty(revisionRoot); revisionRoot = getenv('NARPS_REVISION_ROOT'); end
    if nargin < 3 || isempty(spmDir); spmDir = getenv('SPM_DIR'); end
    if isempty(revisionRoot) || isempty(spmDir)
        error('Provide revisionRoot and spmDir, or set NARPS_REVISION_ROOT and SPM_DIR.');
    end
    if exist(fullfile(spmDir, 'spm.m'), 'file') ~= 2
        error('SPM25 was not found at: %s', spmDir);
    end
    addpath(spmDir);
    spm('Defaults', 'fMRI');
    spm_jobman('initcfg');

    %---------------------------------------------------------
    % Input and output directories
    %---------------------------------------------------------
    correctedDcmRoot = fullfile(revisionRoot, '01_run_level_DCM_files');
    correctedPebDir = fullfile(revisionRoot, '04_within_subject_PEB');
    
    %---------------------------------------------------------
    % Handle input arguments
    %---------------------------------------------------------
    if nargin == 0 || isempty(subjects)
        % Earlier default retained for reference.
        % subj_idx = [57];

        % All 48 subjects included in the corrected analysis.
        subj_idx = [ ...
             1  2  3  4  5  6  8  9 10 11 13 14 15 16 17 18 19 ...
            22 24 25 26 27 29 30 32 33 35 36 37 38 39 40 41 ...
            43 44 45 46 47 49 50 51 52 53 54 55 56 57 58];
        subjects = arrayfun(@(x) sprintf('sub-%03d', x), subj_idx, 'UniformOutput', false);
        fprintf('Using all 48 analyzed subjects.\n');
        
    elseif isnumeric(subjects)
        % Range of subject numbers - find which ones exist
        fprintf('Searching for subjects in range %d to %d...\n', min(subjects), max(subjects));
        subjectList = {};
        for subNum = subjects
            subName = sprintf('sub-%03d', subNum);
            % Earlier input directory and filename retained for reference.
            % subDir = fullfile(firstLvlRoot, subName);
            % dcm01 = fullfile(subDir, 'DCM_run01.mat');

            subDir = fullfile(correctedDcmRoot, subName);
            if exist(subDir, 'dir')
                % Check whether the first corrected DCM exists.
                dcm01 = fullfile(subDir, 'DCM_corrected_run01.mat');
                if exist(dcm01, 'file')
                    subjectList{end+1} = subName;
                end
            end
        end
        fprintf('Found %d subjects with DCM results in the specified range.\n', numel(subjectList));
        subjects = subjectList;
        
    elseif ischar(subjects) || isstring(subjects)
        % Single subject string
        subjects = {char(subjects)};
        
    elseif ~iscell(subjects)
        error('Input must be a subject string, cell array of subjects, numeric range, or empty for default');
    end
    
    if isempty(subjects)
        fprintf('No subjects found to process.\n');
        return;
    end
    
    fprintf('Running within-subject PEB for %d subjects\n\n', numel(subjects));
    
    %---------------------------------------------------------
    % Create corrected within-subject PEB output directory
    %---------------------------------------------------------
    % Earlier code retained for reference.
    % secondLvlDir = fullfile(rootDir, 'secondlevel');

    if ~exist(correctedPebDir, 'dir')
        mkdir(correctedPebDir);
    end
    
    %---------------------------------------------------------
    % Process each subject
    %---------------------------------------------------------
    successCount = 0;
    errorCount = 0;
    
    for s = 1:numel(subjects)
        sub_id = char(subjects{s});
        
        fprintf('============================================\n');
        fprintf('  Subject %s (%d of %d)\n', sub_id, s, numel(subjects));
        fprintf('============================================\n');
        
        try
            % Earlier call retained for reference.
            % process_single_subject_peb(sub_id, rootDir, secondLvlDir);
            process_single_subject_peb(sub_id, correctedDcmRoot, correctedPebDir);
            fprintf('SUCCESS: Completed PEB for %s\n\n', sub_id);
            successCount = successCount + 1;
        catch ME
            fprintf('ERROR processing PEB for %s: %s\n', sub_id, ME.message);
            fprintf('Skipping to next subject...\n\n');
            errorCount = errorCount + 1;
            continue;
        end
    end
    
    fprintf('============================================\n');
    fprintf('Within-subject PEB analyses completed.\n');
    fprintf('Successfully processed: %d subjects\n', successCount);
    fprintf('Errors encountered: %d subjects\n', errorCount);
    fprintf('============================================\n\n');
end

function process_single_subject_peb(sub_id, correctedDcmRoot, correctedPebDir)
% Process within-subject PEB for a single subject

    % Earlier input directory retained for reference.
    % subjFirstLvlDir = fullfile(rootDir, 'firstlevel', sub_id);
    subjFirstLvlDir = fullfile(correctedDcmRoot, sub_id);
    
    if ~exist(subjFirstLvlDir, 'dir')
        error('First-level directory NOT found: %s', subjFirstLvlDir);
    end
    
    %---------------------------------------------------------
    % Collect corrected run-wise DCMs.
    %---------------------------------------------------------
    DCM_files = cell(4,1);
    for r = 1:4
        % Earlier filename retained for reference.
        % thisDCM = fullfile(subjFirstLvlDir, sprintf('DCM_run%02d.mat', r));
        thisDCM = fullfile(subjFirstLvlDir, sprintf('DCM_corrected_run%02d.mat', r));
        if ~exist(thisDCM, 'file')
            error('Missing DCM file for run %02d: %s', r, thisDCM);
        end

        dcmCheck = load(thisDCM, 'DCM');
        if ~isfield(dcmCheck, 'DCM') || ...
                ~isfield(dcmCheck.DCM, 'Ep') || ...
                ~isfield(dcmCheck.DCM, 'Cp') || ...
                ~isfield(dcmCheck.DCM, 'F') || ...
                ~isnumeric(dcmCheck.DCM.Cp) || ...
                any(~isfinite(dcmCheck.DCM.Cp(:))) || ...
                ~isscalar(dcmCheck.DCM.F) || ...
                ~isfinite(dcmCheck.DCM.F)
            error('Corrected DCM was not successfully inverted: %s', thisDCM);
        end
        DCM_files{r} = thisDCM;
    end
    
    fprintf('  Using DCMs:\n');
    for r = 1:4
        fprintf('    %s\n', DCM_files{r});
    end
    
    %---------------------------------------------------------
    % Specify PEB model across runs (within subject)
    %---------------------------------------------------------
    M         = [];
    M.X       = ones(numel(DCM_files), 1);   % intercept-only (mean across runs)
    M.Q       = 'all';
    M.maxit   = 128;
    
    fprintf('  Estimating PEB for %s ...\n', sub_id);
    
    %---------------------------------------------------------
    % Run PEB
    %---------------------------------------------------------
    [PEB, DCM_peb] = spm_dcm_peb(DCM_files, M);
    
    %---------------------------------------------------------
    % Save corrected PEB results without overwriting the original analysis.
    %---------------------------------------------------------
    % Earlier output filename retained for reference.
    % outFile = fullfile(secondLvlDir, sprintf('PEB_%s_runwise.mat', sub_id));
    outFile = fullfile(correctedPebDir, sprintf('PEB_corrected_%s_runwise.mat', sub_id));
    if exist(outFile, 'file')
        error(['Corrected PEB output already exists and will not be overwritten: %s. ' ...
            'Move or rename it explicitly before rerunning this subject.'], outFile);
    end
    save(outFile, 'PEB', 'DCM_peb', '-v6');
    fprintf('  Saved PEB results to:\n    %s\n', outFile);
end
