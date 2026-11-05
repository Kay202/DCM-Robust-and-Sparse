function first_level_glm_multiple(subjects, rootDir)
% FIRST_LEVEL_GLM  Build and estimate 1st-level GLM for NARPS subject(s)
%   subjects: Can be:
%             - Single subject string: 'sub-001'
%             - Cell array of subjects: {'sub-001', 'sub-002', 'sub-003'}
%             - Numeric range: [1 125] or 1:125 (will find sub-001 through sub-125)
%
%   Examples:
%     first_level_glm_multiple('sub-001', rootDir)
%     first_level_glm_multiple({'sub-001', 'sub-002'}, rootDir)

    if nargin < 2 || isempty(rootDir)
        rootDir = getenv('NARPS_ROOT');
    end
    if isempty(rootDir)
        error('Provide rootDir or set the NARPS_ROOT environment variable.');
    end
    derivRoot = fullfile(rootDir, 'derivatives','fmriprep');
    
    % Handle different input types
    if isnumeric(subjects)
        % Range of subject numbers - find which ones exist
        fprintf('Searching for subjects in range %d to %d...\n', min(subjects), max(subjects));
        subjectList = {};
        for subNum = subjects
            subName = sprintf('sub-%03d', subNum);
            subDir = fullfile(derivRoot, subName);
            if exist(subDir, 'dir')
                subjectList{end+1} = subName;
            end
        end
        fprintf('Found %d subjects in the specified range.\n', numel(subjectList));
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
    
    % Process each subject
    successCount = 0;
    errorCount = 0;
    
    for subIdx = 1:numel(subjects)
        sub = char(subjects{subIdx});
        fprintf('\n========================================\n');
        fprintf('Processing subject: %s (%d/%d)\n', sub, subIdx, numel(subjects));
        fprintf('========================================\n');
        
        try
            process_single_subject(sub, rootDir);
            fprintf('SUCCESS: Completed %s\n', sub);
            successCount = successCount + 1;
        catch ME
            fprintf('ERROR processing %s: %s\n', sub, ME.message);
            fprintf('Skipping to next subject...\n');
            errorCount = errorCount + 1;
            continue;
        end
    end
    
    fprintf('\n========================================\n');
    fprintf('Processing complete!\n');
    fprintf('Successfully processed: %d subjects\n', successCount);
    fprintf('Errors encountered: %d subjects\n', errorCount);
    fprintf('========================================\n');
end

function process_single_subject(sub, rootDir)
% Process a single subject - original function logic

    bidsRoot    = rootDir;                                     % where sub-XXX\func\*_events.tsv live
    derivRoot   = fullfile(rootDir, 'derivatives','fmriprep'); % where preproc BOLD + confounds live
    firstLvlRoot = fullfile(rootDir, 'firstlevel');
    % -------------------------------------------------------------------

    if ~exist(firstLvlRoot,'dir'); mkdir(firstLvlRoot); end

    nRuns      = 4;
    subDirFunc = fullfile(derivRoot, sub, 'func');
    outDir     = fullfile(firstLvlRoot, sub);
    if ~exist(outDir,'dir'); mkdir(outDir); end

    % ---------- Build task onsets/durations ----------
    [onsets, durs] = build_conditions_for_subject(sub, bidsRoot, nRuns);

    % ---------- SPM setup ----------
    spm('Defaults','fMRI');
    spm_jobman('initcfg');

    clear matlabbatch
    matlabbatch{1}.spm.stats.fmri_spec.dir = {outDir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units  = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT     = 1;    % TR = 1s (NARPS)
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    for r = 1:nRuns

        %% ---------- Locate preprocessed BOLD ----------
        patBold = sprintf('%s_task-MGT_run-%02d_bold_space-MNI152NLin2009cAsym_preproc.nii*', sub, r);
        bf = dir(fullfile(subDirFunc, patBold));
        if isempty(bf)
            error('Missing BOLD file for %s run %02d in %s', sub, r, subDirFunc);
        end
        boldPath = fullfile(bf(1).folder, bf(1).name);

        % gunzip if needed
        [pth, nm, ext] = fileparts(boldPath);
        if strcmpi(ext,'.gz')
            gunzip(boldPath);
            boldPath = fullfile(pth, nm);  % now *.nii
        end

        V = spm_vol(boldPath);
        nScan = numel(V);
        scans = cell(nScan,1);
        for i = 1:nScan
            scans{i} = sprintf('%s,%d', boldPath, i);
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = scans;

        %% ---------- Task regressors (GainDom, LossDom, Middle) ----------
        % Gain-dominant trials
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).name     = 'GainDom';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).onset    = onsets.GainDom{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).duration = durs.GainDom{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).tmod     = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(1).pmod     = struct('name',{},'param',{},'poly',{});

        % Loss-dominant trials
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).name     = 'LossDom';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).onset    = onsets.LossDom{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).duration = durs.LossDom{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).tmod     = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(2).pmod     = struct('name',{},'param',{},'poly',{});

        % Middle (ambiguous) trials
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(3).name     = 'Middle';
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(3).onset    = onsets.Middle{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(3).duration = durs.Middle{r};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(3).tmod     = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(3).pmod     = struct('name',{},'param',{},'poly',{});

        %% ---------- Confounds (Option A) ----------
        % Motion (6), FD, aCompCor, Cosines, NonSteadyStateOutliers

        patConf = sprintf('%s_task-MGT_run-%02d_bold_confounds*.tsv', sub, r);
        cf = dir(fullfile(subDirFunc, patConf));
        if isempty(cf)
            error('Missing confounds file for %s run %02d in %s', sub, r, subDirFunc);
        end
        confFile = fullfile(cf(1).folder, cf(1).name);

        C = readtable(confFile, 'FileType','text', 'VariableNamingRule','preserve');
        varNames = C.Properties.VariableNames;

        Xconf = [];

        % --- Motion parameters ---
        % Try lowercase (fMRIPrep-style) then uppercase (like your example)
        motionSets = { ...
            {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'}, ...
            {'X','Y','Z','RotX','RotY','RotZ'}};

        motionUsed = false;
        for sIdx = 1:numel(motionSets)
            thisSet = motionSets{sIdx};
            if all(ismember(thisSet, varNames))
                for k = 1:numel(thisSet)
                    Xconf = [Xconf, C.(thisSet{k})];
                end
                motionUsed = true;
                break;
            end
        end

        % --- Motion derivatives ---
        if motionUsed && ~isempty(Xconf)
            Xconf = [Xconf, [zeros(1,size(Xconf,2)); diff(Xconf)]];
        end

        % --- Framewise Displacement ---
        fdCandidates = {'framewise_displacement','FramewiseDisplacement'};
        for k = 1:numel(fdCandidates)
            if ismember(fdCandidates{k}, varNames)
                Xconf = [Xconf, C.(fdCandidates{k})];
                break;
            end
        end

                % --- aCompCor components ---
        for k = 1:numel(varNames)
            nm = varNames{k};
            if startsWith(nm,'aCompCor','IgnoreCase',true) || ...
               startsWith(nm,'a_comp_cor','IgnoreCase',true)
                Xconf = [Xconf, C.(nm)];
            end
        end

        % --- Non-steady-state outliers ---
        for k = 1:numel(varNames)
            nm = varNames{k};
            if contains(nm,'NonSteadyStateOutlier','IgnoreCase',true)
                Xconf = [Xconf, C.(nm)];
            end
        end

        % ---------- CLEAN NaNs / Infs ----------
        if ~isempty(Xconf)
            % Replace any NaN or Inf with 0
            bad = ~isfinite(Xconf);   % true where NaN or Inf
            if any(bad(:))
                fprintf('Run %02d: replacing %d NaN/Inf entries in confounds with 0.\n', ...
                        r, nnz(bad));
                Xconf(bad) = 0;
            end
        end

        % ---------- Match #rows of confounds to #scans ----------
        if isempty(Xconf)
            Xconf = zeros(nScan,0);
        else
            nRows = size(Xconf,1);
            if nRows > nScan
                Xconf = Xconf(1:nScan,:);
            elseif nRows < nScan
                Xconf = [Xconf; zeros(nScan-nRows,size(Xconf,2))];
            end
        end


        confMatFile = fullfile(outDir, sprintf('confounds_run%02d.txt', r));
        dlmwrite(confMatFile, Xconf, 'delimiter','\t');

        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = {confMatFile};
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf       = 128;
    end

    % HRF only, no temporal/dispersion derivatives
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt   = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';

    % ---------- Model estimation ----------
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = ...
        cfg_dep('fMRI model specification: SPM.mat File', ...
        substruct('.','val', '{}',{1},'.','val', '{}',{1},'.','val', '{}',{1}), ...
        substruct('.','spmmat'));

    spm_jobman('run', matlabbatch);

    % ---------- run model specification only ----------
    %spm_jobman('run', matlabbatch);

    % ---------- Model estimation (Classical OLS) ----------
    %clear estbatch
    %estbatch{1}.spm.stats.fmri_est.spmmat = {fullfile(outDir,'SPM.mat')};
    %estbatch{1}.spm.stats.fmri_est.method.Classical = 1;

    %spm_jobman('run', estbatch);
end
