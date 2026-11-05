function multiple_dcm_runs(subjects, rootDir, outputRoot, spmDir)
% CREATE_AND_INVERT_DCM_PER_RUN_NOUFUNCTION  Create and invert DCM for NARPS subject(s)
%   subjects: Can be:
%             - Single subject string: 'sub-001'
%             - Cell array of subjects: {'sub-001', 'sub-002', 'sub-003'}
%             - Numeric range: [1 125] or 1:125 (will find sub-001 through sub-125)
%
%   Examples:
%     multiple_dcm_runs('sub-001', rootDir, outputRoot, spmDir)
%     multiple_dcm_runs({'sub-001', 'sub-002'}, rootDir, outputRoot, spmDir)

    if nargin < 2 || isempty(rootDir); rootDir = getenv('NARPS_ROOT'); end
    if nargin < 3 || isempty(outputRoot)
        revisionRoot = getenv('NARPS_REVISION_ROOT');
        if ~isempty(revisionRoot); outputRoot = fullfile(revisionRoot, '01_run_level_DCM_files'); end
    end
    if nargin < 4 || isempty(spmDir); spmDir = getenv('SPM_DIR'); end
    if isempty(rootDir) || isempty(outputRoot) || isempty(spmDir)
        error(['Provide rootDir, outputRoot, and spmDir, or set NARPS_ROOT, ' ...
               'NARPS_REVISION_ROOT, and SPM_DIR.']);
    end
    firstLvlRoot = fullfile(rootDir, 'firstlevel');
    if exist(fullfile(spmDir, 'spm.m'), 'file') ~= 2
        error('SPM25 was not found at: %s', spmDir);
    end
    addpath(spmDir);
    spm('Defaults', 'fMRI');
    spm_jobman('initcfg');
    
    % Handle different input types
    if isnumeric(subjects)
        % Range of subject numbers - find which ones exist
        fprintf('Searching for subjects in range %d to %d...\n', min(subjects), max(subjects));
        subjectList = {};
        for subNum = subjects
            subName = sprintf('sub-%03d', subNum);
            subDir = fullfile(firstLvlRoot, subName);
            if exist(subDir, 'dir')
                subjectList{end+1} = subName;
            end
        end
        fprintf('Found %d subjects with first-level results in the specified range.\n', numel(subjectList));
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
        fprintf('Processing DCM for subject: %s (%d/%d)\n', sub, subIdx, numel(subjects));
        fprintf('========================================\n');
        
        try
            % Earlier call did not pass the input and output roots.
            % process_single_subject_dcm(sub);
            process_single_subject_dcm(sub, firstLvlRoot, outputRoot);
            fprintf('SUCCESS: Completed DCM for %s\n', sub);
            successCount = successCount + 1;
        catch ME
            fprintf('ERROR processing DCM for %s: %s\n', sub, ME.message);
            fprintf('Skipping to next subject...\n');
            errorCount = errorCount + 1;
            continue;
        end
    end
    
    fprintf('\n========================================\n');
    fprintf('DCM Processing complete!\n');
    fprintf('Successfully processed: %d subjects\n', successCount);
    fprintf('Errors encountered: %d subjects\n', errorCount);
    fprintf('========================================\n');
end

function process_single_subject_dcm(sub, firstLvlRoot, outputRoot)
% Process DCM for a single subject - original function logic

    %% =======================
    % PATHS
    % ========================
    firstLvl = fullfile(firstLvlRoot, sub);
    SPMmat   = fullfile(firstLvl, 'SPM.mat');
    outputDir = fullfile(outputRoot, sub);

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    if ~exist(SPMmat,'file')
        error('SPM.mat not found: %s', SPMmat);
    end

    % SPM is initialized in the main function before SPM.mat is loaded so
    % that MATLAB can reconstruct SPM's nifti and file_array objects.
    % spm('Defaults', 'fMRI');
    % spm_jobman('initcfg');
    load(SPMmat, 'SPM');

    %% Basic info
    RT    = SPM.xY.RT;
    nRuns = numel(SPM.Sess);

    %% ROIs used in DCM
    roiNames = {'vmPFC','vStr','Amygd','aIns'};
    nROI = numel(roiNames);

    %% Loop across sessions/runs
    for r = 1:nRuns

        fprintf('\n==============================\n');
        fprintf('  Creating DCM for %s | Run %02d\n', sub, r);
        fprintf('==============================\n');

        % --------------------------------------------------
        % Load VOIs for this run (only use y and X0)
        % --------------------------------------------------
        Ymat        = [];
        regionNames = roiNames;
        X0          = [];

        for k = 1:nROI
            voiFile = fullfile(firstLvl, ...
                sprintf('VOI_%s_run%d_%d.mat', roiNames{k}, r, r));

            if ~exist(voiFile,'file')
                error('Missing VOI file for %s: %s', roiNames{k}, voiFile);
            end

            tmp = load(voiFile);  % contains xY

            if k == 1
                % determine number of scans from first VOI
                T    = size(tmp.xY.y,1);
                Ymat = zeros(T, nROI);
                % take confounds from first VOI if present
                if isfield(tmp.xY,'X0') && ~isempty(tmp.xY.X0)
                    X0 = tmp.xY.X0;
                end
            end

            Ymat(:,k) = tmp.Y;
        end

        % --------------------------------------------------
        % Build DCM.Y (fMRI data)
        % --------------------------------------------------
        DCM         = struct();
        DCM.n       = nROI;   % number of regions
        DCM.v       = T;      % number of time points

        DCM.Y.y    = Ymat;
        DCM.Y.name = regionNames;
        DCM.Y.X0   = X0;
        DCM.Y.dt   = RT;

        % Minimal xY struct (enough for SPM)
        for k = 1:nROI
            DCM.xY(k).name     = roiNames{k};
            DCM.xY(k).modality = 'fMRI';
        end

        % --------------------------------------------------
        % Build DCM.U from SPM.Sess(r).U
        % --------------------------------------------------
        U_sess  = SPM.Sess(r).U;
        nInputs = numel(U_sess);

        % --------------------------------------------------
        % EARLIER CODE RETAINED FOR THE REVISION AUDIT TRAIL.
        % Do not use: it interpolates the padded SPM input to scan
        % resolution without first removing SPM's 32 leading bins.
        % --------------------------------------------------
        %{
        Umat  = zeros(T, nInputs);
        Uname = cell(1, nInputs);

        for ui = 1:nInputs
            u_vec = full(U_sess(ui).u(:,1));   % main effect regressor
            nU    = numel(u_vec);

            % Try to get input dt; fall back sensibly if missing
            if isfield(U_sess(ui),'dt') && ~isempty(U_sess(ui).dt)
                dt_u = U_sess(ui).dt;    % time step of u_vec
            else
                % Fallback: assume it is already at scan resolution
                dt_u = RT;
            end

            if nU == T
                % Case 1: already one value per scan
                reg = u_vec;
            else
                % Case 2: high-resolution input that must be resampled to scans
                t_u = (0:nU-1) * dt_u;
                t_y = (0:T-1) * RT;
                reg = interp1(t_u, u_vec, t_y, 'linear', 0);
            end

            Umat(:, ui) = reg;
            Uname{ui}   = U_sess(ui).name{1};
        end

        DCM.U.u    = Umat;        % [T x nInputs]
        DCM.U.name = Uname;       % 1 x nInputs cell
        DCM.U.dt   = RT;          % one value per scan
        %}

        % --------------------------------------------------
        % CORRECTED CODE.
        % SPM stores 32 leading microtime bins in each Sess.U input.
        % The standard DCM construction removes them with u(33:end,:)
        % and retains the microtime resolution rather than interpolating
        % the input to one value per scan.
        % --------------------------------------------------
        expectedInputRows  = T * SPM.xBF.T;
        expectedPaddedRows = expectedInputRows + 32;
        Umat  = zeros(expectedInputRows, nInputs);
        Uname = cell(1, nInputs);

        for ui = 1:nInputs
            if ~isfield(U_sess(ui),'dt') || isempty(U_sess(ui).dt)
                error('Missing input time step for %s, run %02d, input %d.', ...
                    sub, r, ui);
            end

            u_vec = full(U_sess(ui).u(:,1));   % main effect regressor
            if numel(u_vec) ~= expectedPaddedRows
                error(['Unexpected SPM input length for %s, run %02d, input %d: ' ...
                    'found %d, expected %d (= %d scans x %d bins + 32).'], ...
                    sub, r, ui, numel(u_vec), expectedPaddedRows, ...
                    T, SPM.xBF.T);
            end

            if abs(double(U_sess(ui).dt) - double(SPM.xBF.dt)) > 1e-12
                error(['Input time-step mismatch for %s, run %02d, input %d: ' ...
                    'Sess.U.dt = %.12g but SPM.xBF.dt = %.12g.'], ...
                    sub, r, ui, U_sess(ui).dt, SPM.xBF.dt);
            end

            Umat(:,ui) = u_vec(33:end);        % remove 32-bin SPM padding
            Uname{ui}  = U_sess(ui).name{1};
        end

        DCM.U.u    = Umat;                     % [T*SPM.xBF.T x nInputs]
        DCM.U.name = Uname;                    % 1 x nInputs cell
        DCM.U.dt   = SPM.xBF.dt;               % microtime step in seconds

        inputPeriod  = size(DCM.U.u,1) * DCM.U.dt;
        outputPeriod = DCM.v * DCM.Y.dt;
        if abs(inputPeriod - outputPeriod) > 1e-8
            error(['Input/output period mismatch for %s, run %02d: ' ...
                'input period %.12g s, output period %.12g s.'], ...
                sub, r, inputPeriod, outputPeriod);
        end

        % --------------------------------------------------
        % Specify Connectivity
        % --------------------------------------------------
        % A: full connectivity
        DCM.a = ones(nROI);

        % B: one modulation per input
        DCM.b = zeros(nROI, nROI, nInputs);
        %for bi = 1:(nInputs-1)
        %    DCM.b(:,:,bi) = eye(nROI);  % self-modulation
        %end

        DCM.b(2,1,1) = 1;   % vmPFC → vStr
        DCM.b(1,2,1) = 1;   % vStr → vmPFC

        % --- LossDom (input 2): negative affect / salience (H2, H3) ---
        % Amygdala ↔ aIns
        DCM.b(4,3,2) = 1;   % Amygd → aIns
        DCM.b(3,4,2) = 1;   % aIns → Amygd
        % Salience/affect → valuation
        DCM.b(1,3,2) = 1;   % Amygd → vmPFC
        DCM.b(1,4,2) = 1;   % aIns → vmPFC

        % C: driving inputs → vmPFC (region 1)
        DCM.c = zeros(nROI, nInputs);
        DCM.c(1, :) = 1;

        % --------------------------------------------------
        % DCM options
        % --------------------------------------------------
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 0;
        DCM.options.stochastic = 0;
        DCM.options.centre     = 1;

        % --------------------------------------------------
        % Save specification
        % --------------------------------------------------
        % Earlier output name retained for reference. Using it here would
        % overwrite the originally fitted DCM.
        % outFile = fullfile(firstLvl, sprintf('DCM_run%02d.mat', r));

        % Save the corrected fit under a new name.
        %outFile = fullfile(firstLvl, sprintf('DCM_corrected_run%02d.mat', r));
        outFile = fullfile(outputDir, sprintf('DCM_corrected_run%02d.mat', r));
        if exist(outFile,'file')
            error(['Corrected output already exists and will not be overwritten: %s. ' ...
                'Move or rename it explicitly before rerunning this subject.'], outFile);
        end
        save(outFile, 'DCM');
        fprintf('Saved DCM spec: %s\n', outFile);

        % --------------------------------------------------
        % Invert DCM
        % --------------------------------------------------
        fprintf('Inverting DCM for run %d...\n', r);
        DCM = spm_dcm_estimate(DCM);

        save(outFile, 'DCM');
        fprintf('✓ Finished inversion for run %02d\n', r);

    end

    fprintf('\nAll DCMs created & inverted successfully for subject %s.\n', sub);

end
