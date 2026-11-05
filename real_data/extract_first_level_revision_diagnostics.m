%% Extract first-level information requested by Reviewer 3
% Run this script from the directory containing folders 01--11.
% It reads existing files only and writes results to a new timestamped folder.

clear;
clc;

%% -------- EDIT THESE SETTINGS -------------------------------------------
projectRoot = pwd;
spmDir = getenv('SPM_DIR');
if isempty(spmDir)
    error('Set the SPM_DIR environment variable to the SPM25 installation directory.');
end

% Use the same participant-information file and columns used in the
% task-version sensitivity analysis. Leave participantFile empty only if
% task version will be merged later in R.
participantFile = '';
subjectColumn = 'participant_id';
taskVersionColumn = 'task_version';

dcmRoot = fullfile(projectRoot, '01_run_level_DCM_files');
glmRoot = fullfile(projectRoot, '02_first_level_GLM');
pebRoot = fullfile(projectRoot, '04_within_subject_PEB');

%% -----------------------------------------------------------------------
addpath(spmDir);
spm('defaults', 'FMRI');

requiredFolders = {dcmRoot, glmRoot, pebRoot};
for k = 1:numel(requiredFolders)
    assert(isfolder(requiredFolders{k}), 'Folder not found: %s', requiredFolders{k});
end

stamp = datestr(now, 'yyyymmdd_HHMMSS');
outDir = fullfile(projectRoot, ['13_revision_diagnostics_' stamp]);
mkdir(outDir);

metadata = read_metadata(participantFile, subjectColumn, taskVersionColumn);
parameterSpecs = requested_parameter_specs();

conditionRows = struct([]);
designRows = struct([]);
fitRows = struct([]);
runVarianceRows = struct([]);
pebVarianceRows = struct([]);

subjectDirs = dir(fullfile(dcmRoot, 'sub-*'));
subjectDirs = subjectDirs([subjectDirs.isdir]);
assert(~isempty(subjectDirs), 'No sub-* folders found under %s', dcmRoot);

for s = 1:numel(subjectDirs)
    subject = string(subjectDirs(s).name);
    taskVersion = lookup_task_version(metadata, subject);

    spmFile = find_one_file(fullfile(glmRoot, char(subject)), 'SPM.mat');
    S = load(spmFile, 'SPM');
    assert(isfield(S, 'SPM'), 'SPM structure not found in %s', spmFile);
    SPM = S.SPM;

    dcmFiles = dir(fullfile(dcmRoot, char(subject), 'DCM_corrected_run*.mat'));
    assert(~isempty(dcmFiles), 'No corrected run-level DCM files found for %s', subject);
    [~, order] = sort({dcmFiles.name});
    dcmFiles = dcmFiles(order);

    for f = 1:numel(dcmFiles)
        dcmFile = fullfile(dcmFiles(f).folder, dcmFiles(f).name);
        runNumber = parse_run_number(dcmFiles(f).name);
        assert(runNumber <= numel(SPM.Sess), 'Run %d is absent from SPM.Sess for %s', runNumber, subject);

        D = load(dcmFile);
        DCM = find_struct_with_fields(D, {'Ep', 'Cp', 'U'});

        % Original first-level condition counts for this subject and run.
        for c = 1:numel(SPM.Sess(runNumber).U)
            U = SPM.Sess(runNumber).U(c);
            conditionName = string(U.name{1});
            row = struct( ...
                'subject', subject, ...
                'task_version', taskVersion, ...
                'run', runNumber, ...
                'condition', conditionName, ...
                'n_events', numel(U.ons));
            conditionRows = append_struct(conditionRows, row);
        end

        % DCM input rank and first-level GLM ranks.
        Umat = full(double(DCM.U.u));
        inputRank = numeric_rank(Umat);
        nonzeroInputColumns = sum(any(abs(Umat) > 1e-12, 1));

        runRows = SPM.Sess(runNumber).row;
        XrunAll = full(SPM.xX.X(runRows, :));
        taskCols = unique(SPM.Sess(runNumber).col(:)');
        taskCols = taskCols(taskCols >= 1 & taskCols <= size(SPM.xX.X, 2));
        XrunTask = full(SPM.xX.X(runRows, taskCols));

        row = struct( ...
            'subject', subject, ...
            'task_version', taskVersion, ...
            'run', runNumber, ...
            'n_dcm_input_columns', size(Umat, 2), ...
            'n_nonzero_dcm_input_columns', nonzeroInputColumns, ...
            'dcm_input_rank', inputRank, ...
            'n_glm_task_columns', size(XrunTask, 2), ...
            'glm_task_rank', numeric_rank(XrunTask), ...
            'n_glm_all_columns', size(XrunAll, 2), ...
            'glm_run_rank', numeric_rank(XrunAll));
        designRows = append_struct(designRows, row);

        % Run-level DCM fit diagnostics.
        [varianceExplained, residualRMS] = dcm_fit_quality(DCM);
        Cp = full(double(DCM.Cp));
        Cp = (Cp + Cp') / 2;
        cpEigenvalues = eig(Cp);
        cpFinite = all(isfinite(Cp(:)));
        cpPositive = cpFinite && all(cpEigenvalues > -1e-10);

        freeEnergy = NaN;
        if isfield(DCM, 'F') && isscalar(DCM.F)
            freeEnergy = double(DCM.F);
        end

        row = struct( ...
            'subject', subject, ...
            'task_version', taskVersion, ...
            'run', runNumber, ...
            'free_energy', freeEnergy, ...
            'variance_explained_percent', varianceExplained, ...
            'residual_rms', residualRMS, ...
            'posterior_covariance_finite', cpFinite, ...
            'posterior_covariance_positive_semidefinite', cpPositive, ...
            'posterior_covariance_min_eigenvalue', min(cpEigenvalues), ...
            'posterior_covariance_rcond', rcond(Cp));
        fitRows = append_struct(fitRows, row);

        % Run-level prior-to-posterior variance changes for the 22 parameters.
        if isfield(DCM, 'M') && isfield(DCM.M, 'pE') && isfield(DCM.M, 'pC')
            pE = DCM.M.pE;
            pC = full(double(DCM.M.pC));
            for p = 1:numel(parameterSpecs)
                spec = parameterSpecs(p);
                idx = structure_parameter_index(pE, spec);
                if idx <= size(pC, 1) && idx <= size(Cp, 1)
                    priorVariance = pC(idx, idx);
                    posteriorVariance = Cp(idx, idx);
                    row = variance_row(subject, taskVersion, runNumber, spec.name, ...
                        priorVariance, posteriorVariance, 'run_level_DCM');
                    runVarianceRows = append_struct(runVarianceRows, row);
                end
            end
        end
    end

    % Within-subject PEB prior-to-posterior variance changes.
    pebFile = find_subject_peb(pebRoot, subject);
    P = load(pebFile);
    PEB = find_struct_with_fields(P, {'Ep', 'Cp'});
    if isfield(PEB, 'M') && isfield(PEB.M, 'pC')
        pebCp = full(double(PEB.Cp));
        pebPrior = full(double(PEB.M.pC));
        pebNames = get_peb_parameter_names(PEB, parameterSpecs);
        for p = 1:numel(parameterSpecs)
            idx = find(strcmp(normalize_parameter_names(pebNames), ...
                normalize_parameter_names({parameterSpecs(p).name})), 1);
            % Some SPM PEB files decorate names with field or cell indices.
            % If the PEB contains exactly the requested 22 parameters and a
            % decorated name cannot be matched, retain their established order.
            if isempty(idx) && size(PEB.Cp, 1) == numel(parameterSpecs)
                idx = p;
            end
            if ~isempty(idx) && idx <= size(pebPrior, 1) && idx <= size(pebCp, 1)
                row = variance_row(subject, taskVersion, NaN, parameterSpecs(p).name, ...
                    pebPrior(idx, idx), pebCp(idx, idx), 'within_subject_PEB');
                pebVarianceRows = append_struct(pebVarianceRows, row);
            end
        end
    else
        warning('Prior covariance was not found in within-subject PEB file %s', pebFile);
    end
end

writetable(struct2table(conditionRows), fullfile(outDir, 'subject_run_condition_counts.csv'));
writetable(struct2table(designRows), fullfile(outDir, 'subject_run_design_ranks.csv'));
writetable(struct2table(fitRows), fullfile(outDir, 'run_level_dcm_fit_diagnostics.csv'));
writetable(struct2table(runVarianceRows), fullfile(outDir, 'run_level_variance_changes.csv'));
writetable(struct2table(pebVarianceRows), fullfile(outDir, 'within_subject_peb_variance_changes.csv'));

fprintf('\nWrote revision diagnostics to:\n%s\n', outDir);
fprintf('Send the five CSV files in this folder for manuscript-ready numerical text.\n');


%% --------------------------- HELPERS ------------------------------------
function metadata = read_metadata(fileName, subjectColumn, taskVersionColumn)
metadata = table();
if strlength(string(fileName)) == 0
    return;
end
assert(isfile(fileName), 'Participant-information file not found: %s', fileName);
metadata = readtable(fileName, 'FileType', 'text', 'VariableNamingRule', 'preserve');
assert(ismember(subjectColumn, metadata.Properties.VariableNames), ...
    'Subject column %s not found in %s', subjectColumn, fileName);
assert(ismember(taskVersionColumn, metadata.Properties.VariableNames), ...
    'Task-version column %s not found in %s', taskVersionColumn, fileName);
metadata.subject_key = normalize_subject_ids(metadata.(subjectColumn));
metadata.task_version_key = string(metadata.(taskVersionColumn));
end

function taskVersion = lookup_task_version(metadata, subject)
taskVersion = "";
if isempty(metadata)
    return;
end
key = normalize_subject_ids(subject);
idx = find(metadata.subject_key == key, 1);
if ~isempty(idx)
    taskVersion = metadata.task_version_key(idx);
end
end

function keys = normalize_subject_ids(x)
x = string(x);
keys = strings(size(x));
for k = 1:numel(x)
    token = regexp(char(x(k)), 'sub-\d+', 'match', 'once');
    if isempty(token)
        keys(k) = strtrim(x(k));
    else
        keys(k) = string(token);
    end
end
end

function fileName = find_one_file(folder, pattern)
files = dir(fullfile(folder, '**', pattern));
assert(~isempty(files), 'No %s found under %s', pattern, folder);
assert(numel(files) == 1, 'Expected one %s under %s, found %d', pattern, folder, numel(files));
fileName = fullfile(files(1).folder, files(1).name);
end

function fileName = find_subject_peb(pebRoot, subject)
% The exported subject names use the exact form PEB_sub-XXX_runwise.
% Search for that exact file before considering any broader fallback.
expectedName = sprintf('PEB_%s_runwise.mat', char(subject));
files = dir(fullfile(pebRoot, '**', expectedName));

if isempty(files)
    allFiles = dir(fullfile(pebRoot, '**', '*.mat'));
    names = lower(string({allFiles.name}));
    keep = contains(names, lower(subject)) & ...
           contains(names, 'peb') & ...
           contains(names, 'runwise') & ...
          ~contains(names, 'backup') & ...
          ~contains(names, 'old') & ...
          ~contains(names, 'original');
    files = allFiles(keep);
end

assert(~isempty(files), 'No runwise within-subject PEB file found for %s', subject);

if numel(files) > 1
    candidatePaths = string(fullfile({files.folder}, {files.name}));
    error('Multiple runwise PEB files remain for %s:\n%s', ...
        subject, strjoin(candidatePaths, newline));
end

fileName = fullfile(files(1).folder, files(1).name);
end

function output = find_struct_with_fields(loaded, requiredFields)
names = fieldnames(loaded);
for k = 1:numel(names)
    candidate = loaded.(names{k});
    if isstruct(candidate) && all(isfield(candidate, requiredFields))
        output = candidate;
        return;
    end
end
error('No structure with required fields was found in the loaded MAT file.');
end

function runNumber = parse_run_number(fileName)
token = regexp(fileName, 'run(\d+)', 'tokens', 'once');
assert(~isempty(token), 'Could not parse run number from %s', fileName);
runNumber = str2double(token{1});
end

function r = numeric_rank(X)
if isempty(X)
    r = 0;
    return;
end
s = svd(double(X), 'econ');
if isempty(s)
    r = 0;
else
    tolerance = max(size(X)) * eps(max(s));
    r = sum(s > tolerance);
end
end

function [pve, residualRMS] = dcm_fit_quality(DCM)
pve = NaN;
residualRMS = NaN;
if ~isfield(DCM, 'R')
    return;
end
R = double(DCM.R);
residualRMS = sqrt(mean(R(isfinite(R)).^2));
if isfield(DCM, 'Y') && isstruct(DCM.Y) && isfield(DCM.Y, 'y')
    Y = double(DCM.Y.y);
    if isequal(size(Y), size(R))
        Ycentered = Y - mean(Y, 1, 'omitnan');
        yv = Ycentered(isfinite(Ycentered));
        rv = R(isfinite(R));
        if ~isempty(yv) && ~isempty(rv) && sum(yv.^2) > 0
            pve = 100 * (1 - sum(rv.^2) / sum(yv.^2));
        end
    end
end
end

function specs = requested_parameter_specs()
specs = struct('name', {}, 'field', {}, 'indices', {});
for j = 1:4
    for i = 1:4
        specs(end + 1) = struct('name', sprintf('A(%d,%d)', i, j), ...
            'field', 'A', 'indices', [i, j]); %#ok<AGROW>
    end
end
b = [2 1 1; 1 2 1; 1 3 2; 4 3 2; 1 4 2; 3 4 2];
for k = 1:size(b, 1)
    specs(end + 1) = struct('name', sprintf('B(%d,%d,%d)', b(k,1), b(k,2), b(k,3)), ...
        'field', 'B', 'indices', b(k,:)); %#ok<AGROW>
end
end

function idx = structure_parameter_index(template, spec)
zeroTemplate = spm_unvec(zeros(size(spm_vec(template))), template);
if strcmp(spec.field, 'A')
    zeroTemplate.A(spec.indices(1), spec.indices(2)) = 1;
else
    zeroTemplate.B(spec.indices(1), spec.indices(2), spec.indices(3)) = 1;
end
idx = find(abs(spm_vec(zeroTemplate)) > 0.5, 1);
assert(~isempty(idx), 'Could not locate parameter %s in the DCM parameter vector', spec.name);
end

function row = variance_row(subject, taskVersion, runNumber, parameter, priorVariance, posteriorVariance, source)
ratio = NaN;
reduction = NaN;
informationGain = NaN;
if isfinite(priorVariance) && priorVariance > 0 && isfinite(posteriorVariance) && posteriorVariance >= 0
    ratio = posteriorVariance / priorVariance;
    reduction = 1 - ratio;
    if posteriorVariance > 0
        informationGain = 0.5 * log(priorVariance / posteriorVariance);
    end
end
row = struct( ...
    'subject', subject, ...
    'task_version', taskVersion, ...
    'run', runNumber, ...
    'parameter', string(parameter), ...
    'source', string(source), ...
    'prior_variance', priorVariance, ...
    'posterior_variance', posteriorVariance, ...
    'posterior_to_prior_ratio', ratio, ...
    'variance_reduction_fraction', reduction, ...
    'information_gain_nats', informationGain);
end

function names = get_peb_parameter_names(PEB, specs)
if isfield(PEB, 'Pnames') && numel(PEB.Pnames) == size(PEB.Cp, 1)
    names = cellstr(string(PEB.Pnames));
elseif size(PEB.Cp, 1) == numel(specs)
    names = {specs.name};
else
    error('Could not map the within-subject PEB covariance to the 22 requested parameters.');
end
end

function y = normalize_parameter_names(x)
y = regexprep(string(x), '\s+', '');
y = regexprep(y, '[{}]', '');
y = regexprep(y, '^A1\(', 'A(');
y = regexprep(y, '^B1\(', 'B(');
end

function array = append_struct(array, row)
if isempty(array)
    array = row;
else
    array(end + 1) = row;
end
end
