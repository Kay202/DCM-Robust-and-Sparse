function output_directory = audit_first_level_dcm_inputs(root_directory)
% Compare event onsets, SPM session inputs, and saved DCM inputs.
%
% Place this file beside folders 01_run_level_DCM_files,
% 02_first_level_GLM, and 03_event_files, then run:
%
%   audit_first_level_dcm_inputs
%
% Alternatively, provide the parent directory explicitly:
%
%   audit_first_level_dcm_inputs('C:\path\to\revision_package')
%
% Existing files are never changed. Results are written to a new,
% timestamped directory each time the function is run.

if nargin < 1 || isempty(root_directory)
    root_directory = fileparts(mfilename('fullpath'));
end
root_directory = char(root_directory);

dcm_root = fullfile(root_directory, '01_run_level_DCM_files');
spm_root = fullfile(root_directory, '02_first_level_GLM');
event_root = fullfile(root_directory, '03_event_files');
assert(isfolder(dcm_root), 'Missing folder: %s', dcm_root);
assert(isfolder(spm_root), 'Missing folder: %s', spm_root);
assert(isfolder(event_root), 'Missing folder: %s', event_root);

output_directory = unique_output_directory(root_directory);
fprintf('Writing new audit files to:\n%s\n\n', output_directory);

listing = dir(fullfile(dcm_root, 'sub-*'));
subjects = sort({listing([listing.isdir]).name});
assert(~isempty(subjects), 'No sub-* folders found in %s.', dcm_root);

condition_names = {'GainDom', 'LossDom', 'Middle'};
number_of_rows = numel(subjects) * 4 * numel(condition_names);
result = cell(number_of_rows, 29);
row = 0;

for subject_number = 1:numel(subjects)
    subject = subjects{subject_number};
    fprintf('Auditing %s (%d/%d)\n', subject, subject_number, numel(subjects));

    spm_file = fullfile(spm_root, subject, 'SPM.mat');
    assert(isfile(spm_file), 'Missing file: %s', spm_file);
    loaded = load(spm_file, 'SPM');
    SPM = loaded.SPM;

    for run_number = 1:4
        event_file = fullfile(event_root, subject, sprintf( ...
            '%s_task-MGT_run-%02d_events.txt', subject, run_number));
        dcm_file = fullfile(dcm_root, subject, ...
            sprintf('DCM_run%02d.mat', run_number));
        assert(isfile(event_file), 'Missing file: %s', event_file);
        assert(isfile(dcm_file), 'Missing file: %s', dcm_file);

        events = readtable(event_file, 'FileType', 'text', ...
            'VariableNamingRule', 'preserve');
        onset = table_column(events, 'onset');
        duration = table_column(events, 'duration');
        gain = table_column(events, 'gain');
        loss = table_column(events, 'loss');

        % This exactly reproduces build_conditions_for_subject.m.
        ratio = gain ./ max(loss, 0.01);
        indices = cell(1, 3);
        indices{1} = find(ratio >= 2);
        indices{2} = find(ratio <= 0.5);
        indices{3} = setdiff((1:numel(onset))', ...
            [indices{1}; indices{2}]);

        loaded = load(dcm_file, 'DCM');
        DCM = loaded.DCM;
        scan_count = double(SPM.nscan(run_number));
        scan_rt = double(SPM.xY.RT);
        microtime_bins = double(SPM.xBF.T);
        spm_dt = double(SPM.xBF.dt);
        spm_units = char(SPM.xBF.UNITS);

        for condition_number = 1:3
            row = row + 1;
            condition = condition_names{condition_number};
            event_onsets = onset(indices{condition_number});
            event_durations = duration(indices{condition_number});

            spm_index = condition_index(SPM.Sess(run_number).U, condition);
            dcm_index = name_index(DCM.U.name, condition);
            assert(~isnan(spm_index), ...
                '%s run %d: %s is absent from SPM.Sess.U.', ...
                subject, run_number, condition);
            assert(~isnan(dcm_index), ...
                '%s run %d: %s is absent from DCM.U.name.', ...
                subject, run_number, condition);

            U = SPM.Sess(run_number).U(spm_index);
            spm_vector = full(double(U.u(:, 1)));
            dcm_vector = full(double(DCM.U.u(:, dcm_index)));
            spm_vector = spm_vector(:);
            dcm_vector = dcm_vector(:);

            % Directly compare event onsets with the onsets saved in SPM.
            spm_onsets = double(U.ons(:));
            if strcmpi(spm_units, 'scans')
                spm_onsets = spm_onsets * scan_rt;
            end
            if numel(spm_onsets) == numel(event_onsets)
                onset_difference = max_or_zero(abs( ...
                    sort(spm_onsets) - sort(event_onsets)));
            else
                onset_difference = NaN;
            end

            % SPM stores 32 leading microtime bins. First verify that the
            % event file reproduces this stored vector, including the padding.
            expected_padded = stimulus_vector(event_onsets, event_durations, ...
                numel(spm_vector), spm_dt, 32);
            [padded_lag, padded_correlation] = best_lag( ...
                expected_padded, spm_vector, 64);

            % SPM's spm_dcm_U removes these bins using U.u(33:end,:).
            assert(numel(spm_vector) > 32, ...
                '%s run %d has an unexpectedly short SPM input.', ...
                subject, run_number);
            unpadded_spm = spm_vector(33:end);
            expected_unpadded = stimulus_vector(event_onsets, ...
                event_durations, numel(unpadded_spm), spm_dt, 0);
            [unpadded_lag, unpadded_correlation] = best_lag( ...
                expected_unpadded, unpadded_spm, 64);

            % Compare DCM.U with both the earlier manual construction and the
            % offset-corrected SPM input. The saved DCM may be at scan or
            % microtime resolution, so both cases are handled.
            legacy_difference = NaN;
            corrected_difference = NaN;
            corrected_dcm_lag = NaN;
            corrected_dcm_correlation = NaN;
            legacy_dcm_lag = NaN;
            legacy_dcm_correlation = NaN;
            dcm_is_centered = isfield(DCM, 'options') && ...
                isfield(DCM.options, 'centre') && logical(DCM.options.centre);

            if numel(dcm_vector) == scan_count
                representation = 'scan-resolution';
                legacy_projection = project_to_scans( ...
                    spm_vector, spm_dt, scan_count, scan_rt);
                corrected_projection = project_to_scans( ...
                    unpadded_spm, spm_dt, scan_count, scan_rt);
                if dcm_is_centered
                    legacy_projection = legacy_projection - mean(legacy_projection);
                    corrected_projection = corrected_projection - mean(corrected_projection);
                end
                legacy_difference = max(abs(dcm_vector - legacy_projection));
                corrected_difference = max(abs( ...
                    dcm_vector - corrected_projection));
                [legacy_dcm_lag, legacy_dcm_correlation] = best_continuous_lag( ...
                    legacy_projection, dcm_vector, 8);
                [corrected_dcm_lag, corrected_dcm_correlation] = best_continuous_lag( ...
                    corrected_projection, dcm_vector, 8);
                corrected_dcm_lag_seconds = corrected_dcm_lag * scan_rt;
            elseif numel(dcm_vector) == numel(unpadded_spm)
                representation = 'SPM-microtime';
                corrected_projection = unpadded_spm;
                if dcm_is_centered
                    corrected_projection = corrected_projection - mean(corrected_projection);
                end
                corrected_difference = max(abs(dcm_vector - corrected_projection));
                [corrected_dcm_lag, corrected_dcm_correlation] = ...
                    best_continuous_lag(corrected_projection, dcm_vector, 64);
                corrected_dcm_lag_seconds = corrected_dcm_lag * spm_dt;
            else
                representation = 'unexpected-length';
                corrected_dcm_lag_seconds = NaN;
            end

            dcm_dt = NaN;
            if isfield(DCM.U, 'dt') && ~isempty(DCM.U.dt)
                dcm_dt = double(DCM.U.dt);
            end

            likely_32_bin_offset = ...
                strcmp(representation, 'scan-resolution') && ...
                isfinite(corrected_dcm_lag_seconds) && ...
                abs(abs(corrected_dcm_lag_seconds) - 32 * spm_dt) ...
                    <= max(scan_rt / 2, spm_dt) && ...
                corrected_dcm_correlation >= 0.80 && ...
                (legacy_difference < corrected_difference || ...
                 (legacy_dcm_lag == 0 && ...
                  legacy_dcm_correlation >= corrected_dcm_correlation));

            other_misalignment = ...
                numel(spm_onsets) ~= numel(event_onsets) || ...
                (isfinite(onset_difference) && onset_difference > spm_dt) || ...
                abs(padded_lag) > 1 || abs(unpadded_lag) > 1 || ...
                (isfinite(corrected_dcm_lag_seconds) && ...
                 abs(corrected_dcm_lag_seconds) > spm_dt && ...
                 ~likely_32_bin_offset);

            result(row, :) = {subject, run_number, condition, ...
                numel(event_onsets), numel(spm_onsets), onset_difference, ...
                first_name(U.name), first_name(DCM.U.name(dcm_index)), ...
                spm_units, scan_count, microtime_bins, spm_dt, dcm_dt, ...
                numel(spm_vector), numel(dcm_vector), representation, ...
                padded_lag, padded_lag * spm_dt, padded_correlation, ...
                unpadded_lag, unpadded_lag * spm_dt, unpadded_correlation, ...
                corrected_dcm_lag_seconds, corrected_dcm_correlation, ...
                legacy_difference, corrected_difference, ...
                likely_32_bin_offset, other_misalignment, ...
                likely_32_bin_offset || other_misalignment};
        end
    end
end

variable_names = { ...
    'subject', 'run', 'condition', 'event_count', 'spm_onset_count', ...
    'maximum_event_vs_spm_onset_difference_seconds', ...
    'spm_condition_name', 'dcm_condition_name', 'spm_units', ...
    'scan_count', 'microtime_bins_per_scan', 'spm_dt_seconds', ...
    'dcm_dt_seconds', 'spm_input_length', 'dcm_input_length', ...
    'dcm_representation', 'event_vs_padded_spm_lag_bins', ...
    'event_vs_padded_spm_lag_seconds', ...
    'event_vs_padded_spm_correlation', ...
    'event_vs_unpadded_spm_lag_bins', ...
    'event_vs_unpadded_spm_lag_seconds', ...
    'event_vs_unpadded_spm_correlation', ...
    'corrected_input_vs_dcm_lag_seconds', ...
    'corrected_input_vs_dcm_correlation', ...
    'dcm_vs_legacy_projection_max_abs_difference', ...
    'dcm_vs_corrected_projection_max_abs_difference', ...
    'likely_32_microtime_bin_offset', ...
    'other_input_misalignment', ...
    'input_misalignment_flag'};

audit_table = cell2table(result, 'VariableNames', variable_names);
writetable(audit_table, fullfile(output_directory, 'dcm_input_audit.csv'));

summary_table = table( ...
    height(audit_table), ...
    sum(audit_table.input_misalignment_flag), ...
    sum(audit_table.likely_32_microtime_bin_offset), ...
    sum(audit_table.other_input_misalignment), ...
    sum(strcmp(audit_table.dcm_representation, 'scan-resolution')), ...
    sum(strcmp(audit_table.dcm_representation, 'SPM-microtime')), ...
    'VariableNames', {'rows_audited', 'rows_flagged', ...
    'rows_with_likely_32_bin_offset', 'rows_with_other_misalignment', ...
    'scan_resolution_DCM_rows', 'microtime_DCM_rows'});
writetable(summary_table, ...
    fullfile(output_directory, 'dcm_input_audit_summary.csv'));
save(fullfile(output_directory, 'dcm_input_audit.mat'), ...
    'audit_table', 'summary_table', 'root_directory');

fprintf('\nAudit complete.\n');
disp(summary_table);
fprintf('Detailed results: %s\n', ...
    fullfile(output_directory, 'dcm_input_audit.csv'));
end


function output_directory = unique_output_directory(root_directory)
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
base = fullfile(root_directory, ['dcm_input_audit_' stamp]);
output_directory = base;
number = 0;
while isfolder(output_directory) || isfile(output_directory)
    number = number + 1;
    output_directory = sprintf('%s_%02d', base, number);
end
mkdir(output_directory);
end


function x = table_column(T, requested_name)
names = T.Properties.VariableNames;
index = find(strcmpi(names, requested_name), 1);
assert(~isempty(index), 'Event file has no %s column.', requested_name);
x = T.(names{index});
if ~isnumeric(x)
    x = str2double(string(x));
end
x = double(x(:));
assert(all(isfinite(x)), 'Column %s contains missing values.', requested_name);
end


function index = condition_index(U, requested_name)
index = NaN;
for i = 1:numel(U)
    if strcmp(normalize(first_name(U(i).name)), normalize(requested_name))
        index = i;
        return
    end
end
end


function index = name_index(names, requested_name)
index = NaN;
for i = 1:numel(names)
    if strcmp(normalize(first_name(names(i))), normalize(requested_name))
        index = i;
        return
    end
end
end


function value = first_name(value)
if iscell(value)
    value = value{1};
end
if isstring(value)
    value = value(1);
end
value = char(value);
end


function value = normalize(value)
value = lower(regexprep(first_name(value), '[^a-zA-Z0-9]', ''));
end


function stimulus = stimulus_vector(onsets, durations, n, dt, padding)
temporary = zeros(n + 128, 1);
if isempty(onsets)
    stimulus = temporary(1:n);
    return
end
amplitude = ones(numel(onsets), 1);
if all(durations == 0)
    amplitude = amplitude / dt;
end
on_index = round(onsets / dt) + padding + 1;
off_index = round(durations / dt) + on_index + 1;
for i = 1:numel(on_index)
    on = max(on_index(i), 1);
    off = max(off_index(i), 1);
    if on <= numel(temporary)
        temporary(on) = temporary(on) + amplitude(i);
    end
    if off <= numel(temporary)
        temporary(off) = temporary(off) - amplitude(i);
    end
end
temporary = cumsum(temporary);
stimulus = temporary(1:n);
end


function projected = project_to_scans(x, dt, scan_count, scan_rt)
input_time = (0:(numel(x) - 1))' * dt;
scan_time = (0:(scan_count - 1))' * scan_rt;
projected = interp1(input_time, x, scan_time, 'linear', 0);
end


function [best_shift, best_correlation] = best_lag(expected, observed, max_shift)
expected = double(abs(expected(:)) > threshold(expected));
observed = double(abs(observed(:)) > threshold(observed));
best_shift = NaN;
best_correlation = -Inf;
max_shift = min([max_shift, numel(expected) - 2, numel(observed) - 2]);

for shift = -max_shift:max_shift
    if shift >= 0
        n = min(numel(expected), numel(observed) - shift);
        x = expected(1:n);
        y = observed((1 + shift):(shift + n));
    else
        advance = -shift;
        n = min(numel(expected) - advance, numel(observed));
        x = expected((1 + advance):(advance + n));
        y = observed(1:n);
    end
    if n < 3 || std(x) == 0 || std(y) == 0
        correlation = -Inf;
    else
        x = x - mean(x);
        y = y - mean(y);
        correlation = sum(x .* y) / sqrt(sum(x .^ 2) * sum(y .^ 2));
    end
    if correlation > best_correlation
        best_correlation = correlation;
        best_shift = shift;
    end
end
if ~isfinite(best_correlation)
    best_correlation = NaN;
end
end


function [best_shift, best_correlation] = best_continuous_lag(expected, observed, max_shift)
expected = double(expected(:));
observed = double(observed(:));
best_shift = NaN;
best_correlation = -Inf;
max_shift = min([max_shift, numel(expected) - 2, numel(observed) - 2]);

for shift = -max_shift:max_shift
    if shift >= 0
        n = min(numel(expected), numel(observed) - shift);
        x = expected(1:n);
        y = observed((1 + shift):(shift + n));
    else
        advance = -shift;
        n = min(numel(expected) - advance, numel(observed));
        x = expected((1 + advance):(advance + n));
        y = observed(1:n);
    end
    if n < 3 || std(x) == 0 || std(y) == 0
        correlation = -Inf;
    else
        x = x - mean(x);
        y = y - mean(y);
        correlation = sum(x .* y) / sqrt(sum(x .^ 2) * sum(y .^ 2));
    end
    if correlation > best_correlation
        best_correlation = correlation;
        best_shift = shift;
    end
end
if ~isfinite(best_correlation)
    best_correlation = NaN;
end
end


function value = threshold(x)
value = max(abs(double(x(:)))) * 1e-10;
end


function value = max_or_zero(x)
if isempty(x)
    value = 0;
else
    value = max(x);
end
end
