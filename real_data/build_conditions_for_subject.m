function [onsets, durs] = build_conditions_for_subject(sub, bidsRoot, nRuns)

onsets.GainDom  = cell(nRuns,1);
onsets.LossDom  = cell(nRuns,1);
onsets.Middle   = cell(nRuns,1);
durs.GainDom    = cell(nRuns,1);
durs.LossDom    = cell(nRuns,1);
durs.Middle     = cell(nRuns,1);

for r = 1:nRuns
    evFile = fullfile(bidsRoot, sub, 'func', ...
        sprintf('%s_task-MGT_run-%02d_events.txt', sub, r));
    T = readtable(evFile, 'FileType','text');

    gain  = T.gain;
    loss  = T.loss;
    onset = T.onset;
    dur   = T.duration;              % should be 4 s in NARPS, but use file

    ratio = gain ./ max(loss, 0.01); % avoid division by zero

    ixGain = find(ratio >= 2);
    ixLoss = find(ratio <= 0.5);
    ixMid  = setdiff(1:numel(onset), [ixGain; ixLoss]);

    onsets.GainDom{r} = onset(ixGain);
    durs.GainDom{r}   = dur(ixGain);

    onsets.LossDom{r} = onset(ixLoss);
    durs.LossDom{r}   = dur(ixLoss);

    onsets.Middle{r}  = onset(ixMid);
    durs.Middle{r}    = dur(ixMid);
end
end
