cfg = get_project_config();

clear; clc;

edf_file = fullfile(cfg.dataRoot, 'Master', 'Subjects', 's1', 'ali-shahab1.edf');
ipmFile  = fullfile(cfg.dataRoot, 'ipm2_standard_64.ced');

EEG = pop_biosig(edf_file, 'channels', 1:64, 'importevent','off');
EEG = eeg_checkset(EEG);

edfLabels = strings(64,1);
for i = 1:64
    if isfield(EEG,'chanlocs') && numel(EEG.chanlocs) >= i && ...
            isfield(EEG.chanlocs,'labels') && ~isempty(EEG.chanlocs(i).labels)
        edfLabels(i) = string(EEG.chanlocs(i).labels);
    else
        edfLabels(i) = "Ch" + i;
    end
end

% Load IPM montage labels
locs = readlocs(ipmFile);
ipmLabels = string({locs(1:64).labels})';

MapT = table((1:64)', edfLabels, ipmLabels, ...
    'VariableNames', {'ChannelIndex','EDF_Label_After_popbiosig','IPM2_Label'});

disp(MapT)

outPath = fullfile(cfg.outputRoot, 'edf_vs_ipm2_channel_map.csv');
writetable(MapT, outPath);

fprintf('Saved: %s\n', outPath);
disp('Important channels:')
disp(MapT([10 13 32 59 63],:))

%% 

locs = readlocs(fullfile(cfg.dataRoot, 'ipm2_standard_64.ced'));

labels = string({locs.labels})';

X = [locs.X]';
Y = [locs.Y]';
Z = [locs.Z]';

T = table((1:64)', labels(1:64), X(1:64), Y(1:64), Z(1:64), ...
    isfinite(X(1:64)) & isfinite(Y(1:64)) & isfinite(Z(1:64)), ...
    'VariableNames', {'ChannelIndex','Label','X','Y','Z','HasXYZ'});

disp(T)
fprintf('Channels with valid XYZ: %d / 64\n', sum(T.HasXYZ));
