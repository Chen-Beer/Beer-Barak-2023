%% load data

directory_files = dir(dir_);
files_names = {directory_files.name};
files_ = {files_names{3:end}};
file = [dir_,files_{(contains(files_,'.h5'))}];
data = McsHDF5.McsData(file);
timeStream = data.Recording{1}.TimeStampStream{1};
Elabels = timeStream.Info.Label;
SP = cell(length(Elabels),1);
for e = 1:length(Elabels)
    if length(Elabels{e}) < 3
        Elabels{e} = strcat(Elabels{e}(1),'0',Elabels{e}(end));
    end
    SP{e} = timeStream.TimeStamps{e};
end
N = length(Elabels);

if control_label
    directory_files_ = dir('../data/38428_17.11/4 exp/');
    files_names_ = {directory_files_.name};
    files__ = {files_names_{3:end}};
    data_ = McsHDF5.McsData(['../data/38428_17.11/4 exp/',files__{(contains(files__,'.h5'))}]);
    stimulation1 = data_.Recording{1}.EventStream{1}.Events{1}(1,:);
    stimulation2 = data_.Recording{1}.EventStream{1}.Events{5}(1,:);
    stimulation3 = data_.Recording{1}.EventStream{1}.Events{9}(1,:);
else
    stimulation1 = data.Recording{1}.EventStream{1}.Events{1}(1,:);
    stimulation2 = data.Recording{1}.EventStream{1}.Events{5}(1,:);
    stimulation3 = data.Recording{1}.EventStream{1}.Events{9}(1,:);
end

nn = min(min(length(stimulation1),length(stimulation2)),length(stimulation3));
stimulation = zeros(nn,1);
idx = 1;
for i = 1:nn
    stimulation(idx) = stimulation1(i);
    stimulation(idx+1) = stimulation2(i);
    stimulation(idx+2) = stimulation3(i);
    idx = idx + 3;
end

electrodes_names = Elabels;

spont_on = stimulation3(round_cycles:round_cycles:end) + 1e6;
spont_off = spont_on + 60e6;