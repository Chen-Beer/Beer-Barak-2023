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
stimulation1 = data.Recording{1}.EventStream{1}.Events{1}(1,:);
stimulation2 = data.Recording{1}.EventStream{1}.Events{5}(1,:);
stimulation3 = data.Recording{1}.EventStream{1}.Events{9}(1,:);

stimulation = sort([stimulation1 stimulation2 stimulation3]);

diff_stimulation = diff(stimulation);
[~,ii] = find(diff_stimulation>60e6);
spont_off = stimulation(ii+1);
spont_on = spont_off - 60e6;

electrodes_names = Elabels;

