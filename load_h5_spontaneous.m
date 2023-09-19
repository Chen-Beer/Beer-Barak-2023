%% load data

directory_files = dir(dir_);
files_names = {directory_files.name};
files_ = {files_names{3:end}};
file = [dir_,files_{(contains(files_,'.h5'))}];
data = McsHDF5.McsData(file);
timeStream = data.Recording{1}.TimeStampStream{1};
Elabels = timeStream.Info.Label;
activity_raw = cell(length(Elabels),1);
for e = 1:length(Elabels)
    if length(Elabels{e}) < 3
        Elabels{e} = strcat(Elabels{e}(1),'0',Elabels{e}(end));
    end
    activity_raw{e} = timeStream.TimeStamps{e};
end
N = length(Elabels);
electrodes_names = Elabels;
T = double(max([activity_raw{:}]));