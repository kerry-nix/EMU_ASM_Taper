metaData = table2struct(readtable('metadata.csv','Delimiter',',','Format','auto'));
for i = 1:length(metaData)
    data = fetch_ieeg_info(metaData(i).filename, login.usr, login.pwd);
    metaData(i).duration = data.duration;
    metaData(i).channels = data.chLabels;
    metaData(i).fs = data.fs;
    metaData(i).chanstr = strjoin(metaData(i).channels');
end

function data = fetch_ieeg_info(fname, login_name, pwfile)

try
    session = IEEGSession(fname, login_name, pwfile);
    channelLabels = session.data.channelLabels(:,1);
    % get fs
    data.fs = session.data.sampleRate;

    % Get ch labels
    data.chLabels = channelLabels;
    if size(channelLabels,1)>0
        % get duration (convert to seconds)
        data.duration = session.data.rawChannels(1).get_tsdetails.getDuration/1e6;
    else
        data.duration = nan;
    end
catch ME
    data.fs = [];
    data.chLabels = {};
    data.duration = [];
end
end