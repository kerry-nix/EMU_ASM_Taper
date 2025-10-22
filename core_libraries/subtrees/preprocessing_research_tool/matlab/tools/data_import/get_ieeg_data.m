function data = get_ieeg_data(fname, login_name, pwfile, times, varargin)

% This function download data from ieeg portal with filename, time range,
% and user information provided
% Inputs:
%       fname: name of to be downloaded file
%       login_name: user name of ieeg portal account
%       pwfile: password file of ieeg portal account
%       times: time range of to be pulled data within file, in seconds
%       extras: whether to assign other data attributes to output var
%       selecElecs: cell array of electrode list of channel lables
%       or array of int channel indices to be included 
%       ignoreElecs: cell array of electrode list of channel lables
%       or array of int channel indices to be ignored
% Output:
%       data: a struct of field 'fs', 'values', 'file_name', 'chLabels',
%       'duration', and 'ann'

%% Added 1/22/23, Haoer

selecElecs = {};
ignoreElecs = {};
extras = 1;
assert(mod(nargin,2)==0,'CNTtools:invalidInput','Optional inputs should be paired.');
para = varargin;
while length(para) >= 2
    name = para{1};
    val = para{2};
    para = para(3:end);
    switch name
        case 'selecElecs'
            selecElecs = val;
        case 'ignoreElecs'
            ignoreElecs = val;
        case 'extras'
            extras = val;
    end
end

%% Added 11/12, Haoer

% Assert input type
assert((isstring(fname) || ischar(fname)),'CNTtools:invalidInputType','Invalid input type.')
assert(isa(times,'numeric'),'CNTtools:invalidInputType','Invalid input type.')

% remove potential blanks/quotes
fname = strip(fname);
fname = strip(fname,'"');
fname = strip(fname,"'");

% If with meta data file, metaData can be fetched through the tool/fetch_metadata 
% function, don't know how to effectively fetch yet
%if exist('meta_data.csv','file')
    
    % Import metaData
    %metaData = readtable('meta_data.csv','Delimiter',',','Format','auto');% read meta data
    
    % test if filename valid, obtain other info from metadata
    %assert(ismember(fname,metaData.filename),'CNTtools:invalidFileName','Invalid filename.')
    %cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
    %assert(~isempty(metaData.chans{cellfun(cellfind(fname),metaData.filename)}),'CNTtools:emptyFile','No channels.')
    %dura = metaData.duration(cellfun(cellfind(fname),metaData.filename));

%else
    % without meta data file, in case ob different data base
attempt = 1;

while attempt < 100
    try
        session = IEEGSession(fname, login_name, pwfile); % check if can fetch this file
        channelLabels = session.data.channelLabels; % get channel info
        break
    catch ME
        if contains(ME.message,'Authentication')
            throw(MException('CNTtools:invalidLoginInfo','Invalid login info.'));
        elseif contains(ME.message,'No snapshot with name')
            throw(MException('CNTtools:invalidFileName','Invalid filename.'));
        elseif contains(ME.message,'positive integers or logical values')
            throw(MException('CNTtools:emptyFile','No channels.'));
        elseif contains(ME.message,'503') || contains(ME.message,'504') || ...
                contains(ME.message,'502') || contains(ME.message,'500')
            attempt = attempt + 1;
            fprintf('Failed to retrieve ieeg.org data, trying again (attempt %d)\n',attempt); 
        else
            throw(ME);
        end
    end
end
%end
% check file is not empty
% assert(~isempty(channelLabels),'CNTtools:emptyFile','No channels.'); 
dura = session.data.rawChannels(1).get_tsdetails.getDuration/1e6;  % get duration info
% Delete session
session.delete;        
% test if time range valid
assert(times(2) > times(1),'CNTtools:invalidTimeRange','Stop before start.')
assert(times(1) >= 0 && times(2) <= dura, 'CNTtools:invalidTimeRange', 'Time outrange.')
%% Added 02/05/2023, Haoer
[channelLabels,~,~] = clean_labels(channelLabels);
allChannelLabels = channelLabels;
numElecs = length(allChannelLabels);
channelIDs = [1:numElecs];
if ~isempty(selecElecs)
    if ~iscell(selecElecs)
        try
            if istable(selecElecs)
                selecElecs = table2cell(selecElecs); 
            % convert to cell if string or char format
            elseif isstring(selecElecs) || ischar(selecElecs)
                selecElecs = cellstr(selecElecs); 
            elseif isa(selecElecs,'numeric')
                selecElecs = num2cell(selecElecs);
            end
        catch ME
            throw(MException('CNTtools:invalidInputType','Electrode labels should be a cell array.'))
        end
    end
    if isa(selecElecs{1},'numeric')
        channelIDs = selecElecs(cellfun(@(x) x >= 1 & x <= numElecs, selecElecs));
        if length(channelIDs) < length(selecElecs)
            warning("CNTtools:invalidChannelID, invalid channels ignored.");
        end
        channelLabels = allChannelLabels(cellfun(@(x) x, channelIDs));
        channelIDs = cell2mat(channelIDs);
    elseif isstring(selecElecs{1}) || ischar(selecElecs{1})
        [selecElecs,~,~] = clean_labels(selecElecs);
        [~, channelIDs] = ismember(selecElecs, allChannelLabels);
        channelIDs(channelIDs == 0) = [];
        channelLabels = allChannelLabels(channelIDs);
        %channelIDs = num2cell(channelIDs);
        if length(channelIDs) < length(selecElecs)
            warning("CNTtools:invalidChannelID, invalid channels ignored.");
        end
    end
end
if ~isempty(ignoreElecs)
    if ~iscell(ignoreElecs)
        try
            if istable(ignoreElecs)
                ignoreElecs = table2cell(ignoreElecs); 
            % convert to cell if string or char format
            elseif isstring(selecElecs) || ischar(ignoreElecs)
                ignoreElecs = cellstr(ignoreElecs); 
            elseif isa(ignoreElecs,'numeric')
                ignoreElecs = num2cell(ignoreElecs);
            end
        catch ME
            throw(MException('CNTtools:invalidInputType','Electrode labels should be a cell array.'))
        end
    end
    if isa(ignoreElecs{1},'numeric')
        channelIDs(cellfun(@(x) x >= 1 & x <= numElecs, ignoreElecs)) = [];
        %channelIDs = num2cell(channelIDs);
        if length(channelIDs) > numElecs - length(ignoreElecs)
            warning("CNTtools:invalidChannelID, invalid channels ignored.");
        end
        channelLabels = allChannelLabels(channelIDs);
    elseif isstring(ignoreElecs{1}) || ischar(ignoreElecs{1})
        [ignoreElecs,~,~] = clean_labels(ignoreElecs);
        channelLabels = allChannelLabels(cellfun(@(x) ~ismember(x,ignoreElecs), allChannelLabels));
        [~, channelIDs] = ismember(channelLabels, allChannelLabels);
        %channelIDs = num2cell(channelIDs);
        if length(channelIDs) > numElecs - length(ignoreElecs)
            warning("CNTtools:invalidChannelID, invalid channels ignored.");
        end
    end
end
% -----end of added code-----

%% original code

attempt = 1;

% Wrap data pulling attempts in a while loop
while attempt < 100
    try
        session = IEEGSession(fname, login_name, pwfile);
        % channelLabels = session.data.channelLabels;
        nchs = size(channelLabels,1);

        % get fs
        data.fs = session.data.sampleRate;

        % Convert times to indices
        run_idx = round(times(1)*data.fs):round(times(2)*data.fs);

        if ~isempty(run_idx)
            try 
                values = session.data.getvalues(run_idx,channelIDs);
            catch ME
                % Break the number of channels in half to avoid wacky server errors
                values1 = session.data.getvalues(run_idx,channelIDs(1:floor(nchs/4)));
                values2 = session.data.getvalues(run_idx,channelIDs(floor(nchs/4)+1:floor(2*nchs/4)));
                values3 = session.data.getvalues(run_idx,channelIDs(floor(2*nchs/4)+1:floor(3*nchs/4))); 
                values4 = session.data.getvalues(run_idx,channelIDs(floor(3*nchs/4)+1:nchs)); 
    
                values = [values1,values2,values3,values4];
            end
        else
            values = [];
        end

        data.values = values;

        if extras == 1

            % get file name
            data.file_name = session.data.snapName;

            % Get ch labels
            data.chLabels = channelLabels;

            % get duration (convert to seconds)
            data.duration = session.data.rawChannels(1).get_tsdetails.getDuration/1e6;

            % Get annotations
            n_layers = length(session.data.annLayer);

            for ai = 1:n_layers
                a = session.data.annLayer(ai).getEvents(0);
                n_ann = length(a);
                for i = 1:n_ann
                    event(i).start = a(i).start/(1e6);
                    event(i).stop = a(i).stop/(1e6); % convert from microseconds
                    event(i).type = a(i).type;
                    event(i).description = a(i).description;
                end
                ann.event = event;
                ann.name = session.data.annLayer(ai).name;
                data.ann(ai) = ann;
            end 

        end
        
        % break out of while loop
        break
        
    % If server error, try again (this is because there are frequent random
    % server errors).
    catch ME
        if contains(ME.message,'503') || contains(ME.message,'504') || ...
                contains(ME.message,'502') || contains(ME.message,'500')
            attempt = attempt + 1;
            fprintf('Failed to retrieve ieeg.org data, trying again (attempt %d)\n',attempt); 
        else
            throw(ME)
            error('Non-server error'); 
        end
        
    end
end

%% Delete session
session.delete;
clearvars -except data

end