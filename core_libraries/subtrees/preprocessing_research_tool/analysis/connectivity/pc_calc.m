function all_pc = pc_calc(values,fs,tw)

%{
Calculate pearson connectivity between each pair of channels
%}

nchs = size(values,2);

%% Define time windows
iw = tw*fs;
window_start = 1:iw:size(values,1);

% remove dangling window
if window_start(end) + iw > size(values,1)
    window_start(end) = [];
end
nw = length(window_start);


%% initialize output vector
all_pc = nan(nchs*(nchs-1)/2,nw);

%% Calculate pc for each window
% Loop over time windows - note I have parallelized this step!
for i = 1:nw
    
    % Define the time clip
    clip = values(window_start:window_start+iw,:);
    
    % Default set it to be nans
    pc = nan(nchs,nchs);
    
    % Set elements along the diagonal to zero
    pc(logical(eye(nchs))) = 0;
    
    
    for ich = 1:nchs
        
        % Skip it if all nans
        if sum(~isnan(clip(:,ich))) == 0
            continue
        end
        
        for jch = 1:ich-1 % check that this is the right num to loop through
            
            % Skip if all nans
            if sum(~isnan(clip(:,jch))) == 0
                continue
            end
            
            % pearson correlation between chs ich and jch
            r = corr(clip(:,ich),clip(:,jch));
            
            pc(ich,jch) = r;
            pc(jch,ich) = r;
        end
    end
    
    if 0
        figure
        imagesc(pc(:,:))
        pause
        close(gcf)
    end
    
    %% unwrap the pc matrix into a one dimensional vector for storage
    all_pc(:,i) = wrap_or_unwrap_adjacency_fc_toolbox(pc);
    
    
end

%% Average the network over all time windows
avg_pc = nanmean(all_pc,2);

end