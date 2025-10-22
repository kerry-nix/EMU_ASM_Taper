function avg_pc = pc_vector_calc(values,fs,tw)

all_at_once = 0;
nchs = size(values,2);


if all_at_once
    avg_pc = corrcoef(values);
    avg_pc = wrap_or_unwrap_adjacency_fc_toolbox(avg_pc);
else
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


        pc = corrcoef(clip);
        pc(logical(eye(size(pc)))) = 0;

        %% unwrap the pc matrix into a one dimensional vector for storage
        all_pc(:,i) = wrap_or_unwrap_adjacency_fc_toolbox(pc);

    end

    %% Average the network over all time windows
    avg_pc = nanmean(all_pc,2);

end

end