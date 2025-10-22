function avg_pc = new_pearson_calc(values,fs,tw,do_tw)

nchs = size(values,2);


if ~do_tw
    avg_pc = corrcoef(values);
else
    %% Define time windows
    iw = round(tw*fs);
    window_start = 1:iw:size(values,1);

    % remove dangling window
    if window_start(end) + iw > size(values,1)
        window_start(end) = [];
    end
    nw = length(window_start);


    %% initialize output vector
    all_pc = nan(nchs,nchs,nw);

    %% Calculate pc for each window
    % Loop over time windows - note I have parallelized this step!
    for i = 1:nw

        % Define the time clip
        clip = values(window_start:window_start+iw,:);


        pc = corrcoef(clip);
        pc(logical(eye(size(pc)))) = 0;

        %% unwrap the pc matrix into a one dimensional vector for storage
        all_pc(:,:,i) = pc;

    end

    %% Average the network over all time windows
    avg_pc = nanmean(all_pc,3);

end

end