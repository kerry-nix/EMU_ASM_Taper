function ai = calc_ai(labels,thing,name,labels_plot,uni,last_dim,...
    which_thing,subplot_path,do_plots)

%{
This function takes a list of features for different electrode contacts,
finds the mesial temporal ones, and calculates a single asymmetry index
(for each frequency band) representing the L-R difference.
%}

which_elecs = {'A','B','C'};
which_lats = {'L','R'};

%% Misc
% replace '-' with '--'
if isempty(labels)
    ai = nan(1,last_dim);
    return
end
labels = cellfun(@(x) strrep(x,'-','--'),labels,'uniformoutput',false);

%% Get numeric info
number = cellfun(@(x) (regexp(x,'\d*','Match')),labels,'uniformoutput',false);
% replace empty with really high one
number(cellfun(@isempty,number)) = {{'9999'}};

% Get the first number per label (e.g., LA4--LA5 -> 4)
number = cellfun(@(x) str2num(x{1}),number);

%% Get which electrode
% get the letters
letters = cellfun(@(x) regexp(x,'[a-zA-Z]+','Match'),labels,'uniformoutput',false);
letters(cellfun(@isempty,letters)) = {{'zzzzz'}};
letters = cellfun(@(x) x{1},letters,'uniformoutput',false);

maxn = 12; % up to 12 contacts per electrode
nmt = length(which_elecs);

%% calculate the AI measurement
% initialize
intra = nan(nmt,maxn,2,last_dim); % n elecs, ncontacts, L and R, nfreq
intra_labels = cell(nmt,maxn,2);

% which electrodes
for i = 1:nmt
    
    % which laterality
    for j = 1:2
        curr_elec = [which_lats{j},which_elecs{i}];

        % which of the 12 contacts
        for k = 1:maxn

            % if this label exists
            if sum(ismember(labels,sprintf('%s%d',curr_elec,k)))>0
                intra_labels(i,k,j) = {sprintf('%s%d',curr_elec,k)};
            end
        end
        
        if uni == 1
           
            % Can do this on individual contact level
            for k = 1:maxn
       
                % Find the contacts matching this contact
                matching_contacts = strcmp(letters,curr_elec) & number == k;

                % should just have one match
                assert(sum(matching_contacts) <= 1)
                       
                % Calculate "intra" for these contacts, just the
                % thing for those contacts
                curr_intra = nanmean(thing(matching_contacts,:,:),1);
                
                % Fill
                intra(i,k,j,:) = curr_intra;
    
            end
          

        elseif uni == 0 % bivariate measures, will do intra-electrode connectivity
                
            %
            % Measure intra-electrode connectivity
            curr_intra = nanmean(thing(strcmp(letters,curr_elec),strcmp(letters,curr_elec),:),[1 2]);
            
            % Fill, repeating across all electrodes
            intra(i,:,j,:) = repmat(curr_intra,1,maxn,1,1);
            %}

            %{
            for k = 1:maxn
                % Find the contacts matching this contact
                matching_contacts = strcmp(letters,curr_elec) & number == k;

                % Measure average connectivity within the electode for this
                % contact
                curr_intra = nanmean(thing(matching_contacts,strcmp(letters,curr_elec),:),[1 2]);

                % Fill
                intra(i,k,j,:) = curr_intra;

            end
            %}
           
        end

    end

end

%% Take AI
ai = (intra(:,:,1,:)-intra(:,:,2,:))./((intra(:,:,1,:)+intra(:,:,2,:)));
ai = (squeeze(nanmean(ai,[1 2 3])))';
        
%% Average across ms  
if last_dim == 1
    ai = nanmean(ai);
end

%% Plot
if do_plots
    for d = 1:last_dim
        curr_intra = intra(:,:,:,d);
        curr_ai = ai(d);
        curr_thing = sprintf('%s_%d',which_thing{1},d);
        show_ai_electrodes(curr_intra,curr_ai,which_elecs,which_lats,name,subplot_path,curr_thing,labels_plot)
    end
end

%% Error checking
% I would expect that the labels with symmetric coverage should have values
% for intra (unless nan for other reasons)
if 0
    a = permute(intra_labels,[2,3,1]);
    b = permute(intra,[2,3,1]);
    table(a(:,:,1),b(:,:,1))
    table(a(:,:,2),b(:,:,2))
    table(a(:,:,3),b(:,:,3))
    ai
    pause
end

%% Checking nelecs
if 0
    a = permute(intra_labels,[2,3,1]);
    table(a(:,:,1))
    table(a(:,:,2))
    table(a(:,:,3))
    ai
    pause
end

end