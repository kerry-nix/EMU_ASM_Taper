function dist_matrix = pseudo_distance_mt(labels)

%% Assumptions about distances (obv not perfect)
% I got these by plotting electrode locs for a single patient and measuring
% distances
lr_gap = 30;
inter_elec_gap = 5;
inter_contact_gap = 5;

nchs = length(labels);
dist_matrix = nan(nchs,nchs);

for ich = 1:nchs
    for jch = ich+1:nchs
        label1 = labels{ich};
        label2 = labels{jch};

        if isempty(label1) || isempty(label2)
            continue
        end

        %% Get laterality gap
        if contains(label1,'L') && contains(label2,'L')
            lat_gap = 0;
        elseif contains(label1,'R') && contains(label2,'R')
            lat_gap = 0;
        elseif contains(label1,'L') && contains(label2,'R')
            lat_gap = lr_gap;
        elseif contains(label1,'R') && contains(label2,'L')
            lat_gap = lr_gap;
        else
            error('what');
        end

        %% Get inter-electrode gap
        % Get 2nd part of label
        if ~ismember(label1(2),{'A','B','C'}) || ~ismember(label2(2),{'A','B','C'})
            error('what');
        end
        if strcmp(label1(2),label2(2))
            ie_gap = 0;
        elseif (strcmp(label1(2),'A') && strcmp(label1(2),'C')) && (strcmp(label1(2),'C') && strcmp(label1(2),'A'))
            ie_gap = inter_elec_gap*2;
        else
            ie_gap = inter_elec_gap;
        end

        %% Get inter-contact gap
        % Get the numerical portion of the label
        number1 = regexp(label1,'\d*','Match');
        number2 = regexp(label2,'\d*','Match');
        ic_gap = abs(str2num(number1{1})-str2num(number2{1}))*inter_contact_gap;

        %% Total gap
        total_gap = lat_gap + ie_gap + ic_gap;
        dist_matrix(ich,jch) = total_gap;
        dist_matrix(jch,ich) = total_gap;

    end
end

end