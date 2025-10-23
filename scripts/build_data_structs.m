% Script written by Nina Ghosn, PhD and Kerry Nix, University of Pennsylvania

addpath('/Volumes/USERS/knix/ASM_Taper_Sz/aed_dose_modeling/figures_code')
addpath('/Volumes/USERS/knix/ASM_Taper_Sz/DATA');
addpath('/Volumes/USERS/knix/ASM_Taper_Sz/aed_dose_modeling/')
addpath('/Volumes/USERS/knix/ASM_Taper_Sz/ASM_spikes_analysis/')
addpath('/Volumes/USERS/knix/ASM_Taper_Sz/helper code')

%close all; clear;

% load the generated ASM curves
load('EMU_cohort_asm_curves.mat')

data_spreadsheet = 'mar_and_timings_08142024.xlsx';
all_meds = readtable(data_spreadsheet);
cohort_info = readtable(data_spreadsheet);
mrn_to_ptid = readtable(data_spreadsheet,'Sheet','MRN_to_HUP');
events = readtable(data_spreadsheet,'Sheet','EVENTS');

% preserve "NaN" values of GTCs history in pat_history
opts = detectImportOptions(data_spreadsheet, 'Sheet', 'PAT_HISTORY','VariableNamesRange', '1:1');
opts = setvaropts(opts, 'GTCs', 'Type', 'double');
pat_history = readtable(data_spreadsheet, opts);

%%
mrn_to_ptid.HUP = str2double(mrn_to_ptid.HUP);
no_data_inds = cellfun(@length,mrn_to_ptid.MRN_DEIDENT) <1;
mrn_to_ptid(no_data_inds,:)=[];

% check which patients are common
cohort_info_asm = readtable('HUP_implant_dates.xlsx');
ptIDs_asm = cohort_info_asm.ptID;
ptIDs =unique(cohort_info.MRN_DEIDENT);

for i = 1:length(ptIDs)
    % how many unique ptIDs just before the '_'
    ptID = strsplit(ptIDs{i},'_');
    ptID=ptID(1);
    ptIDs_unique(i) = ptID;
end
ptIDs_unique =unique(ptIDs_unique);

is_common = false(height(mrn_to_ptid),1);
for i=1:length(is_common)
    is_common(i) = any(mrn_to_ptid.HUP(i) == ptIDs_asm);
end
mrn_to_ptid.is_common = is_common;


%% filter patients by exclusion criteria

%% 1. no medication data
meds_inds = cellfun(@iscell,all_dose_curves);
all_dose_curves = all_dose_curves(meds_inds);
all_Hr = all_Hr(meds_inds);
all_med_names = all_med_names(meds_inds);
ptIDs = ptIDs(meds_inds);


%% 2. no seizures and list seizure type
implant_date =  datetime('01-Jan-2000');

all_sz_types = table();
types =  unique(events.event_type);
types{1}='none';
all_sz_types.type = types;
sz_counts = zeros(length(types),1);

no_sz = false(length(ptIDs),1);
is_intracranial = false(length(ptIDs),1);

for i = 1:length(no_sz)
    [szs,is_implant] = get_sz_inds_nlp(ptIDs(i),events);
    if ~isempty(szs)
        none_inds = cellfun(@isempty,szs.sz_types);
        szs.sz_types(none_inds) = {'none'};
        szs(days(szs.sz_times- implant_date) > 30,:) = [];
        for s = 1:height(szs)
            ind = contains(types,szs.sz_types(s));
            sz_counts(ind) = sz_counts(ind)+1;
        end
    end
    no_sz(i)= isempty(szs);
    is_intracranial(i) = is_implant;

end

all_sz_types.sz_counts = sz_counts;


 %% Uncomment to exclude patients without any events
all_dose_curves_sz = all_dose_curves(~no_sz);
all_Hr = all_Hr(~no_sz);
all_med_names = all_med_names(~no_sz);
ptIDs = ptIDs(~no_sz);
is_intracranial = is_intracranial(~no_sz);


%% use scalp - exclude intracranial patients
use_scalp=1;
if use_scalp==1
    all_dose_curves_sz = all_dose_curves_sz(~is_intracranial);
    all_Hr = all_Hr(~is_intracranial);
    all_med_names = all_med_names(~is_intracranial);
    ptIDs = ptIDs(~is_intracranial);
end


%%
all_pts_taper_info = cell(length(ptIDs),1);
for ipt = 1:length(ptIDs)
    % patient medication administration
    ptID = ptIDs(ipt);
    pt_weight = weights(ipt);
    disp(ptID);

    [med_names, pt_meds, explant_dates, implant_dates] = parse_MAR_nlp(ptID, all_meds);
    
    % skip patients whose med list is missing or empty
    if isempty(pt_meds) || (istable(pt_meds) && height(pt_meds)==0)
        all_pts_taper_info{ipt} = [];
        continue;
    end

    [taper_info, any_taper] = get_taper_info(pt_meds, med_names);

    all_pts_taper_info{ipt} = taper_info;
end


%% Build data structures

% Get which seizures are followed by ativan and the drug levels
all_seizures =table(); % add ptID, seizure ID, preictal AED load, and 
% binary yes/no for ativan admin w/in 1hr, time from admission to event

aed_decrease = zeros(3,length(ptIDs));
all_drug_curves = containers.Map('KeyType', 'char', 'ValueType', 'any'); % key = ptID, value = drug_sum
all_times = [];
seizure_offsets = [];
start_ind =1;
all_implant_dates=NaT(1,length(ptIDs));
all_explant_dates=NaT(1,length(ptIDs));
initial_load=NaN(length(ptIDs),1);
num_meds = NaN(length(ptIDs),1);
time_to_first_seizure = zeros(length(ptIDs),1);

for ipt=1:length(ptIDs)

    ptID = ptIDs(ipt);

    %get drug curve to grab daily levels
    [~,meds,explant_date,implant_date] = parse_MAR_nlp(ptID,all_meds);
    med_names = all_med_names{ipt};

    all_implant_dates(ipt) =  implant_date(1);
    all_explant_dates(ipt) = explant_date(1);


    % get the total AED dose over time
    drugs =zeros(length(med_names),round(600*60)); %450 hours of EMU stay in minutes
    all_dstarts = [];
    for i =1:length(med_names)
        drug=all_dose_curves_sz{ipt}{i};
        
        if ~isempty(drug)
            drug=drug./(nanmean(drug(1:60))); %normalize each drug curve to the average of first hour
            dStart = round(all_Hr{ipt}{i}(1)*60)-1;
            all_dstarts = [all_dstarts dStart];
            drugs(i,dStart+1:dStart+length(drug))=drug;
        end
    end
    drug_sum=nansum(drugs,1);
    drug_sum=drug_sum./(length(med_names)); %normalize for number of drugs
    drug_sum(drug_sum==0) =NaN; %not include all zeros in average, and in histogram

    % store drug curves per patient
    all_drug_curves(ptID{1}) = drug_sum;

    % med curves are sampled at one point per minute
    t0 = min(all_dstarts); %earliest time in min a drug was administered
    t1 = t0 + (24*60);
    t2 = t1 + (24*60);
    %t3 = t2 + (24*60);

    day1_decrease = (nanmean(drug_sum(t0:t1))-nanmean(drug_sum(t1:t2)))./ (nanmean(drug_sum(t0:t1)));
    initial_load(ipt) = nanmean(drug_sum(t0:t1)); % the average drug load on day 1;
    num_meds(ipt) = length(med_names);
    %day2_decrease = (nanmean(drug_sum(t0:t1))-nanmean(drug_sum(t2:t3)))./ (nanmean(drug_sum(t0:t1)));

    aed_decrease(1,ipt) = day1_decrease;
    %aed_decrease(1,ipt) = day2_decrease;

    % get seizure times
    [szs,is_implant] = get_sz_inds_nlp(ptID,events);
    szs(round(szs.sz_ind) > length(drug_sum),:) =[];
    seizure_inds = round(szs.sz_ind);

    % convert seizure times to indices, so to minutes
    time_to_first_seizure(ipt)=seizure_inds(1); %the first seizure

    % get the pre-ictal load before each seizure
    % med curves are sampled at one point per minute
    pt_preictal_1hr = zeros(1,length(seizure_inds));
    for i=1:length(seizure_inds)
        if seizure_inds(i) < 60
            pt_preictal_1hr(i) = nanmean(drug_sum(1:seizure_inds(i)));
        else
            pt_preictal_1hr(i) = nanmean(drug_sum(seizure_inds(i)-60:seizure_inds(i)));
        end
    end

    end_ind = start_ind + length(seizure_inds);

    %Get ativan times
    ativan_inds = strcmp(meds.medication,'lorazepam');
    times = meds.admin_time(ativan_inds);

    seizure_inds = round(seizure_inds./60); % covert to hours to compare to med admin
    time_to_closest_sz = nan(length(seizure_inds),2);
    for n=1:height(seizure_inds)
        sz_diffs = seizure_inds(n)-times;
        before_ativan = sz_diffs< 0;
        if ~isempty(times)
            if ~isempty(before_ativan) && ~(sum(before_ativan)==0)
                [mval,~]=min(abs(sz_diffs(before_ativan)));
                ind = find(abs(sz_diffs)==mval);
                time_to_closest_sz(n,:)=[times(ind(1)) seizure_inds(n)];
            end
        end
    end

    warning('off')
    all_seizures.ptID(start_ind:end_ind-1)=ptID;
    all_seizures.seizureEEC(start_ind:end_ind-1) = seizure_inds;
    all_seizures.t_closest_ativan(start_ind:end_ind-1) =  time_to_closest_sz(:,1);
    all_seizures.preictal_aed_load(start_ind:end_ind-1) = pt_preictal_1hr';
    all_seizures.sz_type(start_ind:end_ind-1)=szs.sz_types;

    start_ind = end_ind;
end

%% Build seizure-related data structures

% find ativan administration
all_seizures.ativan_sz = double(all_seizures.t_closest_ativan - all_seizures.seizureEEC <= 4 &  all_seizures.t_closest_ativan - all_seizures.seizureEEC >=0);
unclear = strcmp(all_seizures.sz_type,'');
all_seizures.sz_type(unclear)={'Unclear'};
all_seizures(all_seizures.seizureEEC<0,:)=[];

has_conv = false(length(ptIDs),1);
has_any_conv = false(length(ptIDs),1);
was_GTC = false(length(ptIDs),1);
sz_type_tbl = cell(length(ptIDs),1); % sz event types
first_sz_preictal = zeros(length(ptIDs),1);

%  excluding PNES events
all_seizures_noPNES = all_seizures; % for building data structure excluding PNES
all_seizures_noPNES(strcmp(all_seizures_noPNES.sz_type, 'PNES'), :) = [];

has_conv_noPNES = false(length(ptIDs),1);
first_sz_preictal_noPNES = zeros(length(ptIDs),1);
sz_type_tbl_noPNES =cell(length(ptIDs),1); % sz event types
was_GTC_noPNES = false(length(ptIDs),1);
has_any_conv_noPNES = false(length(ptIDs),1);

for i=1:length(ptIDs)
    pt_inds = find(strcmp(all_seizures.ptID, ptIDs(i)));
    pt_inds_noPNES = find(strcmp(all_seizures_noPNES.ptID, ptIDs(i)));
    if sum(pt_inds)>0
        has_conv(i) = (all_seizures.ativan_sz(pt_inds(1))); % first sz only
        % has_any_conv(i) = any(all_seizures.ativan_sz(pt_inds)); % any sz with ativan
        first_sz_preictal(i) = (all_seizures.preictal_aed_load(pt_inds(1)));
        sz_type_tbl{i} = all_seizures.sz_type{pt_inds(1)};
        was_GTC(i) = strcmp(all_seizures.sz_type{pt_inds(1)}, 'GTC'); %True if event_type is GTC
        % check if ANY seizure for this patient was GTC
        has_any_conv(i) = any(strcmp(all_seizures.sz_type(pt_inds), 'GTC'));
    end

    % Create arrays excluding PNES events
    if sum(pt_inds_noPNES)>0
        has_conv_noPNES(i) = (all_seizures_noPNES.ativan_sz(pt_inds_noPNES(1))); % first sz only
        % has_any_conv(i) = any(all_seizures.ativan_sz(pt_inds)); % any sz with ativan
        first_sz_preictal_noPNES(i) = (all_seizures_noPNES.preictal_aed_load(pt_inds_noPNES(1)));
        sz_type_tbl_noPNES{i} = all_seizures_noPNES.sz_type{pt_inds_noPNES(1)};
        was_GTC_noPNES(i) = strcmp(all_seizures_noPNES.sz_type{pt_inds_noPNES(1)}, 'GTC'); %True if event_type is GTC
        % check if ANY seizure for this patient was GTC
        has_any_conv_noPNES(i) = any(strcmp(all_seizures_noPNES.sz_type(pt_inds_noPNES), 'GTC'));
    end
end


%% Mark clustered seizures and cluster leaders in both datasets
datasets = {'all_seizures','all_seizures_noPNES'};

for d = 1:numel(datasets)
    all_sz = eval(datasets{d});
    
    % initialize new columns
    all_sz.clustered_sz   = zeros(height(all_sz),1); 
    all_sz.cluster_leader = zeros(height(all_sz),1);
    
    for i = 1:length(ptIDs)
        % get this patient's seizures
        pt_inds = find(strcmp(all_sz.ptID, ptIDs(i)));
        if isempty(pt_inds), continue; end
        
        % sort by seizure time
        [sz_times, sort_idx] = sort(all_sz.seizureEEC(pt_inds));
        sorted_inds = pt_inds(sort_idx);
        
        % initialize leader and clustered sz flags
        leader_flags    = zeros(size(sorted_inds));
        clustered_flags = zeros(size(sorted_inds));
        
        j = 1;
        while j <= length(sz_times)
            % start a potential cluster
            cluster_inds = j;
            
            % collect all seizures within 4h of this one
            k = j + 1;
            while k <= length(sz_times) && (sz_times(k) - sz_times(j)) < hours(4)
                cluster_inds(end+1) = k; %#ok<AGROW>
                k = k + 1;
            end
            
            % mark leader + clustered seizures
            if numel(cluster_inds) > 1
                leader_flags(cluster_inds(1))        = 1; % cluster leader
                clustered_flags(cluster_inds(2:end)) = 1; % subsequent in cluster
            end
            
            j = j + 1;
        end
        
        % write back to table
        all_sz.cluster_leader(sorted_inds) = leader_flags;
        all_sz.clustered_sz(sorted_inds)   = clustered_flags;
    end
    assignin('base', datasets{d}, all_sz);
end


%% for all_seizures
% get baseline seizure frequencies
sz_freqs = nan(length(ptIDs),1);
gtcs = nan(length(ptIDs),1);
all_seizures.SzFreq=nan(height(all_seizures),1);
all_seizures.gtcs=nan(height(all_seizures),1);

for i = 1:length(ptIDs)
    ptID = strsplit(ptIDs{i},'_');
    ptID=ptID(1);
    ind = find(strcmp(pat_history.MRN_DEIDENT,ptID));
    if sum(ind)>0
        sz_freqs(i) = pat_history.SzFreq(ind);
        gtcs(i) = pat_history.GTCs(ind);

        % add to all_seizures
        ind2 = find(strcmp(all_seizures.ptID,ptIDs{i}));
        all_seizures.SzFreq(ind2) = pat_history.SzFreq(ind);
        all_seizures.gtcs(ind2) = pat_history.GTCs(ind);
    end
end

% add column to mark convulsive seizures ('GTC')
all_seizures.conv_sz = strcmp(all_seizures.sz_type,'GTC');


%% for all_seizures_noPNES
% get baseline seizure frequencies
sz_freqs = nan(length(ptIDs),1);
gtcs = nan(length(ptIDs),1);
all_seizures_noPNES.SzFreq=nan(height(all_seizures_noPNES),1);
all_seizures_noPNES.gtcs=nan(height(all_seizures_noPNES),1);

for i = 1:length(ptIDs)
    ptID = strsplit(ptIDs{i},'_');
    ptID=ptID(1);
    ind = find(strcmp(pat_history.MRN_DEIDENT,ptID));
    if sum(ind)>0
        sz_freqs(i) = pat_history.SzFreq(ind);
        gtcs(i) = pat_history.GTCs(ind);

        % add to all_seizures
        ind2 = find(strcmp(all_seizures_noPNES.ptID,ptIDs{i}));
        all_seizures_noPNES.SzFreq(ind2) = pat_history.SzFreq(ind);
        all_seizures_noPNES.gtcs(ind2) = pat_history.GTCs(ind);
    end
end

% add column to mark convulsive seizures ('GTC')
all_seizures_noPNES.conv_sz = strcmp(all_seizures_noPNES.sz_type,'GTC');

%% Create tables (tbl=including PNES and tbl_noPNES=excluding PNES)

%average half life [of medications tapered]
tbl = table();
tbl.asm_decrease = aed_decrease(1,:)'; %decrease from day 1 to 2
tbl.ptID = ptIDs;
tbl.has_conv = double(has_conv);
tbl.has_any_conv = double(has_any_conv);
tbl.gtcs = gtcs;
tbl.SzFreq = sz_freqs;
tbl.first_sz_preictal = first_sz_preictal;
tbl.was_GTC = was_GTC;
tbl.sz_type_tbl = sz_type_tbl;
%number of medications 
%find -inf and negative values and make them zero for no decrease
zero_inds = tbl.asm_decrease == -inf | tbl.asm_decrease <0;
tbl.asm_decrease(zero_inds)=0;


% Create tbl without PNES events
tbl_noPNES= table();
tbl_noPNES.asm_decrease = aed_decrease(1,:)'; %decrease from day 1 to 2
tbl_noPNES.ptID = ptIDs;
tbl_noPNES.has_conv = double(has_conv_noPNES);
tbl_noPNES.gtcs = gtcs;
tbl_noPNES.SzFreq = sz_freqs;
tbl_noPNES.first_sz_preictal = first_sz_preictal_noPNES;
tbl_noPNES.was_GTC = was_GTC_noPNES;
tbl_noPNES.sz_type_tbl = sz_type_tbl_noPNES;
tbl_noPNES.has_any_conv = has_any_conv_noPNES;
