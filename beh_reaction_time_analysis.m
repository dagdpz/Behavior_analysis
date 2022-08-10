function [f_input, out_comp, out_stru_ext, unique_pos] = beh_reaction_time_analysis(files_for_input,dates, batching,subject_files,steady,inactivation_sites)
global  GLO
GLO.dates                           = dates;
%% ??
dag_drive_IP='Y:';
GLO.temp_dir                    = [dag_drive_IP 'Data\' files_for_input{1}];
h                               = findstr('Data\', GLO.temp_dir);
GLO.monkey                      = GLO.temp_dir(h+5:end);

if GLO.delete_last
    for idx_batch                   = 1:numel(files_for_input_batch(:,1))
        DeleteLastTrial(files_for_input_batch{idx_batch,1});
    end
end
f_input{1} = subject_files;

%% saccade definition 4, first saccade, 3 last saccade
%% MAIN FUNCTION STRUCTURE, RUN MPA THEN RUN BY RUN ANALYZYS AND THEN BATCH ANALYSYS
correlation_conditions                                      = {'demanded_hand','choice','type','effector','target_side','success'};
parameters_to_correlate                                     = {'lat','dur'};
Sel_all = {'display',0,'nsacc_max',20,'correlation_conditions',correlation_conditions,'parameters_to_correlate',parameters_to_correlate,'runs_as_batches',batching.runs_as_batches_val};

%for trial by trial inspection 
%    Sel_all = {'display',1,'summary',0,'nsacc_max',5,'correlation_conditions',correlation_conditions,'parameters_to_correlate',parameters_to_correlate,'runs_as_batches',batching.runs_as_batches_val};
% 'reach_1st_pos'
% 'reach_1st_pos_in'
% 'reach_pos_at_state_change'

steady_FN=fieldnames(steady);
for f=1:numel(steady_FN)
    Sel_all={Sel_all{:},steady_FN{f},steady.(steady_FN{f})};
end

Sel_input=repmat({Sel_all},numel(subject_files),1);
MA_input=[subject_files, Sel_input]';
[out_comp,~,~]                                          = monkeypsych_analyze_working(MA_input{:});

if GLO.min_hits
    nn=1;
    for n=1:numel(out_comp)
        if  numel([out_comp{n}.binary.success]) > 50
            [out_comp_temp(nn,:)] = out_comp(n);
            nn=nn+1;
        end
    end
    out_comp = out_comp_temp;
end

[out_comp]                                              = reallocate_positions_from_mpa(out_comp);

for n=1:numel(out_comp)
    %add field containing index of aborted trials (logical)
    out_comp{n}.abort_raw_index     =  [[out_comp{n}.states.state_abo]==4 | [out_comp{n}.states.state_abo]==5 | [out_comp{n}.states.state_abo]==9 | [out_comp{n}.states.state_abo]== 10 ...
        | [out_comp{n}.states.state_abo]==6 | [out_comp{n}.states.state_abo]==7 | [out_comp{n}.states.state_abo]==8];
    %add field containing index of successful trials (logical)
    out_comp{n}.success_raw_index   =  [out_comp{n}.states.state_abo]==-1;
    
    % out_comp{n}.success_raw_index   =  [out_comp{n}.task.abort_code];
    % create fields containing saccades and reach informations (just
    % reorganise)
    tmp = {out_comp{n}.raw.x_eye};
    [out_comp{n}.saccades.raw_x] = tmp{:};
    tmp = {out_comp{n}.raw.y_eye};
    [out_comp{n}.saccades.raw_y] = tmp{:};
    tmp = {out_comp{n}.raw.x_hnd};
    [out_comp{n}.reaches.raw_x] = tmp{:};
    tmp = {out_comp{n}.raw.y_hnd};
    [out_comp{n}.reaches.raw_y] = tmp{:};
    tmp = {out_comp{n}.raw.states};
    [out_comp{n}.saccades.raw_states] = tmp{:};
    tmp = {out_comp{n}.raw.states};
    [out_comp{n}.reaches.raw_states] = tmp{:};
    tmp = {out_comp{n}.raw.time_axis};
    [out_comp{n}.saccades.raw_time_axis] = tmp{:};
    tmp = {out_comp{n}.raw.time_axis};
    [out_comp{n}.reaches.raw_time_axis] = tmp{:};
    tmp = {out_comp{n}.task.abort_code};
    [out_comp{n}.reaches.abort_code] = tmp{:};
    
    
    clear tmp
end

%[out_comp,~,~]                                          = monkeypsych_analyze_hand_eye_choice_seperately7(MA_input{:});

%% Unique positions
out_comp_mat                        = vertcat(out_comp{:});
out_comp_saccades                   = vertcat(out_comp_mat.saccades);
out_comp_reaches                    = vertcat(out_comp_mat.reaches);
saccadepositions                    = [out_comp_saccades.tar_pos]-[out_comp_saccades.fix_pos];
reachpositions                      = [out_comp_reaches.tar_pos]-[out_comp_reaches.fix_pos];
% saccadepositions                    = [out_comp_saccades.tar_pos];
% reachpositions                      = [out_comp_reaches.tar_pos];
% reorganise and calculate mean, raw, stfd, and num of hits for the bacth
unique_pos.saccades                 = unique_positions(saccadepositions,1.5);
unique_pos.reaches = unique_positions(reachpositions,1.5);
%get target/fixation radius for reaches and saccades
unique_pos.reaches_tar_rad = nanmean([out_comp_reaches.tar_rad]);
unique_pos.saccades_tar_rad = nanmean([out_comp_saccades.tar_rad]);
unique_pos.reaches_tar_siz = nanmean([out_comp_reaches.tar_siz]);
unique_pos.saccades_tar_siz = nanmean([out_comp_saccades.tar_siz]);
unique_pos.reaches_fix_rad = nanmean([out_comp_reaches.fix_rad]);
unique_pos.saccades_fix_rad = nanmean([out_comp_saccades.fix_rad]);
unique_pos.reaches_fix_siz = nanmean([out_comp_reaches.fix_siz]);
unique_pos.saccades_fix_siz = nanmean([out_comp_saccades.fix_siz]);


[out_str,  out_ini_fix,out_dur_fix, out_ini_abort, out_hnd_abort]              = rt_s_internal_cal(out_comp,unique_pos,inactivation_sites);
[out_stru_ext]                      = external_cal(out_str);

subparameters                       = {'mean','raw','std','num_hits'};
out_stru_ext.reaches.ini_fix.LH=get_external_means_std([out_ini_fix.LH.ini_fix],subparameters);
out_stru_ext.reaches.ini_fix.RH=get_external_means_std([out_ini_fix.RH.ini_fix],subparameters);

out_stru_ext.reaches.dur_fix.LH=get_external_means_std([out_dur_fix.LH.dur_fix],subparameters);
out_stru_ext.reaches.dur_fix.RH=get_external_means_std([out_dur_fix.RH.dur_fix],subparameters);

out_stru_ext.reaches.ini_abort.LH=get_external_means_std([out_ini_abort.LH.abort_code],subparameters);
out_stru_ext.reaches.ini_abort.RH=get_external_means_std([out_ini_abort.RH.abort_code],subparameters);

out_stru_ext.reaches.hnd_switch.LH=get_external_means_std([out_hnd_abort.LH.hnd_switch],subparameters);
out_stru_ext.reaches.hnd_switch.RH=get_external_means_std([out_hnd_abort.RH.hnd_switch],subparameters);

out_stru_ext.reaches.hnd_stay.LH=get_external_means_std([out_hnd_abort.LH.hnd_stay],subparameters);
out_stru_ext.reaches.hnd_stay.RH=get_external_means_std([out_hnd_abort.RH.hnd_stay],subparameters);
end

function  [out_str, out_ini_fix, out_dur_fix, out_ini_abort, out_hnd_abort]= rt_s_internal_cal(out_comp,unique_pos,inactivation_sites)
global GLO

% effector_names                  = {'0','2','3','4',GLO.type_of_free_gaze};
% type_names                      = {'2','3','4'};

effector_names                  = {'0','1','3','4',GLO.type_of_free_gaze};
type_names                      = {'2','3','4'};

side_names                      = {'L','R'};
decision_names                  = {'CH','IN'};
hand_names                      = {'LH','RH'};
sac_rea                         = {'saccades','reaches'};


% Indexes
for idx_batch=1:numel(out_comp)
    idx.reach_LH{idx_batch}                         =   [out_comp{idx_batch,1}.task.demanded_hand]==1;
    idx.reach_RH{idx_batch}                         =   [out_comp{idx_batch,1}.task.demanded_hand]==2;
    
    idx.hnd_switch{idx_batch}                       =   [false (idx.reach_RH{idx_batch}(2:end) & idx.reach_LH{idx_batch}(1:end-1)) |...
        (idx.reach_LH{idx_batch}(2:end) & idx.reach_RH{idx_batch}(1:end-1))];
    idx.hnd_stay{idx_batch}                         =   [false (idx.reach_RH{idx_batch}(2:end) & idx.reach_RH{idx_batch}(1:end-1)) |...
        (idx.reach_LH{idx_batch}(2:end) & idx.reach_LH{idx_batch}(1:end-1))];
    
    eye_horizontal_distance                         =   real([out_comp{idx_batch,1}.saccades.tar_pos] - [out_comp{idx_batch,1}.saccades.fix_pos]);
    hnd_horizontal_distance                         =   real([out_comp{idx_batch,1}.reaches.tar_pos] - [out_comp{idx_batch,1}.reaches.fix_pos]);
    
    if strcmp(inactivation_sites{idx_batch},'R') | strcmp(inactivation_sites{idx_batch},'Nan')
    idx.L_tar{idx_batch}                            =   eye_horizontal_distance <0.01 | hnd_horizontal_distance <0.01;
    idx.R_tar{idx_batch}                            =   eye_horizontal_distance >0.01 | hnd_horizontal_distance >0.01;
    
    idx.L_tar_far{idx_batch}                        =   eye_horizontal_distance <-15 | hnd_horizontal_distance <-15;
    idx.R_tar_far{idx_batch}                        =   eye_horizontal_distance >15 |  hnd_horizontal_distance >15;
    
    idx.L_tar_close{idx_batch}                      =   (eye_horizontal_distance <0.01    & eye_horizontal_distance>=-15) | (hnd_horizontal_distance<0.01    & hnd_horizontal_distance>=-15);
    idx.R_tar_close{idx_batch}                      =   (eye_horizontal_distance >0.01    & eye_horizontal_distance<=15)  | (hnd_horizontal_distance>0.01    & hnd_horizontal_distance<=15);

    elseif strcmp(inactivation_sites{idx_batch},'L') % if left = reverse sides
    idx.R_tar{idx_batch}                            =   eye_horizontal_distance <0.01 | hnd_horizontal_distance <0.01;
    idx.L_tar{idx_batch}                            =   eye_horizontal_distance >0.01 | hnd_horizontal_distance >0.01;
    
    idx.R_tar_far{idx_batch}                        =   eye_horizontal_distance <-15 | hnd_horizontal_distance <-15;
    idx.L_tar_far{idx_batch}                        =   eye_horizontal_distance >15 |  hnd_horizontal_distance >15;
    
    idx.R_tar_close{idx_batch}                      =   (eye_horizontal_distance <0.01    & eye_horizontal_distance>=-15) | (hnd_horizontal_distance<0.01    & hnd_horizontal_distance>=-15);
    idx.L_tar_close{idx_batch}                      =   (eye_horizontal_distance >0.01    & eye_horizontal_distance<=15)  | (hnd_horizontal_distance>0.01    & hnd_horizontal_distance<=15);
    end
    
    idx.CH{idx_batch}                               =   [out_comp{idx_batch,1}.binary.choice]==1;
    idx.IN{idx_batch}                               =   [out_comp{idx_batch,1}.binary.choice]==0;
    idx.success{idx_batch}                          =   [out_comp{idx_batch,1}.binary.success]==1;
    idx.error{idx_batch}                            =   [out_comp{idx_batch,1}.binary.success]==0;
    idx.error_after_success{idx_batch}              =   [false (idx.error{idx_batch}(2:end) & idx.success{idx_batch}(1:end-1))];
    %     idx.incorrect_hand{idx_batch}                   =   idx.error_after_success{idx_batch} & ismember({out_comp{idx_batch,1}.task.abort_code},'ABORT_USE_INCORRECT_HAND');
    idx.incorrect_hand{idx_batch}                   =   ismember({out_comp{idx_batch,1}.task.abort_code},'ABORT_USE_INCORRECT_HAND');
    idx.rew_mod{idx_batch}                          =   [out_comp{idx_batch,1}.task.reward_modulation]==1;
    
    idx.pre_control{idx_batch}                      =   [out_comp{idx_batch,1}.binary.microstim]==0;
    idx.stim{idx_batch}                             =   [out_comp{idx_batch,1}.binary.microstim]==1;
    
    idx.control{idx_batch}                          =   idx.pre_control{idx_batch}             & ~idx.rew_mod{idx_batch};
    
    batch_correlation_conditions                    =   [out_comp{idx_batch}.correlation.conditions];
    
    corr_idx.reach_LH{idx_batch}                    =   [batch_correlation_conditions.demanded_hand]==1;
    corr_idx.reach_RH{idx_batch}                    =   [batch_correlation_conditions.demanded_hand]==2;
    
    corr_idx.CH{idx_batch}                          =   [batch_correlation_conditions.choice]==1;
    corr_idx.IN{idx_batch}                          =   [batch_correlation_conditions.choice]==0;
    
    corr_idx.L_tar{idx_batch}                       =   [batch_correlation_conditions.target_side]==-1;
    corr_idx.R_tar{idx_batch}                       =   [batch_correlation_conditions.target_side]==1;
    
    idx.hnd_abort_tar_acq{idx_batch}                =  ismember({out_comp{idx_batch,1}.task.abort_code},'ABORT_HND_TAR_ACQ_STATE');
    
    for t=1:numel(type_names)
        for e                                       = 1:numel(effector_names)
            batch_effector                          = [out_comp{idx_batch}.task.effector];
            batch_type                              = [out_comp{idx_batch}.task.type];
            type_effector_indicator{t,e}(idx_batch) = any(batch_effector==str2double(effector_names{e}) & batch_type==str2double(type_names{t}));
        end
    end
    
    % ini_fix and dur_fix
    for h = 1:numel(hand_names)
        hand = hand_names{h};
        temp_index_h = idx.(['reach_' hand]){idx_batch};
        out_ini_fix.(hand)(idx_batch)= get_raw_mean_std(out_comp{idx_batch},temp_index_h,temp_index_h,'ini_fix','IN');
        out_ini_abort.(hand)(idx_batch)= get_raw_mean_std(out_comp{idx_batch},temp_index_h,temp_index_h,'abort_code','IN');
        out_hnd_abort.(hand)(idx_batch)= get_hand_switch_errors(idx.incorrect_hand{idx_batch},temp_index_h,idx,idx_batch);
        out_dur_fix.(hand)(idx_batch)= get_raw_mean_std(out_comp{idx_batch},temp_index_h,temp_index_h,'dur_fix','IN');
    end
    
    % runs sessions trials
    [out_comp{idx_batch}.saccades.run]=out_comp{idx_batch}.selected.run;
    [out_comp{idx_batch}.saccades.session]=out_comp{idx_batch}.selected.session;
    [out_comp{idx_batch}.saccades.trial]=out_comp{idx_batch}.selected.trials;
    [out_comp{idx_batch}.reaches.run]=out_comp{idx_batch}.selected.run;
    [out_comp{idx_batch}.reaches.session]=out_comp{idx_batch}.selected.session;
    [out_comp{idx_batch}.reaches.trial]=out_comp{idx_batch}.selected.trials;
end

%% loops
for t = 1:numel(type_names)
    for e = 1:numel(effector_names)
        batches_to_look_at  = find(type_effector_indicator{t,e});
        effector            =  effector_names{e};
        type                =  type_names{t};
        type_effector = ['t_' type '_e_' effector];
        for s_r = 1:numel(sac_rea)
            reach_or_saccade = sac_rea{s_r};
            for b = 1:numel(batches_to_look_at)
                batch = batches_to_look_at(b);
                temp_index_b = [out_comp{batch}.task.effector] == str2double(effector) & [out_comp{batch}.task.type] == str2double(type);
                batch_correlation_conditions = [out_comp{batch}.correlation.conditions];
                temp_index_bcorr = [batch_correlation_conditions.effector] == str2double(effector) & [batch_correlation_conditions.type] ==  str2double(type) & [batch_correlation_conditions.success] ==true;
                for d = 1:numel(decision_names)
                    decision = decision_names{d};
                    temp_index_d = temp_index_b & idx.(decision){batch};
                    for p = 1:numel(unique_pos.(reach_or_saccade))
                        
                        P_index_s   = abs([out_comp{batch}.(reach_or_saccade).tar_pos]-[out_comp{batch}.(reach_or_saccade).fix_pos]-unique_pos.(reach_or_saccade)(p)) <1.5 & temp_index_d & idx.success{batch};
                        P_index_a   = abs([out_comp{batch}.(reach_or_saccade).tar_pos]-[out_comp{batch}.(reach_or_saccade).fix_pos]-unique_pos.(reach_or_saccade)(p)) <1.5 & temp_index_d & ~idx.success{batch};
                        P_index     = abs([out_comp{batch}.(reach_or_saccade).tar_pos]-[out_comp{batch}.(reach_or_saccade).fix_pos]-unique_pos.(reach_or_saccade)(p)) <1.5 & temp_index_d;
                        P_index_t_a   = abs([out_comp{batch}.(reach_or_saccade).tar_pos]-[out_comp{batch}.(reach_or_saccade).fix_pos]-unique_pos.(reach_or_saccade)(p)) <1.5 & temp_index_d & idx.hnd_abort_tar_acq{batch};
                        
                        
                        
                        
                        if strcmp(reach_or_saccade,'reaches')
                            out_str(batch).(type_effector).(reach_or_saccade).(decision).endpoints_per_position_t_a(p,1)    = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_t_a).endpos]      - [out_comp{batch}.(reach_or_saccade)(P_index_t_a).fix_pos]);
                            out_str(batch).(type_effector).(reach_or_saccade).(decision).endpoints_per_position_s(p,1)  = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_s).endpos]    - [out_comp{batch}.(reach_or_saccade)(P_index_s).fix_pos]);
                            out_str(batch).(type_effector).(reach_or_saccade).(decision).endpoints_per_position_a(p,1)  = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_a).endpos]    - [out_comp{batch}.(reach_or_saccade)(P_index_a).fix_pos]);
                            out_str(batch).(type_effector).(reach_or_saccade).(decision).endpoints_per_position(p,1)    = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index).endpos]      - [out_comp{batch}.(reach_or_saccade)(P_index).fix_pos]);
                        else
                            
                            out_str(batch).(type_effector).(reach_or_saccade).(decision).endpoints_per_position_s(p,1)  = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_s).endpos]);
                            out_str(batch).(type_effector).(reach_or_saccade).(decision).endpoints_per_position_a(p,1)  = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_a).endpos]);
                            out_str(batch).(type_effector).(reach_or_saccade).(decision).endpoints_per_position(p,1)    = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index).endpos]);
                        end
                        %                         out_str(batch).(type_effector).(reach_or_saccade).(decision).accuracy_xy(p,1)            = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index).accuracy_xy]);
                        
                    end
                    if GLO.keep_raw_output
                        abort_raw_index = out_comp{batch}.abort_raw_index;
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).abort_raw_x           = {out_comp{batch}.(reach_or_saccade)(temp_index_d & abort_raw_index).raw_x};
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).abort_raw_y           = {out_comp{batch}.(reach_or_saccade)(temp_index_d & abort_raw_index).raw_y};
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).abort_raw_states      = {out_comp{batch}.(reach_or_saccade)(temp_index_d & abort_raw_index).raw_states};
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).abort_raw_time_axis   = {out_comp{batch}.(reach_or_saccade)(temp_index_d & abort_raw_index).raw_time_axis};
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).abort_fix_pos         = [out_comp{batch}.(reach_or_saccade)(temp_index_d & abort_raw_index).fix_pos];
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).abort_tar_pos         = [out_comp{batch}.(reach_or_saccade)(temp_index_d & abort_raw_index).tar_pos];
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).abort_lat             = [out_comp{batch}.(reach_or_saccade)(temp_index_d & abort_raw_index).lat];
                        
                        success_raw_index = out_comp{batch}.success_raw_index;
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).success_raw_x           = {out_comp{batch}.(reach_or_saccade)(temp_index_d & success_raw_index).raw_x};
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).success_raw_y           = {out_comp{batch}.(reach_or_saccade)(temp_index_d & success_raw_index).raw_y};
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).success_raw_states      = {out_comp{batch}.(reach_or_saccade)(temp_index_d & success_raw_index).raw_states};
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).success_raw_time_axis   = {out_comp{batch}.(reach_or_saccade)(temp_index_d & success_raw_index).raw_time_axis};
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).success_fix_pos         = [out_comp{batch}.(reach_or_saccade)(temp_index_d & success_raw_index).fix_pos];
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).success_tar_pos         = [out_comp{batch}.(reach_or_saccade)(temp_index_d & success_raw_index).tar_pos];
                        out_str(batch).(type_effector).(reach_or_saccade).(decision).success_lat             = [out_comp{batch}.(reach_or_saccade)(temp_index_d & success_raw_index).lat];
                    end
                    for s = 1:numel(side_names)
                        side = side_names{s};
                        counterside = side_names{mod(s,2)+1};
                        s_index = temp_index_d & idx.([side '_tar']){batch};
                        cs_index = temp_index_d & idx.([counterside '_tar']){batch};
                        out_str(batch).(type_effector).(reach_or_saccade).([side '_' decision]) = get_raw_mean_std(out_comp{batch},s_index,cs_index,reach_or_saccade,decision);
                        out_str(batch).(type_effector).(reach_or_saccade).([side '_' decision]).abort_code = [{out_comp{batch}.task(s_index).abort_code}];
                        
                    end
                    if str2num(effector)==0
                        continue
                    end
                    for h = 1:numel(hand_names)
                        hand = hand_names{h};
                        temp_index_h = temp_index_d & idx.(['reach_' hand]){batch};
                        if GLO.keep_raw_output
                            abort_raw_index = out_comp{batch}.abort_raw_index;
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).abort_raw_x           = {out_comp{batch}.(reach_or_saccade)(temp_index_h & abort_raw_index).raw_x};
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).abort_raw_y           = {out_comp{batch}.(reach_or_saccade)(temp_index_h & abort_raw_index).raw_y};
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).abort_raw_states      = {out_comp{batch}.(reach_or_saccade)(temp_index_h & abort_raw_index).raw_states};
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).abort_raw_time_axis   = {out_comp{batch}.(reach_or_saccade)(temp_index_h & abort_raw_index).raw_time_axis};
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).abort_fix_pos         = [out_comp{batch}.(reach_or_saccade)(temp_index_h & abort_raw_index).fix_pos];
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).abort_tar_pos         = [out_comp{batch}.(reach_or_saccade)(temp_index_h & abort_raw_index).tar_pos];
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).abort_lat             = [out_comp{batch}.(reach_or_saccade)(temp_index_h & abort_raw_index).lat];
                            
                            success_raw_index = out_comp{batch}.success_raw_index;
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).success_raw_x           = {out_comp{batch}.(reach_or_saccade)(temp_index_h & success_raw_index).raw_x};
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).success_raw_y           = {out_comp{batch}.(reach_or_saccade)(temp_index_h & success_raw_index).raw_y};
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).success_raw_states      = {out_comp{batch}.(reach_or_saccade)(temp_index_h & success_raw_index).raw_states};
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).success_raw_time_axis   = {out_comp{batch}.(reach_or_saccade)(temp_index_h & success_raw_index).raw_time_axis};
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).success_fix_pos         = [out_comp{batch}.(reach_or_saccade)(temp_index_h & success_raw_index).fix_pos];
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).success_tar_pos         = [out_comp{batch}.(reach_or_saccade)(temp_index_h & success_raw_index).tar_pos];
                            out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).success_lat             = [out_comp{batch}.(reach_or_saccade)(temp_index_h & success_raw_index).lat];
                        end
                        
                        
                        for p = 1:numel(unique_pos.(reach_or_saccade))
                            
                            P_index_s   = abs([out_comp{batch}.(reach_or_saccade).tar_pos]-[out_comp{batch}.(reach_or_saccade).fix_pos]-unique_pos.(reach_or_saccade)(p)) <1.5 & temp_index_h & idx.success{batch};
                            P_index_a   = abs([out_comp{batch}.(reach_or_saccade).tar_pos]-[out_comp{batch}.(reach_or_saccade).fix_pos]-unique_pos.(reach_or_saccade)(p)) <1.5 & temp_index_h & ~idx.success{batch};
                            P_index     = abs([out_comp{batch}.(reach_or_saccade).tar_pos]-[out_comp{batch}.(reach_or_saccade).fix_pos]-unique_pos.(reach_or_saccade)(p)) <1.5 & temp_index_h;
                            P_index_t_a   = abs([out_comp{batch}.(reach_or_saccade).tar_pos]-[out_comp{batch}.(reach_or_saccade).fix_pos]-unique_pos.(reach_or_saccade)(p)) <1.5 & temp_index_h & idx.hnd_abort_tar_acq{batch};
                            
                            if strcmp(reach_or_saccade,'reaches')
                                out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).endpoints_per_position_t_a(p,1)    = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_t_a).endpos]      - [out_comp{batch}.(reach_or_saccade)(P_index_t_a).fix_pos]);
                                out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).endpoints_per_position_s(p,1)  = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_s).endpos]    - [out_comp{batch}.(reach_or_saccade)(P_index_s).fix_pos]);
                                out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).endpoints_per_position_a(p,1)  = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_a).endpos]    - [out_comp{batch}.(reach_or_saccade)(P_index_a).fix_pos]);
                                out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).endpoints_per_position(p,1)    = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index).endpos]      - [out_comp{batch}.(reach_or_saccade)(P_index).fix_pos]);
                            else
                                
                                out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).endpoints_per_position_s(p,1)  = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_s).endpos]);
                                out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).endpoints_per_position_a(p,1)  = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index_a).endpos]);
                                out_str(batch).(type_effector).(reach_or_saccade).([decision '_' hand]).endpoints_per_position(p,1)    = get_raw_mean_std_xy([out_comp{batch}.(reach_or_saccade)(P_index).endpos]);
                            end
                            
                        end
                        for s = 1:numel(side_names)
                            side = side_names{s};
                            counterside = side_names{mod(s,2)+1};
                            temp_index_s = temp_index_h & idx.([side '_tar']){batch};
                            temp_index_cs = temp_index_h & idx.([counterside '_tar']){batch};
                            temp_index_scorr = temp_index_bcorr & corr_idx.([side '_tar']){batch} & corr_idx.(decision){batch}  & corr_idx.(['reach_' hand]){batch};
                            
                            struct_part_1 = get_raw_mean_std(out_comp{batch},temp_index_s,temp_index_cs,reach_or_saccade,decision);
                            out_str(batch).(type_effector).(reach_or_saccade).([side '_' decision '_' hand]) = struct_part_1;
                            
                            out_str(batch).(type_effector).(reach_or_saccade).([side '_' decision '_' hand]) = get_corr(out_comp{batch}.correlation,temp_index_scorr,struct_part_1);
                            out_str(batch).(type_effector).(reach_or_saccade).([side '_' decision '_' hand]).abort_code            = [{out_comp{batch}.task(temp_index_s).abort_code}];
                        end
                    end
                end
            end
        end
    end
end
end

function [out_stru_ext] = external_cal(out_str)

%% MEANS, STDS AND RAW DATA PER RUN

subparameters           = {'mean','raw','std','num_hits'};
type_effector_fieldnames=fieldnames(out_str);
for te=1:numel(type_effector_fieldnames)
    type_effector=type_effector_fieldnames{te};
    valid_batches= find(cellfun(@(x) ~isempty(x),{out_str.(type_effector)}));
    out_stru_te=[out_str(valid_batches).(type_effector)];
    if isempty(out_stru_te)
        continue;
    end
    reaches_saccades=fieldnames(out_stru_te);
    for rs=1:numel(reaches_saccades)
        sac_rea=reaches_saccades{rs};
        out_stru_te_sr=[out_stru_te.(sac_rea)];
        conditions=fieldnames(out_stru_te_sr);
        for c=1:numel(conditions)
            condition=conditions{c};
            out_stru_te_sr_con=[out_stru_te_sr.(condition)];
            subconditions=fieldnames(out_stru_te_sr_con);
            for sc=1:numel(subconditions)
                subcondition=subconditions{sc};
                out_sc=[out_stru_te_sr_con.(subcondition)];
                if strcmp(subcondition,'endpoints_per_position') || strcmp(subcondition,'endpoints_per_position_s') || strcmp(subcondition,'endpoints_per_position_a')
                    for p=1:size(out_sc,1)
                        temp=get_external_means_std(out_sc(p,:),subparameters);
                        out_stru_ext.(sac_rea).(subcondition)(p).([type_effector '_' condition])=temp;
                    end
                elseif (any(ismember({'abort_code','abort_raw_x','abort_raw_y','abort_raw_states','abort_raw_time_axis','abort_fix_pos','abort_tar_pos','abort_lat','success_raw_x','success_raw_y','success_raw_states','success_raw_time_axis','success_fix_pos','success_tar_pos', 'success_lat'},subcondition)))
                    
                    out_stru_ext.(sac_rea).(subcondition).([type_effector '_' condition])=out_sc;
                else
                    temp=get_external_means_std(out_sc,subparameters);
                    out_stru_ext.(sac_rea).(subcondition).([type_effector '_' condition])=[temp];
                end
            end
        end
    end
end

end



function unique_pos=unique_positions(saccadepositions,target_pos_precision)
target_positions                    =   unique(saccadepositions(~isnan(saccadepositions)));
n_targets                           =   numel(target_positions);
for t = 1:n_targets
    target_positions(abs(target_positions-target_positions(t))<target_pos_precision) = target_positions(t);
end
unique_pos                          = unique(target_positions);
end


% get_raw_mean_std(out_comp{idx_batch},temp_index_h,temp_index_h,'abort_code','IN');
function out = get_hand_switch_errors(temp_index_all,temp_index_h,idx,idx_batch)%,temp_index_h,idx,idx_batch,{'hnd_switch','hnd_stay'}
parameters={'hnd_switch','hnd_stay'};
for p = 1:numel(parameters)
    par                             = parameters{p};
    tmp_idx                         = temp_index_h & idx.(par){idx_batch};
    variable_of_interest            = tmp_idx(temp_index_all);
    out.(par).mean                  = nanmean(variable_of_interest);
    out.(par).raw                   = variable_of_interest;
    out.(par).std                   = nanstd(variable_of_interest);
    out.(par).num_hits              = sum(~isnan(variable_of_interest));
    
end
variable_of_interest            = 1*(temp_index_all);
out.successful.mean             = nanmean(variable_of_interest);
out.successful.raw              = variable_of_interest;
out.successful.std              = nanstd(variable_of_interest);
out.successful.num_hits         = sum(variable_of_interest);
end




function out = get_raw_mean_std(input,temp_index,counterside_index,reach_or_saccade,decision)
global GLO

Abort_codes_temp= {
    'DUMMY', 1;
    'ABORT_USE_INCORRECT_HAND' ,2};
%      'ABORT_HND_FIX_ACQ_STATE'  ,3;
%      'ABORT_HND_FIX_HOLD_STATE' ,4 ...
%     };

temp_index_all=temp_index;
temp_index=temp_index & ([input.states.state_abo] > 3 | [input.states.state_abo] == -1);
temp_cindex=counterside_index & ([input.states.state_abo] > 3 | [input.states.state_abo] == -1);
temp_index_su=temp_index & [input.binary.success];
temp_cindex_su=temp_cindex & [input.binary.success];
if strcmp(reach_or_saccade,'reaches')
    parameters = {'lat','dur','endpos','tar_pos','run','session','trial','accuracy_xy','accuracy_x','accuracy_y'};
elseif strcmp(reach_or_saccade,'saccades')
    parameters = {'lat','dur','endpos','tar_pos','velocity','run','session','trial','accuracy_xy','accuracy_x','accuracy_y'};
else
    %     parameters = {'ini_fix','abort_code'};reach_or_saccade='reaches';
    parameters = {'ini_fix','abort_code','dur_fix'};reach_or_saccade='reaches';%MP
end
for p = 1:numel(parameters)
    par                             = parameters{p};
    if strcmp(par,'abort_code')
        [~, variable_of_interest]   = ismember({input.(reach_or_saccade)(temp_index_all).(par)},Abort_codes_temp(:,1));
    elseif strcmp(par,'accuracy_xy')
        variable_of_interest        = abs(1*([input.(reach_or_saccade)(temp_index_su).(par)]));%% abs to calculate euclidean distances before averaging, otherwise averaging happens to signed x and y complex and then we get the absolute
    elseif strcmp(par,'accuracy_x')
        variable_of_interest        = real(1*([input.(reach_or_saccade)(temp_index_su).accuracy_xy]));
    elseif strcmp(par,'accuracy_y')
        variable_of_interest        = imag(1*([input.(reach_or_saccade)(temp_index_su).accuracy_xy]));
    else
        variable_of_interest        = 1*([input.(reach_or_saccade)(temp_index_su).(par)]);
    end
    
    out.(par).mean                  = nanmean(variable_of_interest);
    out.(par).raw                   = variable_of_interest;
    out.(par).std                   = nanstd(variable_of_interest);
    out.(par).num_hits              = sum(~isnan(variable_of_interest));
    
    % Accuracy as absolute of x and y accuracy means... !!!
    if GLO.accuracy_as_absolute && strcmp(par,'accuracy_xy')
        variable_of_interest            = 1*([input.(reach_or_saccade)(temp_index_su).(par)]);
        out.(par).mean                  = abs(nanmean(variable_of_interest));
        out.(par).raw                   = abs(variable_of_interest);
        out.(par).std                   = abs(nanstd(variable_of_interest));
        out.(par).num_hits              = sum(~isnan(variable_of_interest));
    end
    
end
% if strcmp(decision,'IN')
variable_of_interest            = 1*(temp_index_su(temp_index));
% else
%     variable_of_interest            = [ones(1,sum(temp_index_su)) zeros(1,sum(temp_cindex_su))];
% end
out.successful.mean             = nanmean(variable_of_interest);
out.successful.raw              = variable_of_interest;
out.successful.std              = nanstd(variable_of_interest);
out.successful.num_hits         = sum(variable_of_interest);



if GLO.only_successful_side_selection
    variable_of_interest_side            = [ones(1,sum(temp_index_su)) zeros(1,sum(temp_cindex_su))];
else
    variable_of_interest_side            = [ones(1,sum(temp_index)) zeros(1,sum(temp_cindex))];
end


out.side_selection.mean  = nanmean(variable_of_interest_side);
out.side_selection.raw              = variable_of_interest_side;
out.side_selection.std              = nanstd(variable_of_interest_side);
out.side_selection.num_hits         = sum(variable_of_interest_side);
end


function out = get_corr(input,temp_index,out)
parameters={'lat_r','lat_p','lat_slo','lat_int','lat_r_residuals','lat_p_residuals','lat_slo_residuals','lat_int_residuals','lat_raw_sac_rea','lat_difference_sac_rea','lat_residuals_sac_rea'};
for p=1:numel(parameters)
    par=parameters{p};
    variable_of_interest          = 1*(([input(temp_index).(par)]));
    %% redundand, yes, but who cares?
    out.(par).mean                = nanmean(variable_of_interest);
    out.(par).raw                 = variable_of_interest;
    out.(par).std                 = nanstd(variable_of_interest);
    out.(par).num_hits            = sum(~isnan(variable_of_interest));
end

end

function out = get_raw_mean_std_xy(input)
out.raw      = 1*(input);
out.mean     = nanmean(input);
out.std      = nanstd(real(input)) + 1i*nanstd(imag(input));
out.num_hits = sum(~isnan(input));
end

function out=get_external_means_std(input,subparameters)
for sp=1:numel(subparameters)
    subparameter=subparameters{sp};
    out.(['mean_of_' subparameter])       = double(nanmean([input.(subparameter)]));
    out.(['raw_of_' subparameter])        = double([input.(subparameter)]);
    out.(['sem_of_' subparameter])        = double(nansem(real([input.(subparameter)]))+1i*nansem(imag([input.(subparameter)])));
    out.(['std_of_' subparameter])        = double(nanstd(real([input.(subparameter)]))+1i*nanstd(imag([input.(subparameter)])));
    out.(['num_hits_of_' subparameter])   = double(sum(~isnan([input.(subparameter)])));
end
end

function y=nansem(x)
if all(isnan(x))
    y=NaN;
else
    nonanx=x(~isnan(x));
    y=std(nonanx) / sqrt(length(nonanx));
end
end