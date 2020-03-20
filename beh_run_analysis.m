global GLO

if ~isempty(GLO.parent_folder) && GLO.one_subject
    GLO.folder_to_save = [GLO.parent_folder, filesep, group{1}{:}, filesep, num2str(dates_subject_in{1}{:})];
end
if ~exist(GLO.folder_to_save,'dir')
    mkdir(GLO.folder_to_save)
end
if ~isfield(GLO,'accuracy_as_absolute')
   GLO.accuracy_as_absolute=0; 
end

for gr=1:numel(group)  %this entire loop defines how to analyse data (batching, i.e run by run, group by group etc)  
    if GLO.one_subject == 1 && (gr>1)
        continue
    end
    %% !!
    
    dates_subject_temp{gr}                  = {};
    subject_files_temp{gr}                  = {};
    switch batching_type{gr}
        case 1 % run by run
            batching{gr}.runs_as_batches_val = 1;
        case {2,3} % session by session or group by group
            batching{gr}.runs_as_batches_val = 0;
    end
    
    dag_drive='Y:';
    da=[];
    for sub=1:numel(group{gr})        
        if GLO.testing_patient
            glo.tmp_dir                    = [dag_drive filesep 'Data\Patient_and_controls_20160212\' group{gr}{sub}];
        else
            glo.tmp_dir                    = [dag_drive filesep 'Data\' group{gr}{sub}];
        end
        
        [~, date_files{gr}{sub}]                 = beh_arrange_trials_eye_hand(glo.tmp_dir, dates_subject_in{gr}{sub}, batching{gr});
        
        for j=1:size(date_files{gr}{sub},1)
            da(j)                   = str2num(date_files{gr}{sub}{j,1}(end-7:end));
        end
        if ~isempty(da)
            dates_subject_temp{gr}          = [dates_subject_temp{gr}; unique(da)];
            subject_files_temp{gr}               = [subject_files_temp{gr}; date_files{gr}(sub)];            
        else
            continue
        end
        clear da        
    end
    
    if batching_type{gr}==1
        dates_subject{gr}={[dates_subject_temp{gr}{:}]};
        subject_files{gr}={vertcat(subject_files_temp{gr}{:})};
        subject_files_ttemp{gr}{1}={};
        % check if runs are repeated in groups with same sessions for different 'subject' for run_by_run batching
        % check if there can be repeated session/run/subject in a group
        if ~GLO.same_day && ~isempty(batching{gr}.runs)
            for sub=1:numel(group{gr})
                idx_run=ismember([subject_files{gr}{1}{:,2}], [batching{gr}.runs{sub}])' & cellfun(@any,strfind(subject_files{gr}{1}(:,1),num2str(dates_subject_in{gr}{sub})));
                subject_files_ttemp{gr}{1} = [subject_files_ttemp{gr}{1};subject_files{gr}{1}(idx_run,:)];
            end
        end
        subject_files{gr}{1}=uniqueRowsCA(subject_files_ttemp{gr}{1});
    end
    if batching_type{gr}==2
        dates_subjects_temp_gr=[dates_subject_temp{gr}{:}];
        [~,unique_date_index]=unique(dates_subjects_temp_gr);
        dates_subject_tmp=dates_subjects_temp_gr(sort(unique_date_index))';
        subject_files_tmp=vertcat(subject_files_temp{gr}{:});
        subject_files{gr}=cell(numel(dates_subject_tmp),1);
        for d=1:numel(dates_subject_tmp)
            dates_subject{gr}{d}=dates_subject_tmp(d);
            subject_files{gr}{d}=subject_files_tmp(cellfun(@(x) ~isempty(strfind(x,num2str(dates_subject_tmp(d)))),subject_files_tmp(:,1)),:);
        end
    end
    if batching_type{gr}==3
        dates_subject{gr}=dates_subject_temp{gr};
        subject_files{gr}=subject_files_temp{gr};
        if ~GLO.same_day && ~isempty(batching{gr}.runs)
            for sub=1:numel(group{gr})
                idx_run=ismember([subject_files{gr}{sub}{:,2}], [batching{gr}.runs{sub}])';
                subject_files{gr}{sub} = subject_files{gr}{sub}(idx_run,:);
            end
        end
        
        if GLO.same_day && ~isempty(batching{gr}.runs)
            for session=1:size(subject_files{gr},1)
                idx_run=ismember([subject_files{gr}{session}{:,2}], [batching{gr}.runs{session}])';
                subject_files{gr}{session} = subject_files{gr}{session}(idx_run,:);
            end
        end        
    end   
end

for gr=1:numel(group)
    if GLO.one_subject == 1 && (gr>1)
        continue
    end
    %here get behavioral analysis from Monkeypsych_analyse and perform
    %some correlation (latency with hands, etc) everything is saved in the structure bactch 
    [batch.files_for_input.(subject_ID{gr}), batch.out_comp.(subject_ID{gr}), batch.out_stru_ext.(subject_ID{gr}), batch.unique_pos.(subject_ID{gr})] ...
        = beh_reaction_time_analysis(group{gr},dates_subject{gr}, batching{gr},subject_files{gr},steady);
end

%% Here repetition for the same group
if GLO.one_subject == 1
    batch.files_for_input.(subject_ID{2})   = batch.files_for_input.(subject_ID{1});
    batch.out_comp.(subject_ID{2})          = batch.out_comp.(subject_ID{1});
    batch.out_stru_ext.(subject_ID{2})      = batch.out_stru_ext.(subject_ID{1});
    batch.unique_pos.(subject_ID{2})        = batch.unique_pos.(subject_ID{1});
end


fn_g = fieldnames(batch.out_stru_ext);
for gr = 1:numel(fn_g)
    Group_temp(1,gr)=batch.out_stru_ext.(fn_g{gr});
    Positions(1,gr)=batch.unique_pos.(fn_g{gr});
end

unique_saccade_positions=unique([Positions.saccades]);
unique_reach_positions=unique([Positions.reaches]);
fieldnames_per_position={'endpoints_per_position','endpoints_per_position_s','endpoints_per_position_a','endpoints_per_position_t_a' };

for g = 1:numel(Group_temp) % I think this loop just check that the structure field exist, if not create it and fill it with NaN
    for FN=fieldnames_per_position
        % saccades
        if ~isfield(Group_temp(g).saccades,FN{:}) 
            continue
        else
            s_current_fieldnames=fieldnames(Group_temp(g).saccades.(FN{:})(1));
            s_NaNcell=repmat({NaN},numel(s_current_fieldnames),1);
            s_fieldname_dummie=[s_current_fieldnames s_NaNcell]';
            for s_p=1:numel(unique_saccade_positions)
                s_idx=ismember(Positions(g).saccades,unique_saccade_positions(s_p));
                if any(s_idx)
                    batch.out_stru_ext.(subject_ID{g}).saccades.(FN{:})(s_p)=Group_temp(g).saccades.(FN{:})(s_idx);
                else
                    batch.out_stru_ext.(subject_ID{g}).saccades.(FN{:})(s_p)=struct(s_fieldname_dummie{:});
                end
            end
        end
        % reaches
        if ~isfield(Group_temp(g).reaches,FN{:})
            continue          
        else
            r_current_fieldnames=fieldnames(Group_temp(g).reaches.(FN{:})(1));
            r_NaNcell=repmat({NaN},numel(r_current_fieldnames),1);
            r_fieldname_dummie=[r_current_fieldnames r_NaNcell]';           
            for r_p=1:numel(unique_reach_positions)
                r_idx=ismember(Positions(g).reaches,unique_reach_positions(r_p));
                if any(r_idx)
                    batch.out_stru_ext.(subject_ID{g}).reaches.(FN{:})(r_p)=Group_temp(g).reaches.(FN{:})(r_idx);
                else
                    batch.out_stru_ext.(subject_ID{g}).reaches.(FN{:})(r_p)=struct(r_fieldname_dummie{:});
                end
            end
        end
    end
end

if numel(group) > 1
oo=create_combined_nan_structure(batch.out_stru_ext.Control, batch.out_stru_ext.Experimental); %create NaN structure
new.Control=oo;
new.Experimental=oo;
else 
oo=create_combined_nan_structure(batch.out_stru_ext.Control, batch.out_stru_ext.Control);
new.Control=oo;
new.Experimental=oo;
end

fn=fieldnames(batch.out_stru_ext); 
for gr=1:numel(fn)  % here concatinate the existing with NaN to make sure no field is missing anywhere ??
    fnn=fieldnames(batch.out_stru_ext.(fn{gr}));
    for j=1:numel(fnn)
        fnnn= fieldnames(batch.out_stru_ext.(fn{gr}).(fnn{j}));
        for within_sub=1:numel(fnnn)
            new.(fn{gr}).(fnn{j}).(fnnn{within_sub})=catstruct(new.(fn{gr}).(fnn{j}).(fnnn{within_sub}), batch.out_stru_ext.(fn{gr}).(fnn{j}).(fnnn{within_sub}));            
         end
    end
end

clear batch.out_stru_ext
batch.out_stru_ext=new;



if GLO.testing_patient
    testing='patient';
else
    testing='groups';
end


% [anova_out] = bh_anova_3_factors(batch,testing);
[batch.stat]=beh_statistics(batch,testing);
[batch.stat.groups, ~]=beh_statistics_control_exp_group(batch,testing);

save([GLO.folder_to_save filesep 'data'],'-struct','batch');

GLO.correlation_mode = steady.correlation_mode;
if GLO.plot_it          == 1
    beh_compare_groups(batch,testing)
end
