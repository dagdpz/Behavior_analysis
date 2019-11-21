function [out_stat] = beh_statistics(batch,testing)
global GLO
fn_g = fieldnames(batch.out_stru_ext);
for i = 1:numel(fn_g)
    Group(1,i)=batch.out_stru_ext.(fn_g{i});
    monkey{1,i} = fn_g{i};
end
batch.GLO.monkey = monkey;
fn_p = fieldnames(batch.unique_pos);
for i = 1:numel(fn_p)
    Positions(1,i)=batch.unique_pos.(fn_p{i});
end

saccades_reaches=fieldnames(Group);
comparissons={'space','hand','choice','effector'};
for sr=1:numel(saccades_reaches)
    sac_rea=saccades_reaches{sr};
    valid_batches= find(cellfun(@(x) ~isempty(x),{Group.(sac_rea)}));
    out_stru_te=[Group(valid_batches).(sac_rea)];
    if isempty(out_stru_te)
        continue;
    end
    parameters=fieldnames(out_stru_te);
    
    for co=1:numel(comparissons)
        comparisson=comparissons{co};
        for pa=1:numel(parameters)
            par=parameters{pa}; %the next if condition defines the parameters to NOT calculate statistics 
            if strcmp(par,'ini_fix') || strcmp(par,'ini_dur') || strcmp(par,'ini_abort') || strcmp(par,'hnd_switch') || strcmp(par,'abort_raw_states') || strcmp(par,'abort_raw_time_axis') || strcmp(par,'abort_raw_x') || strcmp(par,'abort_raw_y') || strcmp(par,'abort_fix_pos') || strcmp(par,'abort_tar_pos') || strcmp(par,'success_raw_states') || strcmp(par,'success_raw_time_axis') || strcmp(par,'success_raw_x') || strcmp(par,'success_raw_y') || strcmp(par,'success_fix_pos') || strcmp(par,'success_tar_pos') || strcmp(par,'abort_code') || strcmp(par,'abort_lat') || strcmp(par,'success_lat') ...
                    || strcmp(par,'abort_trial') || strcmp(par,'abort_run') || strcmp(par,'abort_session') || strcmp(par,'success_trial') || strcmp(par,'success_run') || strcmp(par,'success_session') || strcmp(par,'trial') || strcmp(par,'run') || strcmp(par,'session'), continue, end
            out_stru_te_sr=vertcat(out_stru_te.(par));
            type_effector_fieldnames=fieldnames(out_stru_te_sr);            
            for te=1:numel(type_effector_fieldnames)
                type_effector=type_effector_fieldnames{te};                
                type_effector_counter=get_counter(comparisson,type_effector);
                for g=1:size(Group,2)
                    if strcmp(par,'endpoints_per_position') || strcmp(par,'endpoints_per_position_s') || strcmp(par,'endpoints_per_position_a') || strcmp(par,'endpoints_per_position_t_a') || strcmp(par,'endpos') || strcmp(par,'tar_pos')%% endpos/tar_pos?
                        for p=1:size(out_stru_te_sr,2)
                            if ~isstruct(out_stru_te_sr(1,p).(type_effector))  ||  ~isstruct(out_stru_te_sr(2,p).(type_effector))
                                out_stat.(comparisson)(g).(sac_rea).(par)(p).(type_effector).(subcondition)={NaN NaN NaN NaN};
                                continue
                            end
                            subconditions=fieldnames(out_stru_te_sr(g,p).(type_effector));
                            for sc=1:numel(subconditions)
                                subcondition=subconditions{sc};
                                x=[]; y=[]; hr=NaN; pr=NaN;  hi=NaN; pi=NaN; tab_r=NaN; tab_i=NaN; 
                                if isfield(out_stru_te_sr(g,p),type_effector) && isfield(out_stru_te_sr(g,p),type_effector_counter)
                                x = [out_stru_te_sr(g,p).(type_effector).(subcondition)]';
                                y = [out_stru_te_sr(g,p).(type_effector_counter).(subcondition)]';
                                if ~any(isempty(real(x)) | isempty(real(y)))
                                    if ~all(isnan(real(x))) && ~all(isnan(real(y)))
                                        if GLO.calculate_statististics && ~GLO.one_subject
                                            %if ~strcmp(par,'sac_rea') %| ~strcmp(par,'accuracy')
                                            if (strcmp(testing,'groups')) && ~GLO.parametric_testing
                                                [pr hr tab_r]= ranksum(real(x),real(y));
                                                [pi hi tab_i]= ranksum(imag(x),imag(y));
                                            elseif (strcmp(testing,'groups')) && GLO.parametric_testing
                                                [hr pr cir tab_r]= ttest2(real(x),real(y),0.05,0,1,1);
                                                [hi pi cii tab_i]= ttest2(imag(x),imag(y),0.05,0,1,1);
                                            elseif (strcmp(testing,'patient')) && GLO.parametric_testing
                                                [hr pr cir tab_r]= ttest2(real(x),real(y),0.05,0,1,1);
                                                [hi pi cii tab_i]= ttest2(imag(x),imag(y),0.05,0,1,1);
                                            elseif (strcmp(testing,'patient')) && ~GLO.parametric_testing
                                                [pr hr tab_r]= ranksum(real(x),real(y));
                                                [pi hi tab_i]= ranksum(imag(x),imag(y));
                                            end
                                        else
                                            hr      = 0;
                                            pr      = 1;
                                            hi      = 0;
                                            pi      = 1;
                                            tab_r   = [];
                                            tab_i   = [];
                                        end
                                    end
                                end
                                end
                                out_sc= {x y hr+1i*hi pr+1i*pi tab_r tab_i};
                                out_stat.(comparisson)(g).(sac_rea).(par)(p).(type_effector).(subcondition)=out_sc;
                            end
                        end
                    else
                        
                        out_stru_te_sr_con=[out_stru_te_sr.(type_effector)];
                        subconditions=fieldnames(out_stru_te_sr_con);
                        for sc=1:numel(subconditions)
                            subcondition=subconditions{sc};
                            x=[]; y=[]; h=NaN; pp=NaN; tab=NaN;
                            if isfield(out_stru_te_sr(g,1),type_effector) && isfield(out_stru_te_sr(g,1),type_effector_counter)
                            x = [out_stru_te_sr(g,1).(type_effector).(subcondition)]';
                            y = [out_stru_te_sr(g,1).(type_effector_counter).(subcondition)]';
                            if ~any(isempty(x) | isempty(y))
                                if (~all(isnan(x)) & ~all(isnan(y)))
                                    if GLO.calculate_statististics && ~GLO.one_subject
                                        if isempty(strfind(par,'sac_rea')) 
                                            if (strcmp(par,'successful') && strcmp(subcondition,'raw_of_raw'))
                                                x_n=sum(x==0);
                                                x_s=sum(x==1);
                                                y_n=sum(y==0);
                                                y_s=sum(y==1);
                                                [h pp]= f_exakt_scalars(x_s,x_n,y_s,y_n);
                                            else
                                                if (strcmp(testing,'groups')) && ~GLO.parametric_testing
                                                    [pp h tab]= ranksum(x,y);
                                                elseif (strcmp(testing,'groups')) && GLO.parametric_testing
                                                    [h pp ci tab]= ttest2(x,y,0.05,0,1,1);
                                                elseif (strcmp(testing,'patient')) && GLO.parametric_testing
                                                    [h pp ci tab]= ttest2(x,y,0.05,0,1,1);
                                                elseif (strcmp(testing,'patient')) && ~GLO.parametric_testing
                                                    [pp h tab]= ranksum(x,y);
                                                end
                                            end
                                        end
                                    else
                                        pp= 1; h=0; tab=[];
                                    end
                                end
                            end
                            end
                            out_sc= {x y h pp tab};
                            out_stat.(comparisson)(g).(sac_rea).(par)().(type_effector).(subcondition)=out_sc;
                        end
                    end
                end
                
            end
        end
    end
end
end


function type_effector_counter=get_counter(comparisson,type_effector) % this function find the opposite condition (defined by comparison)
global GLO

type_effector_counter=type_effector;
switch comparisson
    case 'space'
        L_idx=strfind(type_effector,'L_');
        R_idx=strfind(type_effector,'R_');
        if ~isempty(L_idx)
            type_effector_counter(L_idx)='R';
        elseif ~isempty(R_idx)
            type_effector_counter(R_idx)='L';
        end
    case 'hand'
        L_idx=strfind(type_effector,'LH');
        R_idx=strfind(type_effector,'RH');
        if ~isempty(L_idx)
            type_effector_counter(L_idx)='R';
        elseif ~isempty(R_idx)
            type_effector_counter(R_idx)='L';
        end
    case 'choice'
        CH_idx=strfind(type_effector,'CH');
        IN_idx=strfind(type_effector,'IN');
        if ~isempty(CH_idx)
            type_effector_counter(CH_idx:CH_idx+1)='IN';
        elseif ~isempty(IN_idx)
            type_effector_counter(IN_idx:IN_idx+1)='CH';
        end
    case 'effector'
        %% so far, only comparing effector 6 with effector 4
        e_4_idx=strfind(type_effector,'e_4');
        e_fg_idx=strfind(type_effector,['e_' GLO.type_of_free_gaze]);
        if ~isempty(e_4_idx)
            type_effector_counter(e_4_idx+2)=GLO.type_of_free_gaze;
        elseif ~isempty(e_fg_idx)
            type_effector_counter(e_fg_idx+2)='4';
        end
end
end

function [h P_fexakt]=f_exakt_scalars(n_par1_grp1,n_par0_grp1,n_par1_grp2,n_par0_grp2)
n_pars_group1_binary=[ones(n_par1_grp1,1); zeros(n_par0_grp1,1)];
n_pars_group2_binary=[ones(n_par1_grp2,1); zeros(n_par0_grp2,1)];
Status=[ones(numel(n_pars_group1_binary),1); zeros(numel(n_pars_group2_binary),1)];
Binary_data=[n_pars_group1_binary; n_pars_group2_binary];
P_fexakt = fexact(Binary_data,Status);
if P_fexakt >= 0.05
    h = 0;
else
    h = 1;
end
end
