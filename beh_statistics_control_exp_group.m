function [out_stat, anova_out]= beh_statistics_control_exp_group(batch,testing)
global GLO
anova_out=[];
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

for sr=1:numel(saccades_reaches)
    sac_rea=saccades_reaches{sr};
    valid_batches= find(cellfun(@(x) ~isempty(x),{Group.(sac_rea)}));
    out_stru_te=[Group(valid_batches).(sac_rea)];
    if isempty(out_stru_te)
        continue;
    end
    parameters=fieldnames(out_stru_te);
    to_skip={'ini_fix','ini_dur','ini_abort','hnd_switch','abort_raw_states','abort_raw_time_axis','abort_raw_x','abort_raw_y','abort_fix_pos','abort_tar_pos','success_raw_states','success_raw_time_axis',...
    'success_raw_x','success_raw_y','success_fix_pos','success_tar_pos','abort_code','abort_lat','success_lat','completed_raw_states','completed_raw_time_axis','completed_raw_x','completed_raw_y',...
    'completed_fix_pos','completed_tar_pos','completed_lat','abort_trial','abort_run','abort_session','success_trial','success_run','success_session','trial','run','session'};

    per_pos_parameters={'endpoints_per_position','endpoints_per_position_s','endpoints_per_position_a','endpoints_per_position_c','endpoints_per_position_t_a','endpos','tar_pos'};%% endpos/tar_pos?
    for pa=1:numel(parameters)
        par=parameters{pa};
        if  ismember(par,to_skip)
            continue
        end
        out_stru_te_sr=vertcat(out_stru_te.(par));
        type_effector_fieldnames=fieldnames(out_stru_te_sr);
        
        
        %% ANOVAS!!
        %         if ~strcmp(par,'endpoints_per_position') || strcmp(par,'endpoints_per_position_s') || strcmp(par,'endpoints_per_position_a')
        %             anova_type_effectors=unique(cellfun(@(x) x(1:7),type_effector_fieldnames,'uniformoutput',false));
        %             for ate=1:numel(anova_type_effectors)
        %                 type_effector=anova_type_effectors{ate};
        %
        %                 FN.LSLH=[type_effector '_L_IN_LH'];
        %                 FN.LSRH=[type_effector '_L_IN_RH'];
        %                 FN.RSLH=[type_effector '_R_IN_LH'];
        %                 FN.RSRH=[type_effector '_R_IN_RH'];
        %
        %                 if all(ismember({FN.LSLH; FN.LSRH; FN.RSLH; FN.RSRH},type_effector_fieldnames))
        %                     subconditions=fieldnames(out_stru_te_sr(1).(FN.LSLH));
        %                     for sc=1:numel(subconditions)
        %                         if ~GLO.one_subject
        %                             subcondition=subconditions{sc};
        %                             VA(1).LSLH=out_stru_te_sr(1).(FN.LSLH).(subcondition);
        %                             VA(1).LSRH=out_stru_te_sr(1).(FN.LSRH).(subcondition);
        %                             VA(1).RSLH=out_stru_te_sr(1).(FN.RSLH).(subcondition);
        %                             VA(1).RSRH=out_stru_te_sr(1).(FN.RSRH).(subcondition);
        %                             VA(2).LSLH=out_stru_te_sr(2).(FN.LSLH).(subcondition);
        %                             VA(2).LSRH=out_stru_te_sr(2).(FN.LSRH).(subcondition);
        %                             VA(2).RSLH=out_stru_te_sr(2).(FN.RSLH).(subcondition);
        %                             VA(2).RSRH=out_stru_te_sr(2).(FN.RSRH).(subcondition);
        %
        %                             anovainput=[VA(1).LSLH(:);VA(1).LSRH(:);VA(1).RSLH(:);VA(1).RSRH(:);VA(2).LSLH(:); VA(2).LSRH(:);VA(2).RSLH(:);VA(2).RSRH(:)];
        %                             factorgoup=[ones(numel([VA(1).LSLH(:);VA(1).LSRH(:);VA(1).RSLH(:);VA(1).RSRH(:)]),1);...
        %                                 zeros(numel([VA(2).LSLH(:);VA(2).LSRH(:);VA(2).RSLH(:);VA(2).RSRH(:)]),1)];
        %                             factorspace=[ones(numel([VA(1).LSLH(:);VA(1).LSRH(:)]),1);...
        %                                 zeros(numel([VA(1).RSLH(:);VA(1).RSRH(:)]),1);...
        %                                 ones(numel([VA(2).LSLH(:);VA(2).LSRH(:)]),1);...
        %                                 zeros(numel([VA(2).RSLH(:);VA(2).RSRH(:)]),1)];
        %                             factorhand=[ones(numel(VA(1).LSLH(:)),1);zeros(numel(VA(1).LSRH(:)),1);...
        %                                 ones(numel(VA(1).RSLH(:)),1);zeros(numel(VA(1).RSRH(:)),1);...
        %                                 ones(numel(VA(2).LSLH(:)),1);zeros(numel(VA(2).LSRH(:)),1);...
        %                                 ones(numel(VA(2).RSLH(:)),1);zeros(numel(VA(2).RSRH(:)),1)];
        %                             if any(~isnan(anovainput))
        %                                 %if any(iscomplex(anovainput))
        %                                 anova_out.(sac_rea).(par).(type_effector).(subcondition).preal=anovan(real(anovainput),[factorgoup factorspace factorhand],'display','off','model','full');
        %                                 anova_out.(sac_rea).(par).(type_effector).(subcondition).pimag=anovan(imag(anovainput),[factorgoup factorspace factorhand],'display','off','model','full');
        %                                 %                         else
        %                                 %                     anova_out.(sac_rea).(par).(type_effector).(subcondition).preal=anovan(anovainput,[factorgoup factorspace factorhand],'display','off');
        %                                 %                         end
        %
        %                             end
        %                         else
        %                             anova_out.(sac_rea).(par).(type_effector).(subcondition).preal=NaN;
        %                             anova_out.(sac_rea).(par).(type_effector).(subcondition).pimag=NaN;
        %                         end
        %                     end
        %                 end
        %             end
        %         end
        
        for te=1:numel(type_effector_fieldnames)
            type_effector=type_effector_fieldnames{te};
            if ismember(par,per_pos_parameters) || strcmp(par,'endpos') || strcmp(par,'tar_pos') %% endpos? tar_pos?
                for p=1:size(out_stru_te_sr,2)
                    if ~isstruct(out_stru_te_sr(1,p).(type_effector))  |  ~isstruct(out_stru_te_sr(2,p).(type_effector))
                        out_stat.(sac_rea).(par)(p,1).(type_effector).(subcondition)={NaN NaN NaN NaN};
                        continue
                    end
                    out_stru_te_sr_con=[out_stru_te_sr(:,p).(type_effector)];
                    subconditions=fieldnames(out_stru_te_sr_con);
                    for sc=1:numel(subconditions)
                        subcondition=subconditions{sc};
                        x=[]; y=[]; hr=NaN; pr=NaN;  hi=NaN; pi=NaN; tab_r=NaN; tab_i=NaN;
                        x = [out_stru_te_sr_con(1,1).(subcondition)]';
                        y = [out_stru_te_sr_con(1,2).(subcondition)]';
                        x=x(~isnan(x)); y=y(~isnan(y));
                        if ~any(isempty(real(x)) | isempty(real(y)))
                            if ~all(isnan(real(x))) && ~all(isnan(real(y)))
                                if GLO.calculate_statististics && ~GLO.one_subject && numel(x)>2
                                    if isempty(strfind(par,'sac_rea')) %| ~strcmp(par,'accuracy')
                                        if (strcmp(testing,'groups')) && ~GLO.parametric_testing
                                            [pr hr tab_r]= ranksum(real(x),real(y));
                                            [pi hi tab_i]= ranksum(imag(x),imag(y));
                                        elseif (strcmp(testing,'groups')) && GLO.parametric_testing
                                            [hr pr cir tab_r]= ttest2(real(x),real(y),0.05,0,1,1);
                                            [hi pi cii tab_i]= ttest2(imag(x),imag(y),0.05,0,1,1);
                                        elseif (strcmp(testing,'patient')) && (strcmp(subcondition,'raw_of_mean') || strcmp(subcondition,'raw_of_std' ))
                                            [hr pr tab_r]= ttest_1n(real(x),real(y));
                                            [hi pi tab_i]= ttest_1n(imag(x),imag(y));
                                        elseif (strcmp(testing,'patient')) && ~(strcmp(subcondition,'raw_of_mean') || strcmp(subcondition,'raw_of_std' )) && GLO.parametric_testing
                                            [hr pr cir tab_r]= ttest2(real(x),real(y),0.05,0,1,1);
                                            [hi pi cii tab_i]= ttest2(imag(x),imag(y),0.05,0,1,1);
                                        elseif (strcmp(testing,'patient')) && ~(strcmp(subcondition,'raw_of_mean') || strcmp(subcondition,'raw_of_std' )) && ~GLO.parametric_testing
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
                        out_stat.(sac_rea).(par)(p,1).(type_effector).(subcondition)=out_sc;
                    end
                end
            else
                out_stru_te_sr_con=[out_stru_te_sr.(type_effector)];
                subconditions=fieldnames(out_stru_te_sr_con);
                for sc=1:numel(subconditions)
                    subcondition=subconditions{sc};
                    x=[]; y=[]; h=NaN; pp=NaN; tab=NaN;
                    x = [out_stru_te_sr_con(1,1).(subcondition)]';
                    y = [out_stru_te_sr_con(1,2).(subcondition)]';
                    if ~any(isempty(x) | isempty(y))
                        if ~all(isnan(x)) & ~all(isnan(y))
                            if GLO.calculate_statististics && ~GLO.one_subject
                                if isempty(strfind(par,'sac_rea')) %% strcmp(par,'sac_rea')
                                    if strcmp(par,'successful') && strcmp(subcondition,'raw_of_raw')
                                        x_n=sum(x==0);
                                        x_s=sum(x==1);
                                        y_n=sum(y==0);
                                        y_s=sum(y==1);
                                        [h pp]= f_exakt_scalars(x_s,x_n,y_s,y_n);
                                    else
%                                                                                  if (strcmp(testing,'groups')) && ~GLO.parametric_testing
%                                                                                      [pp h tab]= ranksum(x,y);
%                                                                                  elseif (strcmp(testing,'groups')) && GLO.parametric_testing
%                                                                                      [h pp ci tab]= ttest2(x,y,0.05,0,1,1);
                                        if strcmp(testing,'groups') && ~GLO.parametric_testing && strcmp(subcondition,'raw_of_raw') 
                                            [pp h tab]= ranksum(x,y);
                                        elseif (strcmp(testing,'groups')) && ~GLO.parametric_testing && (strcmp(subcondition,'raw_of_mean') || strcmp(subcondition,'raw_of_std' ))
                                            [pp h tab]= signrank(x,y);
                                        elseif (strcmp(testing,'groups')) && GLO.parametric_testing && (strcmp(subcondition,'raw_of_raw')) 
                                            [h pp ci tab]= ttest2(x,y,0.05,0,1,1);
                                        elseif (strcmp(testing,'groups')) && GLO.parametric_testing  && (strcmp(subcondition,'raw_of_mean') || strcmp(subcondition,'raw_of_std' ))
                                            [h pp ci tab]= ttest(x,y,0.05,'both',1);
                                        elseif (strcmp(testing,'patient')) && (strcmp(subcondition,'raw_of_mean') || strcmp(subcondition,'raw_of_std' ))
                                            [h pp tab]= ttest_1n(x,y);
                                        elseif (strcmp(testing,'patient')) && ~(strcmp(subcondition,'raw_of_mean') || strcmp(subcondition,'raw_of_std' )) && GLO.parametric_testing
                                            [h pp ci tab]= ttest2(x,y,0.05,0,1,1);
                                        elseif (strcmp(testing,'patient')) && ~(strcmp(subcondition,'raw_of_mean') || strcmp(subcondition,'raw_of_std' )) && ~GLO.parametric_testing
                                            [pp h tab]= ranksum(x,y);
                                        end
                                    end
                                end
                            else
                                pp= 1; h=0; tab=[];
                            end
                        end
                    end
                    out_sc= {x y h pp tab};
                    out_stat.(sac_rea).(par).(type_effector).(subcondition)=out_sc;
                end
            end
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

