function [anova_out] = bh_anova_3_factors(batch,testing)
global GLO

% fn_g = fieldnames(batch.out_stru_ext);
fn_g = {'Experimental'};
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
    %     parameters=fieldnames(out_stru_te);
    parameters={'accuracy_y';'lat';'dur'};
    
    for pa=1:numel(parameters)
        par=parameters{pa};
        out_stru_te_sr=vertcat(out_stru_te.(par));
        type_effector_fieldnames=fieldnames(out_stru_te_sr);
        anova_type_effectors = {'t_4_e_4';'t_4_e_6'};
        t_e_A=anova_type_effectors{1};
        t_e_B=anova_type_effectors{2};
        
        FN_A.LSLH=[t_e_A '_L_IN_LH'];
        FN_A.LSRH=[t_e_A '_L_IN_RH'];
        FN_A.RSLH=[t_e_A '_R_IN_LH'];
        FN_A.RSRH=[t_e_A '_R_IN_RH'];
        
        FN_B.LSLH=[t_e_B '_L_IN_LH'];
        FN_B.LSRH=[t_e_B '_L_IN_RH'];
        FN_B.RSLH=[t_e_B '_R_IN_LH'];
        FN_B.RSRH=[t_e_B '_R_IN_RH'];
        
        if all(ismember({FN_A.LSLH; FN_A.LSRH; FN_A.RSLH; FN_A.RSRH; FN_B.LSLH; FN_B.LSRH; FN_B.RSLH; FN_B.RSRH},type_effector_fieldnames))
            %                     subconditions=fieldnames(out_stru_te_sr(1).(FN_A.LSLH));
            subconditions={'raw_of_raw'};
            for sc=1:numel(subconditions)
                if ~GLO.one_subject
                    subcondition=subconditions{sc};
                    VA(1).LSLH=out_stru_te_sr(1).(FN_A.LSLH).(subcondition);
                    VA(1).LSRH=out_stru_te_sr(1).(FN_A.LSRH).(subcondition);
                    VA(1).RSLH=out_stru_te_sr(1).(FN_A.RSLH).(subcondition);
                    VA(1).RSRH=out_stru_te_sr(1).(FN_A.RSRH).(subcondition);
                    VA(2).LSLH=out_stru_te_sr(1).(FN_B.LSLH).(subcondition);
                    VA(2).LSRH=out_stru_te_sr(1).(FN_B.LSRH).(subcondition);
                    VA(2).RSLH=out_stru_te_sr(1).(FN_B.RSLH).(subcondition);
                    VA(2).RSRH=out_stru_te_sr(1).(FN_B.RSRH).(subcondition);
                    
                    anovainput=[VA(1).LSLH(:);VA(1).LSRH(:);VA(1).RSLH(:);VA(1).RSRH(:);VA(2).LSLH(:); VA(2).LSRH(:);VA(2).RSLH(:);VA(2).RSRH(:)];
                    factoreffector=[...
                        ones(numel([VA(1).LSLH(:);VA(1).LSRH(:);VA(1).RSLH(:);VA(1).RSRH(:)]),1);...
                        zeros(numel([VA(2).LSLH(:);VA(2).LSRH(:);VA(2).RSLH(:);VA(2).RSRH(:)]),1)];
                    factorspace=[...
                        ones(numel([VA(1).LSLH(:);VA(1).LSRH(:)]),1);...
                        zeros(numel([VA(1).RSLH(:);VA(1).RSRH(:)]),1);...
                        ones(numel([VA(2).LSLH(:);VA(2).LSRH(:)]),1);...
                        zeros(numel([VA(2).RSLH(:);VA(2).RSRH(:)]),1)];
                    factorhand=[...
                        ones(numel(VA(1).LSLH(:)),1);...
                        zeros(numel(VA(1).LSRH(:)),1);...
                        ones(numel(VA(1).RSLH(:)),1);...
                        zeros(numel(VA(1).RSRH(:)),1);...
                        ones(numel(VA(2).LSLH(:)),1);...
                        zeros(numel(VA(2).LSRH(:)),1);...
                        ones(numel(VA(2).RSLH(:)),1);...
                        zeros(numel(VA(2).RSRH(:)),1)];
                    
                    if any(~isnan(anovainput))
                        [anova_out.(sac_rea).(par).(subcondition).preal, anova_out.(sac_rea).(par).(subcondition).tab_real]=anovan(real(anovainput),[factoreffector factorspace factorhand],'varnames',{'effector_type','space','hand'},'display','off','model','full');
                        [anova_out.(sac_rea).(par).(subcondition).pimag, anova_out.(sac_rea).(par).(subcondition).tab_imag]=anovan(imag(anovainput),[factoreffector factorspace factorhand],'varnames',{'effector_type','space','hand'},'display','off','model','full');
                    end
                else
                    [anova_out.(sac_rea).(par).(subcondition).preal, anova_out.(sac_rea).(par).(subcondition).tab_real]=NaN;
                    [anova_out.(sac_rea).(par).(subcondition).pimag, anova_out.(sac_rea).(par).(subcondition).tab_imag]=NaN;
                end
            end
        end
        %             end
    end
end
end
