function [eye_x eye_y time] = bh_eye_traces_eye_hand(varargin)
global GLO

errors=1;
batch_title='succesful';
[target_positions, displacement_types] = center_displacement(trial,keys,'sac');


user_drive=getUserName;
switch user_drive
    case 'a.doming'
        GLO.drive                           = 'L';
    case 'lschneider'
        GLO.drive                           = 'Y';
end

GLO.start_date                      = 20140101;
GLO.end_date                        = 20160101;
GLO.type_to_use                     = 3;
GLO.effector_to_use                 = 0;
GLO.reach_hand                      = 2;
GLO.run_analyze                     = 0;
GLO.save_analyze_output             = 0;
GLO.run_internal_calculation        = 0;
GLO.save_batch_for_later            = 0;
GLO.create_pdf                      = 1;
GLO.append_pdfs                     = 0;
GLO.plot_batches_individually       = 0;
GLO.run_baseline_included_in_ttests = 0;
GLO.stat_to_use                     = 'signed rank';%'ttest_bonf'; % ttest paired, ttest unpaired, signed rank
GLO.ttest_text                      = 1;
GLO.RT_vs_bias_condition            = 'Instructed ipsi'; %'Choice contra','Choice ipsi','Instructed contra','Instructed ipsi'
GLO.low_excentricity_threshold      = 15;
GLO.Labels                          = {'RB','B','-200', '-160', '-120','-80','-40','Go','40','80','100','110','120','130','140','150','160','240','-80 Cue', 'Cue', '80 Cue','-80 Go', 'Go','80 Go'};
GLO.All_stim_states                 = {'fix','fix','fix','fix','fix','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','tar_acq','fix','cue','cue','mem','tar_acq_inv','tar_acq_inv'};
GLO.All_stim_windows                = [-0.200 -0.160 -0.120 -0.080 -0.040 0 0.040 0.080 0.100 0.110 0.120 0.130 0.140 0.150 0.160 0.240 -0.080 0 0.080 -0.080 0 0.080];
GLO.windows_to_use                  = 1;
GLO.folder_to_save                  = strcat(GLO.drive,':', filesep, 'microstim_behavior', filesep, monkey, '_summaries', filesep);
GLO.tests                           = {'Bias','farther'}; % Bias || Endpoints %  all || farther || closer
GLO.excentricity                    = GLO.tests{1,2};
GLO.plot_to_show                    = 6;%[1:7,9:17];
GLO.table                           = 0;
GLO.return                          = 0;
GLO.target_pos_precision            = 2;
errors_only_for_eye_trace_plot      = 1;
GLO.plot_eye_traces                 = 0;
% look up for errors_only_for_eye_trace_plot      = 1; or 0;
if GLO.plot_eye_traces
    GLO.Sel_all = [GLO.Sel_all {'keep_raw_data',1}];
end



STATE.INI_TRI       = 1; % initialize trial
STATE.FIX_ACQ       = 2; % fixation acquisition
STATE.FIX_HOL       = 3; % fixation hold
STATE.TAR_ACQ       = 4; % target acquisition
STATE.TAR_HOL       = 5; % target hold
STATE.CUE_ON        = 6; % cue on
STATE.MEM_PER       = 7; % memory period
STATE.DEL_PER       = 8; % delay period
STATE.TAR_ACQ_INV   = 9; % target acquisition invisible
STATE.TAR_HOL_INV   = 10; % target hold invisible
STATE.MAT_ACQ       = 11; % target acquisition in sample to match
STATE.MAT_HOL       = 12; % target acquisition in sample to match
STATE.MAT_ACQ_MSK   = 13; % target acquisition in sample to match
STATE.MAT_HOL_MSK   = 14; % target acquisition in sample to match
STATE.SEN_RET       = 15; % return to sensors for poffenberger
STATE.ABORT         = 19;
STATE.SUCCESS       = 20;
STATE.REWARD        = 21;
STATE.ITI           = 50;
STATE.CLOSE         = 99;


GLO.line_stim                       = GLO.All_stim_windows([GLO.windows_to_use])*1000;
GLO.train_duration                  = 200;
% GLO.Labels([1 2 GLO.windows_to_use+2])
Side_chosen_BL_binary                       =   [];
Side_chosen_stim_binary                     =   [];
num_hits                                    =   [];
num_hits_temp                               =   [];

idx=0;
for k = 1:numel(out_comp)
    if out_comp{k}.emptyflag==0
        idx=idx+1;
        notemptybatches(idx)=k;
    end
end

switch GLO.effector_to_use
    case 0
        sac_or_rea='saccades';
    case 4
        sac_or_rea='reaches';
end

lat_80_to_cue=[];
for k = notemptybatches
    switch GLO.effector_to_use
        case 0
            horizontal_distance = real([out_comp{k,1}.saccades.tar_pos] - [out_comp{k,1}.saccades.fix_pos]);
            vertical_distance   = imag([out_comp{k,1}.saccades.tar_pos] - [out_comp{k,1}.saccades.fix_pos]);
        case 4
            horizontal_distance = real([out_comp{k,1}.reaches.tar_pos] - [out_comp{k,1}.reaches.fix_pos]);
            vertical_distance   = imag([out_comp{k,1}.reaches.tar_pos] - [out_comp{k,1}.reaches.fix_pos]);
    end
    switch GLO.excentricity
        case 'all'
            idx_L{k}                            =   horizontal_distance<0;
            idx_R{k}                            =   horizontal_distance>0;
            
        case 'farther'
            
            idx_L{k}                            =   horizontal_distance<-GLO.low_excentricity_threshold;
            idx_R{k}                            =   horizontal_distance>GLO.low_excentricity_threshold;
        case 'closer'
            
            idx_L{k}                            =   horizontal_distance<0 & horizontal_distance>=-15;
            idx_R{k}                            =   horizontal_distance>0 & horizontal_distance<=15;
    end
    
    idx_direct{k}                               =   [out_comp{k,1}.task.type]==2;
    idx_memory{k}                               =   [out_comp{k,1}.task.type]==3;
    
    
    idx_choice{k}                               =   [out_comp{k,1}.binary.choice]==1;
    idx_instructed{k}                           =   [out_comp{k,1}.binary.choice]==0;
    
    idx_success{k}                              =   [out_comp{k,1}.binary.success]==1;
    idx_error{k}                                =   [out_comp{k,1}.binary.success]==0 & ~ismember({out_comp{k,1}.task.abort_code},'ABORT_JAW') & [out_comp{k,1}.states.state_abo]>2; %!!
    
    idx_stim{k}                                 =   [out_comp{k,1}.binary.microstim]==1;
    idx_baseline{k}                             =   [out_comp{k,1}.binary.microstim]==0;
    idx_baseline_cho{k}                         =   [out_comp{k,1}.binary.microstim]==0 & idx_choice{k};
    idx_baseline_ins{k}                         =   [out_comp{k,1}.binary.microstim]==0 & idx_instructed{k};
    
    idx_baseline_succ{k}                        =   idx_baseline{k}  & idx_success{k};
    idx_baseline_cho_succ{k}                    =   [out_comp{k,1}.binary.microstim]==0 & idx_choice{k}  & idx_success{k};
    idx_baseline_ins_succ{k}                    =   [out_comp{k,1}.binary.microstim]==0 & idx_instructed{k}  & idx_success{k};
    
    idx_baseline_unsucc{k}                      =   idx_baseline{k}  & idx_error{k};
    idx_baseline_cho_unsucc{k}                  =   [out_comp{k,1}.binary.microstim]==0 & idx_choice{k}  & idx_error{k};
    idx_baseline_ins_unsucc{k}                  =   [out_comp{k,1}.binary.microstim]==0 & idx_instructed{k}  & idx_error{k};
    
    idx_stim_at_0{k}                            =   [out_comp{k,1}.task.stim_start]==0 & idx_stim{k};
    idx_stim_in_fix{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.FIX_HOL & idx_stim{k};
    idx_stim_in_acq{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.TAR_ACQ & idx_stim{k};
    idx_stim_in_cue{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.CUE_ON & idx_stim{k};
    idx_stim_in_mem{k}                          =   [out_comp{k,1}.task.stim_state]==STATE.MEM_PER & idx_stim{k};
    idx_stim_in_acq_inv{k}                      =   [out_comp{k,1}.task.stim_state]==STATE.TAR_ACQ_INV & idx_stim{k};
    
    if errors
        idx_base.all=idx_baseline_unsucc;
        idx_base.cho=idx_baseline_cho_unsucc;
        idx_base.ins=idx_baseline_ins_unsucc;
    else
        idx_base.all=idx_baseline_succ;
        idx_base.cho=idx_baseline_cho_succ;
        idx_base.ins=idx_baseline_ins_succ;
    end
    
    fn= fieldnames(idx_base);
    
    
    switch GLO.effector_to_use
        case 0
            tar_positions=[out_comp{k,1}.saccades.tar_pos]-[out_comp{k,1}.saccades.fix_pos];
            nct_positions=[out_comp{k,1}.saccades.nct_pos]-[out_comp{k,1}.saccades.fix_pos];
        case 4
            tar_positions=[out_comp{k,1}.reaches.tar_pos]-[out_comp{k,1}.reaches.fix_pos];
            nct_positions=[out_comp{k,1}.reaches.nct_pos]-[out_comp{k,1}.reaches.fix_pos];
    end
    
    
    for pos=1:numel(target_positions)
        idx_tar_pos{pos}=abs(tar_positions-tar_positions(pos))<=GLO.target_pos_precision;
        for q=1:numel(fn)
            idx_input_baseline=find(idx_base.(fn{q}){k} & idx_tar_pos{pos});
            for idx_trial=1: numel(idx_input_baseline)
                sac_pos{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.saccades(idx_input_baseline(idx_trial)).endpos;
                idx_acq_onset=find((idx_direct{k}(idx_input_baseline(idx_trial)) & [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==4]) |(idx_memory{k}(idx_input_baseline(idx_trial)) &   [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==6]));
                idx_go=find((idx_direct{k}(idx_input_baseline(idx_trial)) & [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==4]) |(idx_memory{k}(idx_input_baseline(idx_trial)) &   [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==9]));
                idx_state=(idx_direct{k}(idx_trial) & ([out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==3] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==4]  | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==5])) |...
                    (idx_memory{k}(idx_trial) & ( [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==4] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==5] | ...
                    [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==6] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==7] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==9] | [out_comp{k,1}.raw(idx_input_baseline(idx_trial)).states==10]));
                eye_x{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.raw(idx_input_baseline(idx_trial)).x_eye(idx_state) ;
                eye_y{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.raw(idx_input_baseline(idx_trial)).y_eye(idx_state);
                idx_acq_onset=idx_acq_onset(1);
                time{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.raw(idx_input_baseline(idx_trial)).time_axis(idx_state)-out_comp{k,1}.raw(idx_input_baseline(idx_trial)).time_axis(idx_acq_onset);
                ini_ms{k,1}.(fn{q}).position(pos).window(1).trial{idx_trial}=out_comp{k,1}.task(idx_input_baseline(idx_trial)).ini_mis;
            end
            if isempty(idx_input_baseline)
                eye_x{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                eye_y{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                time{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                RT{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                ini_ms{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
                sac_pos{k,1}.(fn{q}).position(pos).window(1).trial{1}=NaN;
            end
            
        end
    end
    
    GLO.fontsize=20;
    
    col                                                         = jet(numel(GLO.Labels));
    eye_y_offset                                                = 0;
    col(1,:)                                                    = [0.2 0.2 0.2];
    col(2,:)                                                    = [0.5 0.5 0.5];
    if GLO.type_to_use                 == 3
        col                                                    = [0.2 0.2 0.2;0.5 0.5 0.5;0.859 0.275 0.6; 0.65 0.247 0.6; 0.35 0.3412 0.6471; 0.122 0.255 0.604];
    end

    subplot_assignment=[1:15];
    for q=1:numel(fn)
        plot_1_title='Raw eye traces x versus time';
        summary_1                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_1_title);
        for pos=1:numel(target_positions)
            
            subplot(5,3,subplot_assignment(pos)) % not correctly assigned yet
            hold on
            for idx_shifting_windows = 1
                for idx_trial=1: numel(eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial)
                    if idx_shifting_windows>1 && isnan(ini_ms{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(end)*1000)
                        plot(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}*1000,eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial},'color','k');
                    else
                        plot(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}*1000,eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial},'color',col(idx_shifting_windows+1,:));
                    end
                end
            end
            if all([idx_memory{k}])
                set(gca,'xlim', [-200,2500],'ylim', [-32,32],'FontSize',GLO.fontsize-6);
            else
                set(gca,'xlim', [-200,500],'ylim', [-32,32],'FontSize',GLO.fontsize-6);
            end
            for idx_shifting_windows = 1
                line([GLO.line_stim(idx_shifting_windows), GLO.line_stim(idx_shifting_windows)+GLO.train_duration], [-2-idx_shifting_windows -2-idx_shifting_windows ],'color',col(idx_shifting_windows+2,:),'linewidth',2,'LineStyle', ':')
                line([GLO.line_stim(idx_shifting_windows), GLO.line_stim(idx_shifting_windows)], [0 -2-idx_shifting_windows ],'color',col(idx_shifting_windows+2,:),'linewidth',2,'LineStyle', ':')
                line([GLO.line_stim(idx_shifting_windows)+GLO.train_duration, GLO.line_stim(idx_shifting_windows)+GLO.train_duration], [0 -2-idx_shifting_windows ],'color',col(idx_shifting_windows+2,:),'linewidth',2,'LineStyle', ':')
            end
            xlabel('Time relative to GO [ms]', 'fontsize', GLO.fontsize-4);
            ylabel('Horizontal eye position [deg]', 'fontsize', GLO.fontsize-4);
            title(['Position ' num2str(round(target_positions(pos)))], 'fontsize', GLO.fontsize-2);
        end
        mtit(summary_1,  [plot_1_title, ' ', batch_title, ' ', fn{q}, ' batch ', num2str(k), ' plotting errors? ' num2str(errors) ], 'xoff', -0.0, 'yoff', 0.04, 'color', [0 0 0], 'fontsize', GLO.fontsize,'Interpreter', 'none');
        
        if GLO.create_pdf
            export_fig([GLO.folder_to_save batch_title, ' ', fn{q}, ' ' plot_1_title ' batch ' num2str(k) , ' plotting errors ' num2str(errors) ], '-pdf','-transparent') % pdf by run
            close(gcf)
        end
        
        plot_2_title='Raw eye traces y versus x';
        summary_2                                                   = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_2_title);
        hold on
        for pos=1:numel(target_positions)
            
            angles=[0:pi/100:2*pi];
            circle_x=cos(angles);
            circle_y=sin(angles);
            stepsize_quadrants=pi/2;
            all_phis=[-pi:stepsize_quadrants:pi];
            stepsize_plot=pi/60;
            
            for idx_shifting_windows = 1
                for idx_trial=1: numel(eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial)
                    if idx_shifting_windows==1 && isnan(ini_ms{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(end)*1000)
                        plot(eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}>=0),...
                            eye_y{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}>=0),'color','k');
                        
                    else
                        plot(eye_x{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}>=0),...
                            eye_y{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}(time{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}>=0),'color',col(idx_shifting_windows+1,:));
                    end
                    scatter(real(sac_pos{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}),imag(sac_pos{k,1}.(fn{q}).position(pos).window(idx_shifting_windows).trial{idx_trial}),'ro');
                end
            end

            center=[real(target_positions(pos)) imag(target_positions(pos))];
            hold on
            cl=plot(5*circle_x+center(1),(5*circle_y+center(2))+eye_y_offset,'-');
            set(cl,'Color','r','LineWidth',.5)
            cross=plot(center(1),center(2)+eye_y_offset,'MarkerSize',10,'Marker','+','MarkerEdgeColor','r','MarkerFaceColor','r');
            
        end
        
        
        axis equal
        set(gca,'xlim', [-35,35],'ylim', [5,35],'FontSize',GLO.fontsize-6)
        xlabel('Eye x position', 'fontsize', GLO.fontsize-4);
        ylabel('Eye y position', 'fontsize', GLO.fontsize-4);
        mtit(summary_2,  [plot_2_title, ' ', batch_title, ' ', fn{q}, ' batch ', num2str(k), ' plotting errors? ' num2str(errors) ], 'xoff', -0.0, 'yoff', 0.04, 'color', [0 0 0], 'fontsize', GLO.fontsize,'Interpreter', 'none');
        
        if GLO.create_pdf
            export_fig([GLO.folder_to_save batch_title, ' ', fn{q}, ' ' plot_2_title ' batch ' num2str(k), ' plotting errors ' num2str(errors) ], '-pdf','-transparent') % pdf by run
            close(gcf)
        end
    end
    
end
a=1;
