function beh_compare_groups(batch,testing)
global GLO
close all




Plot_settings.saccades.effectors    = GLO.saccades_effectors ;
Plot_settings.reaches.effectors     = GLO.reaches_effectors ;
Plot_settings.types                 = GLO.types_to_plot ;

Plot_settings.saccades.n_rows       = numel(Plot_settings.types);
Plot_settings.saccades.n_columns    = numel(Plot_settings.saccades.effectors);
Plot_settings.reaches.n_rows        = numel(Plot_settings.types);
Plot_settings.reaches.n_columns     = numel(Plot_settings.reaches.effectors);

Plot_settings.saccades.effectors_raw_xy     = GLO.saccades.effectors_raw_xy;
Plot_settings.reaches.effectors_raw_xy      = GLO.reaches.effectors_raw_xy;
Plot_settings.saccades.n_columns_raw_xy     = numel(Plot_settings.saccades.effectors_raw_xy);
Plot_settings.reaches.n_columns_raw_xy      = numel(Plot_settings.reaches.effectors_raw_xy);

Plot_settings.colors.IN             = [.5 .5 .5];
Plot_settings.colors.CH             = [0 0 0];
Plot_settings.colors.RH             = [0 1 0;0 0.5 0];
Plot_settings.colors.LH             = [0 1 1;0 0 1];
Plot_settings.colors.sac            = [1 0 0;0.5 0 0];
Plot_settings.multiplier            = [1, 0.5];

Plot_settings.markers.B             = 'o';
Plot_settings.markers.B_R           = '>';
Plot_settings.markers.B_L           = '<';

Plot_settings.markersize            = 10;
Plot_settings.linewidth             = 2;
Plot_settings.linestyle.CH          = '-';
Plot_settings.linestyle.IN          = ':';


fn_g = fieldnames(batch.out_stru_ext);

if numel(batch.unique_pos) == 1
    batch.unique_pos.Experimental = batch.unique_pos.Control;
end

for i = 1:numel(fn_g)
    Group_temp(1,i)=batch.out_stru_ext.(fn_g{i});
    monkey{1,i} = fn_g{i};
    Positions(1,i)=batch.unique_pos.(fn_g{i});
    
    %in order to plot target on raw eye traces when dissociated reaches
    if isempty(Positions(1,i).saccades) & strcmp(Plot_settings.saccades.effectors,'4')
        Positions(1,i).saccades = Positions(1,i).reaches;
        Positions(1,i).saccades_tar_rad = Positions(1,i).reaches_tar_rad;
        Positions(1,i).saccades_tar_siz = Positions(1,i).reaches_tar_siz;
    end
    
end



GLO.subject = monkey;

Group = Group_temp;
unique_saccade_positions    =   unique([Positions.saccades]);
unique_reach_positions      =   unique([Positions.reaches]);

for i = 1:numel(Group)
    Positions(1,i).saccades =   unique_saccade_positions;
    Positions(1,i).reaches  =   unique_reach_positions;
end



% print_out = [GLO.subject{1} '(a)' ' Vs ' GLO.subject{2} '(b)'];
print_out = [''];

reaches_saccades =fieldnames(Group);
reaches_saccades = {'saccades'}; % to avoid empty plots when only 1 effector

for rs=1:numel(reaches_saccades)
    sac_rea=reaches_saccades{rs};
    if strcmp(sac_rea,'reaches')
        %         parameters={'lat','dur','endpoints_per_position','endpoints_per_position_s','endpoints_per_position_a','accuracy_xy','precision_xy','successful','lat_residuals_sac_rea','lat_raw_sac_rea','lat_r_residuals','lat_slo_residuals','lat_int_residuals','lat_difference_sac_rea','ini_fix','ini_abort','abort_raw_x', 'success_raw_x'};
        parameters={'ini_abort','lat','dur','ini_fix','dur_fix','endpoints_per_position','endpoints_per_position_s','endpoints_per_position_a','endpoints_per_position_t_a','accuracy_xy','successful','side_selection','abort_raw_x','success_raw_x'};
        
        eye_or_hand_evaluated = ' Hand';
    else
                  parameters={'lat','dur','endpoints_per_position','endpoints_per_position_s','endpoints_per_position_a','accuracy_xy','successful','velocity','abort_raw_x','success_raw_x','side_selection'};
%          parameters={'side_selection'};
        eye_or_hand_evaluated = ' Eye';
    end
    for s_p=1:numel(parameters);
        par=parameters{s_p};
        
        switch par
            case 'lat';                         par_title = 'Reaction time to target';
            case 'dur';                         par_title = 'Movement time to target';
            case 'endpoints_per_position';      par_title = 'accuracy all';
            case 'endpoints_per_position_s';    par_title = 'accuracy successful only';
            case 'endpoints_per_position_a';    par_title = 'accuracy aborted only';
            case 'endpoints_per_position_t_a';  par_title = 'accuracy hand target aborted only';
            case 'velocity';                    par_title = 'velocity';
            case 'lat_r';                       par_title = 'latency corr';
            case 'lat_r_residuals';             par_title = 'latency corr, residuals';
            case 'lat_slo';                     par_title = 'latency corr slope';
            case 'lat_slo_residuals';           par_title = 'latency corr slope, residuals';
            case 'lat_int';                     par_title = 'latency corr intercept';
            case 'lat_int_residuals';           par_title = 'latency corr intercept, residuals';
            case 'dur_r';                       par_title = 'duration corr';
            case 'dur_r_residuals';             par_title = 'duration corr, residuals';
            case 'accuracy_xy_r';               par_title = 'accuracy corr';
            case 'accuracy_xy_r_residuals';     par_title = 'accuracy corr, residuals';
            case 'successful';                  par_title = 'success rate';
            case 'lat_raw_sac_rea';             par_title = 'correlation raw saccade reaches';
            case 'lat_residuals_sac_rea';       par_title = 'correlation residuals saccade reaches';
            case 'lat_difference_sac_rea';      par_title = 'latency difference reach minus saccade';
            case 'abort_raw_x';                 par_title = 'Error raw traces';
            case 'success_raw_x';               par_title = 'Success raw traces';
            case 'ini_fix';                     par_title = 'Reaction time to fixation';
            case 'dur_fix';                     par_title = 'Movement time to fixation';
            case 'ini_abort';                   par_title = 'Abort by using incorrect hand';
            case 'accuracy_xy';                 par_title = 'euclidean distance to target';
            case 'precision_xy';                par_title = 'euclidean cloud spread';
            case 'side_selection';              par_title = 'side selection';
        end
        
        %         print_out2  = [print_out ',' eye_or_hand_evaluated];
        print_out2  = [eye_or_hand_evaluated];
        
        %         idx_p   = findstr(par, 'residuals');
        %         if ~isempty(idx_p)
        %             par_sig = [par(1:idx_p-3) 'p' par(idx_p-1:end)];
        %         else
        par_sig = 'lat_p_residuals';
        %         end
        precision = 0;
        if strcmp(par,'precision_xy'), precision = 1; par='accuracy_xy'; end
        
        if (~strcmp(par,'ini_fix')&& ~strcmp(par,'dur_fix') && ~strcmp(par,'lat_residuals_sac_rea') && ~strcmp(par,'endpoints_per_position') && ~strcmp(par,'endpoints_per_position_s') && ~strcmp(par,'endpoints_per_position_a')&& ~strcmp(par,'endpoints_per_position_t_a')  && ~strcmp(par,'ini_abort') && ~strcmp(par,'lat_raw_sac_rea') ...
                && ~strcmp(par,'abort_raw_states') && ~strcmp(par,'abort_raw_time_axis') && ~strcmp(par,'abort_raw_x') && ~strcmp(par,'abort_raw_y') ...
                && ~strcmp(par,'success_raw_states') && ~strcmp(par,'success_raw_time_axis') && ~strcmp(par,'success_raw_x') && ~strcmp(par,'success_raw_y')) ...
                && (any(ismember(1,GLO.summary)) || any(ismember(2,GLO.summary)) || any(ismember(10,GLO.summary)) || any(ismember(11,GLO.summary)) || any(ismember(-1,GLO.summary)))
            
            if (any(ismember(1,GLO.summary)) || any(ismember(-1,GLO.summary)))
                % ERROR BAR FIGURES
                plot_title                                                = [' Summary 1, ' print_out2 ', ' par_title ];
                summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
                subplot_indexes=[];
                clear ylim_sp;
                for t=1:numel(Plot_settings.types)
                    type=Plot_settings.types{t};
                    if isfield(Group(1).(sac_rea),par) && (any(ismember(1,GLO.summary)) || any(ismember(-1,GLO.summary)))
                        subplot(Plot_settings.(sac_rea).n_rows,1,t);
                        if ~isfield(Group(1).(sac_rea),par_sig), Group(1).(sac_rea).(par_sig)=NaN; end
                        if ~isfield(Group(2).(sac_rea),par_sig), Group(2).(sac_rea).(par_sig)=NaN; end
                        isdata= temp_means_bars(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,type,Plot_settings.(sac_rea).effectors,Positions(1).(sac_rea),Plot_settings,par,batch.stat,Group(1).(sac_rea).(par_sig),Group(2).(sac_rea).(par_sig),precision);
                        if isdata
                            subplot_indexes=[subplot_indexes t];
                        end
                    end
                end
                if isfield(Group(1).(sac_rea),par)
                    for sp=1:numel(subplot_indexes)
                        subplot(Plot_settings.(sac_rea).n_rows,1,subplot_indexes(sp));
                        ylim_sp(sp,:)=get(gca,'ylim');
                    end
                    for sp=1:numel(subplot_indexes)
                        subplot(Plot_settings.(sac_rea).n_rows,1,subplot_indexes(sp));
                        set(gca,'ylim',[min(ylim_sp(:,1)) max(ylim_sp(:,2))]);
                    end
                end
                if strcmp(par, 'successful') || strcmp(par, 'side_selection')
                    set(gca,'ylim',[0 1]);
                end
                title_and_save(summary_figure,plot_title);
                
                if strcmp(par,'endpoints_per_position') || strcmp(par,'endpoints_per_position_s') || strcmp(par,'endpoints_per_position_a')|| strcmp(par,'endpoints_per_position_t_a') || strcmp(par,'successful')
                    continue
                end
            end
            
            
            
            if (any(ismember(11,GLO.summary)) || any(ismember(-1,GLO.summary))) %&& (strcmp(par,'lat') || strcmp(par,'successful'))
                % Session by session dots plot
                plot_title                                                = [' Summary 11, ' print_out2 ', ' par_title ];
                summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
                subplot_indexes=[];
                clear ylim_sp;
                for t=1:numel(Plot_settings.types)
                    type=Plot_settings.types{t};
                    if isfield(Group(1).(sac_rea),par) && (any(ismember(11,GLO.summary)) || any(ismember(-1,GLO.summary)))
                        subplot(Plot_settings.(sac_rea).n_rows,1,t);
                        isdata= temp_means_bars_consecutive(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,type,Plot_settings.(sac_rea).effectors,Positions(1).(sac_rea),Plot_settings,par,batch.stat,precision);
                        if isdata
                            subplot_indexes=[subplot_indexes t];
                        end
                    end
                end
                if isfield(Group(1).(sac_rea),par)
                    for sp=1:numel(subplot_indexes)
                        subplot(Plot_settings.(sac_rea).n_rows,1,subplot_indexes(sp));
                        ylim_sp(sp,:)=get(gca,'ylim');
                    end
                    for sp=1:numel(subplot_indexes)
                        subplot(Plot_settings.(sac_rea).n_rows,1,subplot_indexes(sp));
                        set(gca,'ylim',[min(ylim_sp(:,1)) max(ylim_sp(:,2))]);
                    end
                end
                title_and_save(summary_figure,plot_title);
                
                if strcmp(par,'endpoints_per_position') || strcmp(par,'endpoints_per_position_s') || strcmp(par,'endpoints_per_position_a') || strcmp(par,'endpoints_per_position_t_a') ||strcmp(par,'successful')
                    continue
                end
            end
            
            
            
            
            if (any(ismember(10,GLO.summary)) || any(ismember(-1,GLO.summary)))
                % CHOICE - INSTRUCTED FIGURES
                plot_title                                                = [' Summary 10, ' print_out2 ', ' par_title ' CH - IN '];
                summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
                subplot_indexes=[];
                clear ylim_sp;
                for t=1:numel(Plot_settings.types)
                    type=Plot_settings.types{t};
                    if isfield(Group(1).(sac_rea),par) && (any(ismember(10,GLO.summary)) || any(ismember(-1,GLO.summary)))
                        subplot(Plot_settings.(sac_rea).n_rows,1,t);
                        isdata= temp_means_bars_ch_in(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,type,Plot_settings.(sac_rea).effectors,Positions(1).(sac_rea),Plot_settings,par,batch.stat,Group(1).(sac_rea).(par_sig),Group(2).(sac_rea).(par_sig));
                        if isdata
                            subplot_indexes=[subplot_indexes t];
                        end
                    end
                end
                if isfield(Group(1).(sac_rea),par)
                    for sp=1:numel(subplot_indexes)
                        subplot(Plot_settings.(sac_rea).n_rows,1,subplot_indexes(sp));
                        ylim_sp(sp,:)=get(gca,'ylim');
                    end
                    for sp=1:numel(subplot_indexes)
                        subplot(Plot_settings.(sac_rea).n_rows,1,subplot_indexes(sp));
                        set(gca,'ylim',[min(ylim_sp(:,1)) max(ylim_sp(:,2))]);
                    end
                end
                title_and_save(summary_figure,plot_title);
                
                if strcmp(par,'endpoints_per_position') || strcmp(par,'endpoints_per_position_s') || strcmp(par,'endpoints_per_position_a') || strcmp(par,'endpoints_per_position_t_a') || strcmp(par,'successful')
                    continue
                end
            end
            
            if any(ismember(2,GLO.summary)) || any(ismember(-1,GLO.summary))
                % HISTOGRAM FIGURES
                plot_title                                                = [' Summary 2, ' print_out2 ', ' par_title ];
                summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
                subplot_indexes_h=[];
                clear ylim_sp;
                for t=1:numel(Plot_settings.types)
                    for e=1:numel(Plot_settings.(sac_rea).effectors)
                        type=Plot_settings.types{t};
                        effector=Plot_settings.(sac_rea).effectors{e};
                        hi((t-1)*2*Plot_settings.(sac_rea).n_columns + e*2-1)=subplot(Plot_settings.(sac_rea).n_rows, Plot_settings.(sac_rea).n_columns*2,(t-1)*2*Plot_settings.(sac_rea).n_columns + e*2-1);
                        if ~isfield(Group(1).(sac_rea),par), Group(1).(sac_rea).(par)=NaN; end
                        if ~isfield(Group(2).(sac_rea),par), Group(2).(sac_rea).(par)=NaN; end
                        isdata= plot_par_histograms(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,type,Plot_settings.(sac_rea).effectors,Positions(1).(sac_rea),Plot_settings,par,precision,e,1);
                        hi((t-1)*2*Plot_settings.(sac_rea).n_columns + e*2)=subplot(Plot_settings.(sac_rea).n_rows, Plot_settings.(sac_rea).n_columns*2,(t-1)*2*Plot_settings.(sac_rea).n_columns + e*2);
                        isdata= plot_par_histograms(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,type,Plot_settings.(sac_rea).effectors,Positions(1).(sac_rea),Plot_settings,par,precision,e,2);
                        if isdata
                            subplot_indexes_h=[subplot_indexes_h t];
                        end
                    end
                end
                
                title_and_save(summary_figure,plot_title);
            end
            
        elseif (strcmp(par,'endpoints_per_position') || strcmp(par,'endpoints_per_position_s') || strcmp(par,'endpoints_per_position_a') || strcmp(par,'endpoints_per_position_t_a')) && (any(ismember(3,GLO.summary)) || any(ismember(-1,GLO.summary))) % && strcmp(effector,GLO.type_of_free_gaze)
            % ACCURACY BAR FIGURES
            plot_title                                                = [' Summary 3, ' print_out2 ', ' par_title ];
            summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            subplot_indexes=[];
            clear ylim_sp;
            for t=1:numel(Plot_settings.types)
                type=Plot_settings.types{t};
                if (strcmp(par,'endpoints_per_position') || strcmp(par,'endpoints_per_position_s') || strcmp(par,'endpoints_per_position_a') || strcmp(par,'endpoints_per_position_t_a')) && isfield(Group(1).(sac_rea),par) && isfield(Group(2).(sac_rea),par) && (any(ismember(3,GLO.summary)) || any(ismember(-1,GLO.summary)))
                    
                    for e=1:numel(Plot_settings.(sac_rea).effectors)
                        effector=Plot_settings.(sac_rea).effectors{e};
                        subplot(Plot_settings.(sac_rea).n_rows,Plot_settings.(sac_rea).n_columns,(t-1)*Plot_settings.(sac_rea).n_columns + e);
                        %                         plot_accuracy_internal(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,type,effector,Positions(1).(sac_rea),Plot_settings,par);
                        plot_accuracy_internal_ellipse(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,type,effector,...
                            Positions(1).(sac_rea),Positions(1).([sac_rea '_tar_rad']), Positions(1).([sac_rea '_tar_siz']),...
                            Plot_settings,par,batch.stat.groups.(sac_rea).(par),testing,precision);
                        axis('equal')
                        set(gca,'ylim',[-30 30],'Xlim',[-30 30])
                        
                        %  axis equal
                        %   set(gca,'xlim',[-30 -18],'ylim', [-4 4]) %4
                        % set(gca,'xlim',[18 30],'ylim', [-4 4]) %3
                        %   set(gca,'xlim',[-18 -6],'ylim', [-4 4]) %2
                        %   set(gca,'xlim',[6 18],'ylim', [-4 4]) %1
                    end
                end
            end
            %             if (strcmp(par,'endpoints_per_position') || strcmp(par,'endpoints_per_position_s') || strcmp(par,'endpoints_per_position_a'))   && isfield(Group(1).(sac_rea),par) && isfield(Group(2).(sac_rea),par)
            %                 for sp=1:numel(subplot_indexes)
            %                     subplot(Plot_settings.(sac_rea).n_rows,Plot_settings.(sac_rea).n_columns,subplot_indexes(sp));
            %                     %                     ylim_sp(sp,:)=get(gca,'ylim');
            %
            %                 end
            %                 for sp=1:numel(subplot_indexes)
            %                     subplot(Plot_settings.(sac_rea).n_rows,Plot_settings.(sac_rea).n_columns,subplot_indexes(sp));
            %                     %                      set(gca,'ylim',[min(ylim_sp(:,1)) max(ylim_sp(:,2))]);
            %
            %                 end
            %             end
            title_and_save(summary_figure,plot_title);
            
        elseif (strcmp(par,'lat_residuals_sac_rea') || strcmp(par,'lat_raw_sac_rea')) && (any(ismember(4,GLO.summary)) || any(ismember(-1,GLO.summary))) % && strcmp(effector,GLO.type_of_free_gaze)
            % CORRELATION FIGURES
            plot_title                                                = [' Summary 4, ' print_out2 ', ' par_title ];
            summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            subplot_indexes_h=[];
            clear ylim_sp;
            for t=1:numel(Plot_settings.types)
                
                for e=1:numel(Plot_settings.(sac_rea).effectors)
                    type=Plot_settings.types{t};
                    effector=Plot_settings.(sac_rea).effectors{e};
                    co((t-1)*2*Plot_settings.(sac_rea).n_columns + e*2-1)=subplot(Plot_settings.(sac_rea).n_rows, Plot_settings.(sac_rea).n_columns*2,(t-1)*2*Plot_settings.(sac_rea).n_columns + e*2-1);
                    if ~isfield(Group(1).(sac_rea),par), Group(1).(sac_rea).(par)=NaN; end
                    if ~isfield(Group(2).(sac_rea),par), Group(2).(sac_rea).(par)=NaN; end
                    plot_par_correlations(Group(1).(sac_rea),Group(2).(sac_rea),sac_rea,type,Plot_settings.(sac_rea).effectors,Positions(1).(sac_rea),Plot_settings,par,e,1);
                    co((t-1)*2*Plot_settings.(sac_rea).n_columns + e*2)=subplot(Plot_settings.(sac_rea).n_rows, Plot_settings.(sac_rea).n_columns*2,(t-1)*2*Plot_settings.(sac_rea).n_columns + e*2);
                    plot_par_correlations(Group(1).(sac_rea),Group(2).(sac_rea),sac_rea,type,Plot_settings.(sac_rea).effectors,Positions(1).(sac_rea),Plot_settings,par,e,2);
                    %                     if isdata
                    %                         subplot_indexes_h=[subplot_indexes_h t];
                    %                     end
                end
            end
            title_and_save(summary_figure,plot_title);
            
        elseif (strcmp(par,'ini_fix') && strcmp('reaches',sac_rea)) && (any(ismember(5,GLO.summary)) || any(ismember(-1,GLO.summary)))
            
            % HAND RELEASE FIGURE
            plot_title                                                = [' Summary 5, ' print_out2 ', ' par_title ];
            summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            hand_release(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,Plot_settings,par);
            title_and_save(summary_figure,plot_title);
            
        elseif (strcmp(par,'dur_fix') && strcmp('reaches',sac_rea)) && (any(ismember(5,GLO.summary)) || any(ismember(-1,GLO.summary)))
            % Fixation Movement time
            plot_title                                                = [' Summary 5, ' print_out2 ', ' par_title ];
            summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            hand_release(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,Plot_settings,par);
            title_and_save(summary_figure,plot_title);
            
            
            
        elseif (strcmp(par,'abort_raw_x') || strcmp(par,'success_raw_x') )&& (any(ismember(6,GLO.summary)) || any(ismember(7,GLO.summary)) || any(ismember(-1,GLO.summary)))
            if (any(ismember(6,GLO.summary)) || any(ismember(-1,GLO.summary))) && GLO.keep_raw_output
                % RAW TRACES FIGURES
                plot_title                                                = [' Summary 6, ' print_out2 ', ' par_title ];
                summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
                clear ylim_sp;
                for t=1:numel(Plot_settings.types)
                    
                    for e=1:numel(Plot_settings.(sac_rea).effectors_raw_xy)
                        type=Plot_settings.types{t};
                        effector=Plot_settings.(sac_rea).effectors_raw_xy{e};
                        ra((t-1)*Plot_settings.(sac_rea).n_columns_raw_xy + e)=subplot(Plot_settings.(sac_rea).n_rows, Plot_settings.(sac_rea).n_columns_raw_xy,  (t-1)*Plot_settings.(sac_rea).n_columns_raw_xy + e);
                        raw_plotting(Group(1).(sac_rea),Group(2).(sac_rea),sac_rea,type,Plot_settings.(sac_rea).effectors_raw_xy,Positions(1).(sac_rea),Plot_settings,par,e,1,...
                            Positions(1).([sac_rea '_tar_rad']), Positions(1).([sac_rea '_tar_siz']),Positions(1).([sac_rea '_fix_siz']),Positions(1).([sac_rea '_fix_rad']));
                        axis('equal')
                        set(gca,'ylim',[-30 30],'Xlim',[-30 30])
                        
                    end
                end
                title_and_save(summary_figure,plot_title);
            end
            
            %             if (any(ismember(6.5,GLO.summary)) || any(ismember(-1,GLO.summary))) && GLO.keep_raw_output
            %                 % RAW TRACES FIGURES (not what it does but not relevant if
            %                 % no saccades information
            %                 plot_title                                                = [' Summary 6.5, ' print_out2 ', ' par_title ];
            %                 summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            %                 clear ylim_sp;
            %                 for t=1:numel(Plot_settings.types)
            %
            %                     for e=1:numel(Plot_settings.(sac_rea).effectors_raw_xy)
            %                         type=Plot_settings.types{t};
            %                         effector=Plot_settings.(sac_rea).effectors_raw_xy{e};
            %                         ra((t-1)*Plot_settings.(sac_rea).n_columns_raw_xy + e)=subplot(Plot_settings.(sac_rea).n_rows, Plot_settings.(sac_rea).n_columns_raw_xy,  (t-1)*Plot_settings.(sac_rea).n_columns_raw_xy + e);
            %                         raw_EH_first_reach_or_saccade(Group(1),Group(2),sac_rea,type,Plot_settings.(sac_rea).effectors_raw_xy,Positions(1),Plot_settings,par,e,1);
            %                         axis('equal')
            %                         set(gca,'ylim',[-30 30],'Xlim',[-30 30])
            %                     end
            %                 end
            %                 title_and_save(summary_figure,plot_title);
            %             end
            %
            %             if (any(ismember(7,GLO.summary)) || any(ismember(-1,GLO.summary))) && GLO.keep_raw_output
            %                 plot_title                                                = [' Summary 7x, ' print_out2 ', ' par_title ];
            %                 current_axis = 'x';
            %                 summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            %                 clear ylim_sp;
            %                 for t=1:numel(Plot_settings.types)
            %
            %                     for e=1:numel(Plot_settings.(sac_rea).effectors_raw_xy)
            %                         type=Plot_settings.types{t};
            %                         effector=Plot_settings.(sac_rea).effectors_raw_xy{e};
            %                         ra((t-1)*Plot_settings.(sac_rea).n_columns_raw_xy + e)=subplot(Plot_settings.(sac_rea).n_rows, Plot_settings.(sac_rea).n_columns_raw_xy,  (t-1)*Plot_settings.(sac_rea).n_columns_raw_xy + e);
            %                         raw_plotting_xy(Group(1).(sac_rea),Group(2).(sac_rea),sac_rea,type,Plot_settings.(sac_rea).effectors_raw_xy,Positions(1).(sac_rea),Plot_settings,par,e,1,current_axis,...
            %                             Positions(1).([sac_rea '_tar_rad']), Positions(1).([sac_rea '_tar_siz']));
            %
            %                         if GLO.testing_patient
            %                             set(gca,'ylim',[-30 30],'Xlim',[-0.1 2])
            %                         else
            %                             set(gca,'ylim',[-30 30],'Xlim',[-0.1 2])
            %                         end
            %                         %                     axis('equal')
            %                     end
            %                 end
            %                 title_and_save(summary_figure,plot_title);
            %
            %
            %                 plot_title                                                = [' Summary 7y, ' print_out2 ', ' par_title ];
            %                 current_axis = 'y';
            %                 summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            %                 clear ylim_sp;
            %                 for t=1:numel(Plot_settings.types)
            %
            %                     for e=1:numel(Plot_settings.(sac_rea).effectors_raw_xy)
            %                         type=Plot_settings.types{t};
            %                         effector=Plot_settings.(sac_rea).effectors_raw_xy{e};
            %                         ra((t-1)*Plot_settings.(sac_rea).n_columns_raw_xy + e)=subplot(Plot_settings.(sac_rea).n_rows, Plot_settings.(sac_rea).n_columns_raw_xy,  (t-1)*Plot_settings.(sac_rea).n_columns_raw_xy + e);
            %                         raw_plotting_xy(Group(1).(sac_rea),Group(2).(sac_rea),sac_rea,type,Plot_settings.(sac_rea).effectors_raw_xy,Positions(1).(sac_rea),Plot_settings,par,e,1,current_axis,...
            %                             Positions(1).([sac_rea '_tar_rad']), Positions(1).([sac_rea '_tar_siz']));
            %                         %                         axis('equal')
            %                         if GLO.testing_patient
            %                             set(gca,'ylim',[-30 30],'Xlim',[-0.1 2])
            %                         else
            %                             set(gca,'ylim',[-30 30],'Xlim',[-0.1 2])
            %                         end
            %                         %
            %                     end
            %                 end
            %                 title_and_save(summary_figure,plot_title);
            %             end
            
        elseif (strcmp(par,'ini_abort') && strcmp('reaches',sac_rea)) && (any(ismember(12,GLO.summary)) || any(ismember(-1,GLO.summary)))
            
            % HAND RELEASE FIGURE
            plot_title                                                = [' Summary 12, ' print_out2 ', ' par_title ];
            summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
            hand_release_abort(Group(1).(sac_rea).(par),Group(2).(sac_rea).(par),sac_rea,Plot_settings,par);
            title_and_save(summary_figure,plot_title);
            
            
        end
    end
    
    if (any(ismember(9,GLO.summary)) || any(ismember(-1,GLO.summary)))
        % ERROR PLOT
        [G Abort_codes] = error_str_2_num(Group(1).(sac_rea).abort_code, Group(2).(sac_rea).abort_code);
        
        plot_title                                                = [' Summary 9, '  'Error summary ' sac_rea];
        summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
        subplot_indexes=[];
        clear ylim_sp;
        for t=1:numel(Plot_settings.types)
            type=Plot_settings.types{t};
            subplot(Plot_settings.(sac_rea).n_rows,1,t);
            errors(G, Abort_codes,sac_rea,type,Plot_settings.(sac_rea).effectors,Positions(1).(sac_rea),Plot_settings,par,batch.stat,Group(1).(sac_rea).(par_sig),Group(2).(sac_rea).(par_sig));
        end
        %         if isfield(Group(1).(sac_rea),par)
        %             for sp=1:numel(subplot_indexes)
        %                 subplot(Plot_settings.(sac_rea).n_rows,1,subplot_indexes(sp));
        %                 ylim_sp(sp,:)=get(gca,'ylim');
        %             end
        %             for sp=1:numel(subplot_indexes)
        %                 subplot(Plot_settings.(sac_rea).n_rows,1,subplot_indexes(sp));
        %                 set(gca,'ylim',[min(ylim_sp(:,1)) max(ylim_sp(:,2))]);
        %             end
        %         end
        set(gca,'ylim',[-0.3 1]) ;
        title_and_save(summary_figure,plot_title);
        
    end
    
    
    
end

%% FILELIST PLOT
if (any(ismember(8,GLO.summary)) || any(ismember(-1,GLO.summary)))
    plot_title                                                = [' Summary 8, ' ' filenames '];
    summary_figure                                            = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    plot_file_names(batch,summary_figure)
    title_and_save(summary_figure,plot_title);
end



end

%% EXTERNAL LOOP PLOTTING FUNCTIONS
function isdata=temp_means_bars(input_group1,input_group2,sac_rea,type,effectors,unique_pos,Plot_settings,par,stat,sig_residuals_group1,sig_residuals_group2,precision)
global GLO

switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='Delayed';
end
hands={'LH','RH'};
sides={'L','R'};
types = {'IN','CH'};
side_labels={'LEFT','RIGHT'};
countersides={'R','L'};
hold on

groups=[input_group1; input_group2];
sig_residuals=[sig_residuals_group1; sig_residuals_group2];
idx=-2;
isdata=false;

for s=1:numel(sides)
    side=sides{s};
    side_label=side_labels{s};
    counter_side=countersides{s};
    if strcmp(sac_rea,'saccades')
        if strcmp(side,'R')
            idx = 4;
        end
    end
    idx=idx+1;
    slh(s)=text(idx+2,1,[side_label ' HEMISPACE']);
    
    for e=1:numel(effectors)
        idx=idx+1;
        effector=effectors{e};
        switch str2double(effector)
            case 0, effector_label='saccades'; case 1, effector_label='free gaze reaches'; case 2, effector_label='joint eye-hand'; case 3, effector_label='disociated saccades'; case 4, effector_label='dissociated reaches'; case 6, effector_label='free gaze with fixation';
        end
        %         elh(numel(effectors)*(s-1)+e)=text(idx+1,1,effector_label);
        type_effector=['t_' type '_e_' effector];
        
        if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
            for h=1:numel(hands)
                hand=['_' hands{h}];idx_bar_start=idx+1;
                [idx e_bar_tmp]=plot_internal_per_condition(idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),par,stat,sig_residuals,precision);
                
                [idx_stat_IN{idx} idx_stat_CH{idx}]=stat_sig(idx,sac_rea,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),par,stat,precision);
                e_bar(idx_bar_start:idx)=e_bar_tmp(idx_bar_start:idx);
            end
        elseif strcmp(sac_rea,'saccades') && str2double(effector)==0
            if strcmp(side,'R')
                idx = 6;
            end
            
            hand=[];idx_bar_start=idx+1;
            [idx e_bar_tmp]=plot_internal_per_condition(idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.sac,par,stat,sig_residuals,precision);
            
            [idx_stat_IN{idx} idx_stat_CH{idx}]=stat_sig(idx,sac_rea,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.sac,par,stat,precision);
            e_bar(idx_bar_start:idx)=e_bar_tmp(idx_bar_start:idx);
        end
        title(type_labels)
    end
end


% if strcmp(sac_rea,'reaches')
%      comp_IN={[1 2] [3 4], [6 7] [8 9], [12 13] [14 15], [17 18] [19 20]};
%      comp_LS={[1 2 3 4] [6 7 8 9]};
%      comp_RS={[[1 2 3 4]+11] [[6 7 8 9]+11]};
%      comp_LH={[1 3 6 8 12 14 17 19]};
%      comp_RH={[1 3 6 8 12 14 17 19]+1};

if  strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
    comp_IN={[1 2] [3 4], [7 8] [9 10]};
    comp_LS={[1 2 3 4] [7 8 9 10]};
    comp_RS={[1 2 3 4] [7 8 9 10]};
    comp_LH={[1 3 7 9]};
    comp_RH={[1 3 7 9]+1};
elseif (strcmp(sac_rea,'saccades') && str2double(effector)==0)
    comp_IN={[1 2]  [7 8]};
    comp_LS={[1 2] [7 8]};
    comp_RS={[1 2] [7 8]};
    comp_LH={[]};
    comp_RH={[]};
    %     comp_IN={[1 2]  [7 8]}; new try
    %     comp_LS={[1 7]};
    %     comp_RS={[2 8]};
    %     comp_LH={[]};
    %     comp_RH={[]};
end


% else
% %     comp_IN={[1 2], [4 5], [6 7], [9 10] [11 12], [15 16], [18 19] [20 21], [23 24] [25 26]};
% %     comp_LS={[1 2] [4 5 6 7] [9 10 11 12]};
% %     comp_RS={[[1 2]+14 [4 5 6 7]+14 [9 10 11 12]+14]};
% %     comp_LH={[1 4 6 9 11 15 18 20 23 25]};
% %     comp_RH={[1 4 6 9 11 15 18 20 23 25]+1};
%
%  comp_IN={[1 2] [3 4], [7 8] [9 10]};
%     comp_LS={[1 2 3 4] [7 8 9 10]};
%     comp_RS={[1 2 3 4] [7 8 9 10]};
%     comp_LH={[1 3 7 9]};
%     comp_RH={[1 3 7 9]+1};
%
% end
comp_IN_bar=comp_IN;


if GLO.calculate_statististics ==1 && GLO.plot_statististics == 1
    if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
        %         comp_IN={[1 2] [3 4], [6 7] [8 9], [12 13] [14 15], [17 18] [19 20]};
        comp_CH=comp_IN;
        ds      = 6;
        dh      = 2;
        dc      = 1;
        
        %         compIN.space={[1 1+ds], [2 2+ds], [3 3+ds], [4 4+ds], [6 6+ds], [7 7+ds], [8 8+ds], [9 9+ds]};
        %         compIN.hand={[1 1+dh], [2 2+dh], [12 12+dh] [13 13+dh], [6 6+dh] [7 7+dh], [17 17+dh] [18 18+dh]};
        %         compIN.choice={[1 1+ds], [2 2+ds], [3 3+ds], [4 4+ds], [12 5+ds], [13 6+ds], [14 7+ds], [15 8+ds]}; % placeholder
        
        compIN.space={[1 1+ds], [2 2+ds], [3 3+ds], [4 4+ds]};
        compIN.hand={[1 1+dh], [2 2+dh],  [7 7+dh],[8 8+dh]};
        compIN.choice={[1 1+ds], [2 2+ds], [3 3+ds], [4 4+ds]}; % placeholder
        
        
    else
        comp_CH=comp_IN;
        ds      = 6;
        dh      = 2;
        dc      = 0;
        
        compIN.space={[1 1+ds], [2 2+ds]};
        compIN.hand={};
        compIN.choice={[1 1+ds], [2 2+ds]}; % placeholder
        
        %         %         comp_IN={[1 2], [4 5], [6 7], [9 10] [11 12], [15 16], [18 19] [20 21], [23 24] [25 26]};
        %         comp_CH=comp_IN;
        %        % ds      = 14;
        %         ds      = 6;
        %         dh      = 2;
        %         dc      = 0;
        % %         compIN.space={[1 1+ds], [2 2+ds], [NaN NaN], [NaN NaN], [4 4+ds], [5 5+ds], [6 6+ds], [7 7+ds], [9 9+ds], [10 10+ds], [11 11+ds], [12 12+ds]};
        % %         compIN.hand={[NaN NaN], [NaN NaN], [NaN NaN], [NaN NaN],[4 4+dh], [5 5+dh], [18 18+dh] [19 19+dh], [9 9+dh] [10 10+dh],[23 23+dh] [24 24+dh]};
        % %         compIN.choice={[1 1+ds], [2 2+ds], [4 4+ds], [5 5+ds], [6 6+ds], [7 7+ds], [9 9+ds], [10 10+ds], [11 11+ds], [12 12+ds]};     % placeholder
        %
        %         compIN.space={[1 1+ds], [2 2+ds], [3 3+ds], [4 4+ds]};
        %         compIN.hand={[1 1+dh], [2 2+dh],  [7 7+dh],[8 8+dh]};
        %         compIN.choice={[1 1+ds], [2 2+ds], [3 3+ds], [4 4+ds]}; % placeholder
        %     end
        %     comp_IN_bar=comp_IN;
    end
    sig_IN=[]; sig_CH=[];
    for l=1:numel(comp_IN)
        sig_IN = [sig_IN idx_stat_IN{comp_IN{l}(2)}{2}];
        sig_CH = [sig_CH idx_stat_CH{comp_CH{l}(2)}{2}];
    end
    comp_IN(isnan(sig_IN))=[];
    sig_IN(isnan(sig_IN))=[];
    sig_IN(sig_IN>0.05)=NaN;
    comp_CH(isnan(sig_CH))=[];
    sig_CH(isnan(sig_CH))=[];
    sig_CH(sig_CH>0.05)=NaN;
end

%% setting limits and ticks
y_lim=get(gca,'ylim');
for s=1:numel(slh)
    pos=get(slh(s),'Position');
    set(slh(s),'Position',[pos(1) y_lim(2)+diff(y_lim)*2/10])
end
% for e=1:numel(elh)
%     pos=get(elh(e),'Position');
%     set(elh(e),'Position',[pos(1) y_lim(2)+diff(y_lim)*1/10])
% end
y_lim(2)=y_lim(2)+diff(y_lim)*3/10;

if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
    xticks=[1 2 3 4, 7 8 9 10];
    %     xticklabels={'LHa' 'LHb' 'RHa' 'RHb' 'LHa' 'LHb' 'RHa' 'RHb' 'LHa' 'LHb' 'RHa' 'RHb'};
    xticklabels={'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' };
else
    xticks=[1 2 7 8];
    xticklabels={'Con' 'Ina' 'Con' 'Ina' };
    %     xticks=[1 2, 4 5 6 7, 9 10 11 12, 15 16, 18 19 20 21, 23 24 25 26];
    %     %     xticklabels={'NHa' 'NHb' 'LHa' 'LHb' 'RHa' 'RHb' 'LHa' 'LHb' 'RHa' 'RHb' 'NHa' 'NHb' 'LHa' 'LHb' 'RHa' 'RHb'};
    %     xticklabels={'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina'};
end

set(gca,'ylim',y_lim,'xlim',[0,max(xticks)+1],'xtick',xticks,'XTickLabel',xticklabels);
ylabel(par);

if GLO.calculate_statististics ==1 && GLO.plot_statististics == 1
    if GLO.only_significant
        comp_IN=comp_IN(~isnan(sig_IN)&sig_IN<0.05);
        sig_IN=sig_IN(~isnan(sig_IN)&sig_IN<0.05);
        comp_CH=comp_CH(~isnan(sig_CH)&sig_CH<0.05);
        sig_CH=sig_CH(~isnan(sig_CH)&sig_CH<0.05);
    end
    s_IN = sigstar_eye_hand(comp_IN,sig_IN,0,Plot_settings.colors.IN);
    s_CH = sigstar_eye_hand(comp_CH,sig_CH,0,Plot_settings.colors.CH);
    
    if ~GLO.only_between_group_comparissons
        %comparissons={'space','hand','choice'};
        comparissons={'hand','space'};
        not_exisent_comparison_counter=0;
        for comps = 1:2
            for co=1:numel(comparissons)
                comparisson=comparissons{co};
                for e=1:numel(effectors)
                    for compc=1:2
                        switch comparisson
                            case 'space'
                                subcondition=['t_' type '_e_' effectors{e} '_L_' types{comps} '_' hands{compc}];
                            case 'hand'
                                subcondition=['t_' type '_e_' effectors{e} '_' sides{compc} '_' types{comps} '_LH'];
                        end
                        if strcmp(sac_rea,'saccades') && str2double(effectors{e})==0
                            subcondition(end-2:end)=[];
                        end
                        for g=1:2
                            if strcmp(sac_rea,'saccades') && str2double(effectors{e})==0 && (strcmp(comparisson,'hand') || compc==2)
                                continue;
                            end
                            
                            if ~isfield(stat.(comparisson)(g).(sac_rea).(par),(subcondition))
                                continue
                            else
                                if precision, totest_raw = 'raw_of_std'; totest_mean = 'raw_of_std'; else totest_raw = 'raw_of_raw'; totest_mean = 'raw_of_mean'; end
                                if strcmp(sac_rea,'saccades') && str2double(effectors{e})==0 %saccades
                                    if (GLO.trial_by_trial) ||  GLO.testing_patient %trial by trial
                                        if GLO.only_significant && (~isnan(stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_raw){4}) && stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_raw){4} < 0.05  )
                                            sigstar_eye_hand(compIN.(comparisson)(e*4+compc*2+g-6),stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_raw){4},0,Plot_settings.colors.(types{comps}));
                                        elseif ~GLO.only_significant
                                            sigstar_eye_hand(compIN.(comparisson)(e*4+compc*2+g-6),stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_raw){4},0,Plot_settings.colors.(types{comps}));
                                        end
                                    else %summaries per batch
                                        if GLO.only_significant && (~isnan(stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_mean){4}) && stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_mean){4} < 0.05 )
                                            sigstar_eye_hand(compIN.(comparisson)(e*4+compc*2+g-6),stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_mean){4},0,Plot_settings.colors.(types{comps}));
                                        elseif ~GLO.only_significant
                                            sigstar_eye_hand(compIN.(comparisson)(e*4+compc*2+g-6),stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_mean){4},0,Plot_settings.colors.(types{comps}));
                                        end
                                    end
                                else %reaches
                                    if (GLO.trial_by_trial) ||  GLO.testing_patient %trial by trial
                                        if GLO.only_significant && (~isnan(stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_raw){4}) && stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_raw){4} < 0.05  )
                                            sigstar_eye_hand(compIN.(comparisson)(e*4+compc*2+g-6-not_exisent_comparison_counter),stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_raw){4},0,Plot_settings.colors.(types{comps}));
                                        elseif ~GLO.only_significant
                                            sigstar_eye_hand(compIN.(comparisson)(e*4+compc*2+g-6-not_exisent_comparison_counter),stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_raw){4},0,Plot_settings.colors.(types{comps}));
                                        end
                                    else %summaries per batch
                                        if GLO.only_significant && (~isnan(stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_mean){4}) && stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_mean){4} < 0.05  )
                                            sigstar_eye_hand(compIN.(comparisson)(e*4+compc*2+g-6-not_exisent_comparison_counter),stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_mean){4},0,Plot_settings.colors.(types{comps}));
                                        elseif ~GLO.only_significant
                                            sigstar_eye_hand(compIN.(comparisson)(e*4+compc*2+g-6-not_exisent_comparison_counter),stat.(comparisson)(g).(sac_rea).(par).(subcondition).(totest_mean){4},0,Plot_settings.colors.(types{comps}));
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
hold on

% set(gca,'layer','top','ygrid','on','yminorgrid','on')

% comp_IN_mat=vertcat(comp_IN_bar{:});
% ho = patch([comp_IN_mat(:,1)-0.5 comp_IN_mat(:,1)+0.5 comp_IN_mat(:,1)+0.5 comp_IN_mat(:,1)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
% set(ho,'facecolor',[0.85 0.85 0.85])
% ha = patch([comp_IN_mat(:,2)-0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
% set(ha,'facecolor',[[0.65,0.65,0.65]])
%
% uistack(ho,'bottom');
% uistack(ha,'bottom');

% comp_S_mat=[horzcat(comp_LS{:})', horzcat(comp_RS{:})'];
% hLS = patch([comp_S_mat(:,1)-0.5 comp_S_mat(:,1)+0.5 comp_S_mat(:,1)+0.5 comp_S_mat(:,1)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
% set(hLS,'facecolor',[0.85 0.85 0.85])
% hRS = patch([comp_IN_mat(:,2)-0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
% set(hRS,'facecolor',[[0.65,0.65,0.65]])
%
%
% uistack(ho,'bottom');
% uistack(ha,'bottom');

end

function isdata=temp_means_bars_consecutive(input_group1,input_group2,sac_rea,type,effectors,unique_pos,Plot_settings,par,stat,precision)
global GLO

switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
hands={'LH','RH'};
sides={'L','R'};
side_labels={'LEFT','RIGHT'};
countersides={'R','L'};
hold on

groups=[input_group1; input_group2];
%  sig_residuals=[sig_residuals_group1; sig_residuals_group2];
idx=-2;
isdata=false;
idx_e=0;
idx_s=0;
for s=1:numel(sides)
    side=sides{s};
    side_label=side_labels{s};
    counter_side=countersides{s};
    idx=idx+1;
    slh(s)=text(idx+2,1,[side_label ' HEMISPACE']);
    idx_e=idx_e+1;
    idx_s=idx_s+1;
    for e=1:numel(effectors)
        idx=idx+1;
        effector=effectors{e};
        switch str2double(effector)
            case 0, effector_label='saccades'; case 1, effector_label='free gaze reaches'; case 2, effector_label='joint eye-hand'; case 3, effector_label='disociated saccades'; case 4, effector_label='dissociated reaches'; case 6, effector_label='free gaze with fixation';
        end
        %         elh(numel(effectors)*(s-1)+e)=text(idx+1,0.95,effector_label);
        type_effector=['t_' type '_e_' effector];
        
        if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
            for h=1:numel(hands)
                hand=['_' hands{h}];idx_bar_start=idx+1;
                [idx e_bar_tmp]=plot_internal_per_condition_consecutive(idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),par,stat,precision);
            end
        elseif strcmp(sac_rea,'saccades') && str2double(effector)==0
            hand=[];idx_bar_start=idx+1;
            [idx e_bar_tmp]=plot_internal_per_condition_consecutive(idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.sac,par,stat,precision);
        end
        title(type_labels)
        idx_e=idx_e+1;
    end
end



%% setting limits and ticks
y_lim=get(gca,'ylim');
for s=1:numel(slh)
    pos=get(slh(s),'Position');
    set(slh(s),'Position',[pos(1) y_lim(2)- diff(y_lim)*1/10])
end
% for e=1:numel(elh)
%     pos=get(elh(e),'Position');
%     set(elh(e),'Position',[pos(1) (y_lim(2)- diff(y_lim)*1/10)])
% end
y_lim(2)=y_lim(2)+diff(y_lim)*3/10;




ylabel(par);

end

function isdata=temp_means_bars_ch_in(input_group1,input_group2,sac_rea,type,effectors,unique_pos,Plot_settings,par,stat,sig_residuals_group1,sig_residuals_group2)
global GLO

switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
hands={'LH','RH'};
sides={'L','R'};
side_labels={'LEFT','RIGHT'};
countersides={'R','L'};
hold on

groups=[input_group1; input_group2];
sig_residuals=[sig_residuals_group1; sig_residuals_group2];
idx=-2;
isdata=false;

for s=1:numel(sides)
    side=sides{s};
    side_label=side_labels{s};
    counter_side=countersides{s};
    idx=idx+1;
    slh(s)=text(idx+2,1,[side_label ' HEMISPACE']);
    
    for e=1:numel(effectors)
        idx=idx+1;
        effector=effectors{e};
        switch str2double(effector)
            case 0, effector_label='saccades'; case 1, effector_label='free gaze reaches'; case 2, effector_label='joint eye-hand'; case 3, effector_label='disociated saccades'; case 4, effector_label='dissociated reaches'; case 6, effector_label='free gaze with fixation';
        end
        elh(numel(effectors)*(s-1)+e)=text(idx+1,1,effector_label);
        type_effector=['t_' type '_e_' effector];
        
        if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
            for h=1:numel(hands)
                hand=['_' hands{h}];idx_bar_start=idx+1;
                [idx e_bar_tmp]=plot_internal_per_condition_ch_in(idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),par,stat,sig_residuals);
                
                %                 [idx_stat_IN{idx} idx_stat_CH{idx}]=stat_sig(idx,sac_rea,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),par,stat,precision);
                e_bar(idx_bar_start:idx)=e_bar_tmp(idx_bar_start:idx);
            end
        elseif strcmp(sac_rea,'saccades') && str2double(effector)==0
            if strcmp(side,'R')
                idx = 6;
            end
            hand=[];idx_bar_start=idx+1;
            [idx e_bar_tmp]=plot_internal_per_condition_ch_in(idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.sac,par,stat,sig_residuals);
            
            %             [idx_stat_IN{idx} idx_stat_CH{idx}]=stat_sig(idx,sac_rea,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.sac,par,stat,precision);
            e_bar(idx_bar_start:idx)=e_bar_tmp(idx_bar_start:idx);
        end
        title(type_labels)
    end
end


% if strcmp(sac_rea,'reaches')
%     comp_IN={[1 2] [3 4], [6 7] [8 9], [12 13] [14 15], [17 18] [19 20]};
comp_IN={[1 2] [3 4], [7 8] [9 10]};

% else
%     comp_IN={[1 2], [4 5], [6 7], [9 10] [11 12], [15 16], [18 19] [20 21], [23 24] [25 26]};
% end
comp_IN_bar=comp_IN;


%% setting limits and ticks
y_lim=get(gca,'ylim');
for s=1:numel(slh)
    pos=get(slh(s),'Position');
    set(slh(s),'Position',[pos(1) y_lim(2)+diff(y_lim)*2/10])
end
for e=1:numel(elh)
    pos=get(elh(e),'Position');
    set(elh(e),'Position',[pos(1) y_lim(2)+diff(y_lim)*1/10])
end
y_lim(2)=y_lim(2)+diff(y_lim)*3/10;

% if strcmp(sac_rea,'reaches')
xticks=[1 2 3 4, 7 8 9 10];
%     xticklabels={'LHa' 'LHb' 'RHa' 'RHb' 'LHa' 'LHb' 'RHa' 'RHb' 'LHa' 'LHb' 'RHa' 'RHb'};
xticklabels={'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina'};
% else
%     xticks=[1 2, 4 5 6 7, 9 10 11 12, 15 16, 18 19 20 21, 23 24 25 26];
%     %     xticklabels={'NHa' 'NHb' 'LHa' 'LHb' 'RHa' 'RHb' 'LHa' 'LHb' 'RHa' 'RHb' 'NHa' 'NHb' 'LHa' 'LHb' 'RHa' 'RHb'};
%     xticklabels={'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina' 'Con' 'Ina'};
% end

set(gca,'ylim',y_lim,'xlim',[0,max(xticks)+1],'xtick',xticks,'XTickLabel',xticklabels);
ylabel(par);

hold on

% set(gca,'layer','top','ygrid','on','yminorgrid','on')

comp_IN_mat=vertcat(comp_IN_bar{:});
ho = patch([comp_IN_mat(:,1)-0.5 comp_IN_mat(:,1)+0.5 comp_IN_mat(:,1)+0.5 comp_IN_mat(:,1)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
set(ho,'facecolor',[0.85 0.85 0.85])
ha = patch([comp_IN_mat(:,2)-0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
set(ha,'facecolor',[[0.65,0.65,0.65]])
uistack(ho,'bottom');
uistack(ha,'bottom');

end

function isdata=plot_par_histograms(input_group1,input_group2,sac_rea,type,effectors,unique_pos,Plot_settings,par,precision,e,s)
global GLO
switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
switch str2double(e)
    case 0, effector_label='saccades'; case 1, effector_label='free gaze reaches'; case 2, effector_label='joint eye-hand'; case 3, effector_label='disociated saccades'; case 4, effector_label='dissociated reaches'; case 6, effector_label='free gaze with fixation';
end
hands={'LH','RH'};
sides={'L','R'};
side_labels={'LEFT','RIGHT'};
countersides={'R','L'};
hold on

groups=[input_group1; input_group2];
idx=-2;
isdata=false;
side=sides{s};
side_label=side_labels{s};
counter_side=countersides{s};
idx=idx+1;

effector=effectors{e};
switch str2double(effector)
    case 0, effector_labels='saccades'; case 1, effector_labels='free gaze reaches'; case 2, effector_labels='joint eye-hand'; case 3, effector_labels='disociated saccades'; case 4, effector_labels='dissociated reaches'; case 6, effector_labels='free gaze with fixation';
end

type_effector=['t_' type '_e_' effector];

if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
    for h=1:numel(hands)
        hand=['_' hands{h}];idx_bar_start=idx+1;
        [idx]=plot_par_histograms_internal(idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),par,precision);
    end
elseif strcmp(sac_rea,'saccades') && str2double(effector)==0
    hand=[];idx_bar_start=idx+1;
    [idx]=plot_par_histograms_internal(idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.sac,par,precision);
end
% set(gca,'xlim',[-inf,inf],'ylim',[-inf,inf])
title([type_labels ' ' effector_labels ' ' side_label])
end

function c=plot_par_correlations(input_group1,input_group2,sac_rea,type,effectors,unique_pos,Plot_settings,par,e,s)
global GLO
switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
hands={'LH','RH'};
sides={'L','R'};
side_labels={'LEFT','RIGHT'};
countersides={'R','L'};
hold on

groups=[input_group1.(par);  input_group2.(par)];

info(1,:).trial=input_group1.trial;
info(2,:).trial=input_group2.trial;
info(1,:).run=input_group1.run;
info(2,:).run=input_group2.run;
info(1,:).session=input_group1.session;
info(2,:).session=input_group2.session;

idx=-2;
side=sides{s};
side_label=side_labels{s};
counter_side=countersides{s};
idx=idx+1;

effector=effectors{e};
switch str2double(effector)
    case 0, effector_labels='saccades'; case 1, effector_labels='free gaze reaches'; case 2, effector_labels='joint eye-hand'; case 3, effector_labels='disociated saccades'; case 4, effector_labels='dissociated reaches'; case 6, effector_labels='free gaze with fixation';
end

type_effector=['t_' type '_e_' effector];
col_temp=[];
if strcmp(sac_rea,'reaches')
    for h=1:numel(hands)
        hand=['_' hands{h}];idx_bar_start=idx+1;
        plot_correlation_internal(groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),par,h,info);
        
    end
end
set(gca,'xlim',[-inf,inf],'ylim',[-inf,inf])
title([type_labels ' ' effector_labels ' ' side_label])
end

function c=errors(G, Abort_codes,sac_rea,type,effectors,unique_pos,Plot_settings,par,stat,sig_residuals_group1,sig_residuals_group2)
global GLO

switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
hands={'LH','RH'};
sides={'L','R'};
side_labels={'LEFT','RIGHT'};
countersides={'R','L'};
hold on

groups=[G(1); G(2)];
idx=-1;
isdata=false;

for s=1:numel(sides)
    side=sides{s};
    side_label=side_labels{s};
    counter_side=countersides{s};
    idx=idx+1;
    slh(s)=text(idx+2,1,[side_label ' HEMISPACE']);
    
    for e=1:numel(effectors)
        idx=idx+1;
        effector=effectors{e};
        switch str2double(effector)
            case 0, effector_label='saccades'; case 1, effector_label='free gaze reaches'; case 2, effector_label='joint eye-hand'; case 3, effector_label='disociated saccades'; case 4, effector_label='dissociated reaches'; case 6, effector_label='free gaze with fixation';
        end
        %         elh(numel(effectors)*(s-1)+e)=text(idx+1,1,effector_label);
        type_effector=['t_' type '_e_' effector];
        
        if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
            for h=1:numel(hands)
                hand=['_' hands{h}];
                idx_bar_start=idx+1;
                [idx e_bar_tmp, Abort_codes]=errors_internal(Abort_codes,idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),par,stat);
                e_bar(idx_bar_start:idx)=e_bar_tmp(idx_bar_start:idx);
            end
        elseif strcmp(sac_rea,'saccades') && str2double(effector)==0
            hand=[];idx_bar_start=idx+1;
            [idx e_bar_tmp, Abort_codes]=errors_internal(Abort_codes,idx,groups,side,counter_side,type_effector,hand,Plot_settings,Plot_settings.colors.sac,par,stat);
            e_bar(idx_bar_start:idx)=e_bar_tmp(idx_bar_start:idx);
        end
        title(type_labels)
    end
end

y_lim = [-0.5,1];

xticklabels='';
set(gca,'XTickLabel',xticklabels);

end

function  hand_release(input_group1,input_group2,sac_rea,Plot_settings,par)
global GLO
if GLO.trial_by_trial
    sub     ='raw_of_raw';
    tmean   ='mean_of_raw';
    tsem    ='std_of_raw';
else
    sub     ='raw_of_mean';
    tmean   ='mean_of_mean';
    tsem    ='sem_of_mean';
end

temp_raw = 'raw_of_mean';
LH1=input_group1.LH; RH1=input_group1.RH;
LH2=input_group2.LH; RH2=input_group2.RH;


if ~GLO.parametric_testing && GLO.calculate_statististics && GLO.plot_statististics
    [ps(1)] = signrank(LH1.(sub)',RH1.(sub)');
    [ps(2)] = signrank(LH2.(sub)',RH2.(sub)');
    [ps(3)] = signrank(LH1.(sub)',LH2.(sub)');
    [ps(4)] = signrank(RH1.(sub)',RH2.(sub)');
elseif GLO.parametric_testing && GLO.calculate_statististics && GLO.plot_statististics
    [~,     ps(1)] = ttest2(LH1.(sub)',RH1.(sub)',0.05,'both',1);
    [~,     ps(2)] = ttest2(LH2.(sub)',RH2.(sub)',0.05,'both',1);
    [~,     ps(3)] = ttest2(LH1.(sub)',LH2.(sub)',0.05,'both',1);
    [~,     ps(4)] = ttest2(RH1.(sub)',RH2.(sub)',0.05,'both',1);
end

if GLO.calculate_statististics && GLO.plot_statististics
    for g =1:size(ps,2)
        if ps(g)>=0.05
            ps(g)=NaN;
        end
    end
end

bars = {'LH1','LH2','RH1','RH2'};

subplot(2,1,1)

for idx=1:numel(bars)
    hold on
    d=eval(bars{idx});
    col1=[Plot_settings.colors.LH(1,:); Plot_settings.colors.LH(2,:);Plot_settings.colors.RH(1,:); Plot_settings.colors.RH(2,:)];
    col2=[Plot_settings.colors.IN];
    if GLO.point_per_batch
        raw_h{idx} = plot(idx,d.(temp_raw),'o','MarkerFaceColor','none','MarkerEdgeColor',col1(idx,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
    else
        e_bar{idx}=errorbar_constant_barsize_working(idx,d.(tmean),d.(tsem),d.(tsem),0.3,Plot_settings.markers.B);
        set(e_bar{idx},'MarkerFaceColor','none','MarkerEdgeColor',col1(idx,:),'Color',col1(idx,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
    end
    
    if GLO.text_in_plot
        t=text(idx+(0.1), d.(tmean)-(d.(tmean)*0), sprintf('%.3f + %.3f',d.(tmean),d.(tsem)));
        set(t,'Color', [0.3 0.3 0.3], 'FontSize', GLO.fontsize_small, 'FontWeight', 'bold')
    end
end

ylabel(par);
xticks_cell              ={[1 2], [3 4]};
xticks                  =[1, 2, 3, 4];
xticklabels         ={'Con' 'Ina' 'Con' 'Ina'};
if GLO.calculate_statististics ==1 && GLO.plot_statististics == 1
    pairs_lab           ={[1 2], [3 4], [1 3], [2 4]};
    pairs_sig            =[ps(3) ps(4) ps(1) ps(2)];
    if GLO.only_significant
        pairs_lab=pairs_lab(~isnan(pairs_sig)&pairs_sig<0.05);
        pairs_sig=pairs_sig(~isnan(pairs_sig)&pairs_sig<0.05);
    end
    sigstar_eye_hand(pairs_lab,pairs_sig,0,Plot_settings.colors.IN);
end
y_lim=get(gca,'ylim');
set(gca,'ylim',y_lim,'xlim',[0,max(xticks)+1],'xtick',xticks,'XTickLabel',xticklabels);
ylabel(par);
% set(gca,'layer','top','ygrid','on','yminorgrid','on')
comp_IN_mat=vertcat(xticks_cell{:});
ho = patch([comp_IN_mat(:,1)-0.5 comp_IN_mat(:,1)+0.5 comp_IN_mat(:,1)+0.5 comp_IN_mat(:,1)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
set(ho,'facecolor',[0.85 0.85 0.85])
ha = patch([comp_IN_mat(:,2)-0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
set(ha,'facecolor',[[0.65,0.65,0.65]])
uistack(ho,'bottom');
uistack(ha,'bottom');

subplot(2,1,2)
lin=['-'];

rt_b=0:0.1:2;
col1=[Plot_settings.colors.LH(1,:); Plot_settings.colors.LH(2,:);Plot_settings.colors.RH(1,:); Plot_settings.colors.RH(2,:)];
col2=[Plot_settings.colors.IN];
if isnan(LH1.mean_of_mean) && isnan(LH2.mean_of_mean) && isnan(RH1.mean_of_mean) && isnan(RH2.mean_of_mean), return, end

if ~GLO.CDF
    hold on
    LH1h = hist(LH1.raw_of_raw,rt_b)/nansum(LH1.raw_of_raw)*100;
    LH2h = hist(LH2.raw_of_raw,rt_b)/nansum(LH2.raw_of_raw)*100;
    RH1h = hist(RH1.raw_of_raw,rt_b)/nansum(RH1.raw_of_raw)*100;
    RH2h = hist(RH2.raw_of_raw,rt_b)/nansum(RH2.raw_of_raw)*100;
    
    hold on
    plot(rt_b,LH1h,'Color',col1(1,:),'Linewidth',GLO.linewidth,'LineStyle',lin)
    plot(rt_b,LH2h,'Color',col1(2,:),'Linewidth',GLO.linewidth,'LineStyle',lin)
    plot(rt_b,RH1h,'Color',col1(3,:),'Linewidth',GLO.linewidth,'LineStyle',lin)
    plot(rt_b,RH2h,'Color',col1(4,:),'Linewidth',GLO.linewidth,'LineStyle',lin)
    
else
    hold on
    if ~GLO.testing_patient
        if ~isempty(LH1.raw_of_raw) && ~all(isnan(LH1.raw_of_raw)), LH1c = cdfplot(LH1.raw_of_raw); set(LH1c,'Color',col1(1,:),'Linewidth',GLO.linewidth,'LineStyle',lin); end
        if ~isempty(LH2.raw_of_raw) && ~all(isnan(LH2.raw_of_raw)), LH2c = cdfplot(LH2.raw_of_raw); set(LH2c,'Color',col1(2,:),'Linewidth',GLO.linewidth,'LineStyle',lin); end
        if ~isempty(RH1.raw_of_raw) && ~all(isnan(RH1.raw_of_raw)), RH1c = cdfplot(RH1.raw_of_raw); set(RH1c,'Color',col1(3,:),'Linewidth',GLO.linewidth,'LineStyle',lin); end
        if ~isempty(RH2.raw_of_raw) && ~all(isnan(RH2.raw_of_raw)), RH2c = cdfplot(RH2.raw_of_raw); set(RH2c,'Color',col1(4,:),'Linewidth',GLO.linewidth,'LineStyle',lin); end
        grid off
    end
end
end



function  hand_release_abort(input_group1,input_group2,sac_rea,Plot_settings,par)
global GLO
% if GLO.trial_by_trial
sub     ='raw_of_raw';
tmean   ='mean_of_raw';
tsem    ='std_of_raw';
% else
%     sub     ='raw_of_mean';
%     tmean   ='mean_of_mean';
%     tsem    ='sem_of_mean';
% end
Abort_code_temp = {
    'Dummy', 1;
    'Wrong hnd' ,2;
    'Hnd Fix Acq'  ,3;
    'Hnd Fix Hol' ,4};


pre_LH1=input_group1.LH; pre_RH1=input_group1.RH;
pre_LH2=input_group2.LH; pre_RH2=input_group2.RH;

LH1 = sum(pre_LH1.(sub)~=0);
RH1 = sum(pre_RH1.(sub)~=0);
LH2 = sum(pre_LH2.(sub)~=0);
RH2 = sum(pre_RH2.(sub)~=0);

cLH1 = sum(pre_LH1.(sub)==0);
cRH1 = sum(pre_RH1.(sub)==0);
cLH2 = sum(pre_LH2.(sub)==0);
cRH2 = sum(pre_RH2.(sub)==0);

trials_LH1 = pre_LH1.num_hits_of_raw;
trials_LH2 = pre_LH2.num_hits_of_raw;

trials_RH1 = pre_RH1.num_hits_of_raw;
trials_RH2 = pre_RH2.num_hits_of_raw;

[h_LH p_LH]=f_exakt_scalars(LH1,cLH1,LH2,cLH2);
[h_RH p_RH]=f_exakt_scalars(RH1,cRH1,RH2,cRH2);

[h_LR1 p_LR1]=f_exakt_scalars(LH1,cLH1,RH1,cRH1);
[h_LR1 p_LR2]=f_exakt_scalars(LH2,cLH2,RH2,cRH2);

if p_LH>=0.05,   p_LH=NaN;    end
if p_RH>=0.05,   p_RH=NaN;    end
if p_LR1>=0.05,  p_LR1=NaN;   end
if p_LR2>=0.05,  p_LR2=NaN;   end

bars = {'LH1','LH2','RH1','RH2'};
total_bars = {'trials_LH1','trials_LH2','trials_RH1','trials_RH2'};
subplot(3,1,1)

for idx=1:numel(bars)
    hold on
    d=eval(bars{idx});
    t=eval(total_bars{idx});
    col1=[Plot_settings.colors.LH(1,:); Plot_settings.colors.LH(2,:);Plot_settings.colors.RH(1,:); Plot_settings.colors.RH(2,:)];
    col2=[Plot_settings.colors.IN];
    %     if GLO.point_per_batch
    raw_h{idx} = plot(idx,d/t,'o','MarkerFaceColor','none','MarkerEdgeColor',col1(idx,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
    %     else
    %         e_bar{idx}=errorbar_constant_barsize_working(idx,(tmean),d/t,d.(tsem),0.3,Plot_settings.markers.B);
    %         set(e_bar{idx},'MarkerFaceColor','none','MarkerEdgeColor',col1(idx,:),'Color',col1(idx,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
    %     end
    
    if GLO.text_in_plot
        t=text(idx+(0.1), d/t, sprintf('%.2f',d/t));
        set(t,'Color', [0.3 0.3 0.3], 'FontSize', GLO.fontsize_small, 'FontWeight', 'bold')
    end
end



ylabel(par);
xticks_cell              ={[1 2], [3 4]};
xticks                  =[1, 2, 3, 4];
xticklabels         ={'Con' 'Ina' 'Con' 'Ina'};
pairs_lab           ={[1 2], [3 4], [1 3], [2 4]};
pairs_sig            =[p_LH p_RH p_LR1 p_LR2];

if GLO.only_significant
    pairs_lab=pairs_lab(~isnan(pairs_sig)&pairs_sig<0.05);
    pairs_sig=pairs_sig(~isnan(pairs_sig)&pairs_sig<0.05);
end
if GLO.calculate_statististics ==1 && GLO.plot_statististics == 1
    sigstar_eye_hand(pairs_lab,pairs_sig,0,Plot_settings.colors.IN);
end
y_lim=get(gca,'ylim');
set(gca,'ylim',y_lim,'xlim',[0,max(xticks)+1],'xtick',xticks,'XTickLabel',xticklabels);
ylabel(par);
% set(gca,'layer','top','ygrid','on','yminorgrid','on')
comp_IN_mat=vertcat(xticks_cell{:});
ho = patch([comp_IN_mat(:,1)-0.5 comp_IN_mat(:,1)+0.5 comp_IN_mat(:,1)+0.5 comp_IN_mat(:,1)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
set(ho,'facecolor',[0.85 0.85 0.85])
ha = patch([comp_IN_mat(:,2)-0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)+0.5 comp_IN_mat(:,2)-0.5]',repmat([-1000 -1000 1000 1000],size(comp_IN_mat,1),1)','r');
set(ha,'facecolor',[[0.65,0.65,0.65]])
uistack(ho,'bottom');
uistack(ha,'bottom');
end

function  raw_EH_first_reach_or_saccade(input_group1,input_group2,sac_rea,type,effectors,unique_pos,Plot_settings,par,e,s)
if strcmp(sac_rea,'saccades'); return; end;
global GLO
switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
hands={'LH','RH'};
hold on
if strcmp(par,'abort_raw_x')
    groups(1,1).x         =input_group1.saccades.abort_raw_x;            groups(2,1).x         =input_group2.saccades.abort_raw_x;
    groups(1,1).y         =input_group1.saccades.abort_raw_y;            groups(2,1).y         =input_group2.saccades.abort_raw_y;
    groups(1,1).states    =input_group1.saccades.abort_raw_states;       groups(2,1).states    =input_group2.saccades.abort_raw_states;
    groups(1,1).time      =input_group1.saccades.abort_raw_time_axis;    groups(2,1).time      =input_group2.saccades.abort_raw_time_axis;
    groups(1,1).fix_pos   =input_group1.saccades.abort_fix_pos;          groups(2,1).fix_pos   =input_group2.saccades.abort_fix_pos;
    groups(1,1).tar_pos   =input_group1.saccades.abort_tar_pos;          groups(2,1).tar_pos   =input_group2.saccades.abort_tar_pos;
    groups(1,1).REACHX    =input_group1.reaches.abort_raw_x;             groups(2,1).REACHX         =input_group2.reaches.abort_raw_x;
elseif strcmp(par,'success_raw_x')
    groups(1,1).x         =input_group1.saccades.success_raw_x;            groups(2,1).x         =input_group2.saccades.success_raw_x;
    groups(1,1).y         =input_group1.saccades.success_raw_y;            groups(2,1).y         =input_group2.saccades.success_raw_y;
    groups(1,1).states    =input_group1.saccades.success_raw_states;       groups(2,1).states    =input_group2.saccades.success_raw_states;
    groups(1,1).time      =input_group1.saccades.success_raw_time_axis;    groups(2,1).time      =input_group2.saccades.success_raw_time_axis;
    groups(1,1).fix_pos   =input_group1.saccades.success_fix_pos;          groups(2,1).fix_pos   =input_group2.saccades.success_fix_pos;
    groups(1,1).tar_pos   =input_group1.saccades.success_tar_pos;          groups(2,1).tar_pos   =input_group2.saccades.success_tar_pos;
    groups(1,1).REACHX    =input_group1.reaches.success_raw_x;             groups(2,1).REACHX    =input_group2.reaches.success_raw_x;
end

idx=-2;
idx=idx+1;
effector=effectors{e};
switch str2double(effector)
    case 0, effector_labels='saccades'; case 1, effector_labels='free gaze reaches'; case 2, effector_labels='joint eye-hand'; case 3, effector_labels='disociated saccades'; case 4, effector_labels='dissociated reaches'; case 6, effector_labels='free gaze with fixation';
end

type_effector=['t_' type '_e_' effector];
if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
    for h=1:numel(hands)
        hand=['_' hands{h}];
        if GLO.instructed_only
            conditions={'IN'};
            lin=['-'];
            col2=[0.5 0.5 0.5];
        elseif GLO.choice_only
            conditions={'CH'};
            lin=[':'];
            col2=[0 0 0];
        else
            conditions={'IN','CH'};
            lin=['-',':'];
            col2=[0.5 0.5 0.5; 0 0 0];
        end
        switch str2num(type_effector(3))
            case 1
                mov_in_state = [2];
            case 2
                mov_in_state = [4,5];
            case 3
                mov_in_state = [9,10,4,5];
            case 4
                mov_in_state = [4,5];
        end
        
        for g=1:size(groups,1)
            current_group=groups(g,:);
            %     id_x=id_x+0.1;
            for d=1:numel(conditions)
                con=conditions{d};
                condition=[type_effector  '_' con hand];
                hold on
                if ~isfield(current_group.x,condition)
                    continue;
                end
                sac_before_reach.trace_idx=[];
                sac_before_reach.trace_idx_L=[];
                sac_before_reach.trace_idx_R=[];
                for jj=1:numel(current_group.x.(condition))
                    if ~iscell(current_group.states.(condition)(jj))
                        continue
                    end
                    idx_in_state = ismember(current_group.states.(condition){jj},mov_in_state);
                    fix_pos=current_group.fix_pos.(condition)(jj);
                    tar_pos=current_group.tar_pos.(condition)(jj);
                    trace_tmp=current_group.x.(condition){jj}(idx_in_state)+1i*current_group.y.(condition){jj}(idx_in_state);
                    
                    if real(tar_pos) > 0.01
                        tar_pos_R = current_group.tar_pos.(condition)(jj);
                        tar_pos_L = NaN + 1i*NaN;
                    elseif real(tar_pos) < 0.01
                        tar_pos_L = current_group.tar_pos.(condition)(jj);
                        tar_pos_R = NaN + 1i*NaN;
                    else
                        tar_pos_R = NaN + 1i*NaN;
                        tar_pos_L = NaN + 1i*NaN;
                        tar_pos = NaN + 1i*NaN;
                    end
                    
                    trace_idx=abs(trace_tmp-tar_pos)<5; %%%% radius==5 !?
                    trace_idx_R = abs(trace_tmp-tar_pos_R)<5;
                    trace_idx_L = abs(trace_tmp-tar_pos_L)<5;
                    space_divided={'trace_idx','trace_idx_L','trace_idx_R'};
                    tar_pos_current = {'tar_pos','tar_pos_L','tar_pos_R'};
                    
                    for s_d = 1:numel(space_divided)
                        
                        if isnan(eval(tar_pos_current{s_d}))
                            sac_before_reach.(space_divided{s_d})(jj)=NaN;
                            continue
                        end
                        
                        
                        Inside=find(diff(eval(space_divided{s_d}))>0);
                        Outside=find(diff((eval(space_divided{s_d})))<0);
                        if isempty(Inside)
                            sac_before_reach.(space_divided{s_d})(jj)=0;
                            continue
                        end
                        if ~isempty(Outside) && Inside(1)>Outside(1)
                            Outside(1)=[];
                        end
                        if isempty(Outside) || Inside(end)>Outside(end)
                            Outside(end+1)=Inside(end)+1000;
                        end
                        t_inside=find((Outside-Inside)>50);
                        
                        notouch_idx=isnan(current_group.REACHX.(condition){jj}(idx_in_state));
                        reach_in_ms=find(diff(notouch_idx)<0,1,'last'); %first or last? diff<0 touch, > release
                        if ~isempty(reach_in_ms) && any(Inside(t_inside) < reach_in_ms)
                            sac_before_reach.(space_divided{s_d})(jj)=1;
                        elseif ~isempty(reach_in_ms) && any(Inside(t_inside) >= reach_in_ms)
                            sac_before_reach.(space_divided{s_d})(jj)=-1;
                        else
                            sac_before_reach.(space_divided{s_d})(jj)=0;
                        end
                        clear Inside Outside
                    end
                    
                end
                n=sum(~isnan(sac_before_reach.trace_idx));
                n_L=sum(~isnan(sac_before_reach.trace_idx_L));
                n_R=sum(~isnan(sac_before_reach.trace_idx_R));
                
                %             text(-75+40*d,g*20+h*10-40,sprintf('g %d, %s S<R: %.0d  S>R:%.0d  No:%.0d  ',g,condition,(sum(sac_before_reach.trace_idx==1)/n)*100,(sum(sac_before_reach.trace_idx==-1)/n)*100,(sum(sac_before_reach.trace_idx==0)/n)*100),'interpreter','none');
                if g==1
                    text(-75+40*d,g*20+h*10-40+5,sprintf('g_LS %d, %s S<R: %.0f  S>R:%.0f  No:%.0f  ',g,condition,(sum(sac_before_reach.trace_idx_L==1)/n_L)*100,(sum(sac_before_reach.trace_idx_L==-1)/n_L)*100,(sum(sac_before_reach.trace_idx_L==0)/n_L)*100),'interpreter','none');
                    text(-75+40*d,g*20+h*10-40+10,sprintf('g_RS %d, %s S<R: %.0f  S>R:%.0f  No:%.0f  ',g,condition,(sum(sac_before_reach.trace_idx_R==1)/n_R)*100,(sum(sac_before_reach.trace_idx_R==-1)/n_R)*100,(sum(sac_before_reach.trace_idx_R==0)/n_R)*100),'interpreter','none');
                end
            end
        end
    end
end

set(gca,'xlim',[-30 30],'ylim',[-30 30])
title([type_labels ' ' effector_labels ])
end


function  raw_plotting(input_group1,input_group2,sac_rea,type,effectors,unique_pos,Plot_settings,par,e,s,target_radius, target_size, fixation_size, fixation_radius)
global GLO
switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
hands={'LH','RH'};
hold on
if strcmp(par,'abort_raw_x')
    groups(1,1).x         =input_group1.abort_raw_x;            groups(2,1).x         =input_group2.abort_raw_x;
    groups(1,1).y         =input_group1.abort_raw_y;            groups(2,1).y         =input_group2.abort_raw_y;
    groups(1,1).states    =input_group1.abort_raw_states;       groups(2,1).states    =input_group2.abort_raw_states;
    groups(1,1).time      =input_group1.abort_raw_time_axis;    groups(2,1).time      =input_group2.abort_raw_time_axis;
    groups(1,1).fix_pos   =input_group1.abort_fix_pos;          groups(2,1).fix_pos   =input_group2.abort_fix_pos;
elseif strcmp(par,'success_raw_x')
    groups(1,1).x         =input_group1.success_raw_x;            groups(2,1).x         =input_group2.success_raw_x;
    groups(1,1).y         =input_group1.success_raw_y;            groups(2,1).y         =input_group2.success_raw_y;
    groups(1,1).states    =input_group1.success_raw_states;       groups(2,1).states    =input_group2.success_raw_states;
    groups(1,1).time      =input_group1.success_raw_time_axis;    groups(2,1).time      =input_group2.success_raw_time_axis;
    groups(1,1).fix_pos   =input_group1.success_fix_pos;          groups(2,1).fix_pos   =input_group2.success_fix_pos;
end

idx=-2;
idx=idx+1;
effector=effectors{e};
switch str2double(effector)
    case 0, effector_labels='saccades'; case 1, effector_labels='free gaze reaches'; case 2, effector_labels='joint eye-hand'; case 3, effector_labels='disociated saccades'; case 4, effector_labels='dissociated reaches'; case 6, effector_labels='free gaze with fixation';
end

%plot target
if ~isempty(unique_pos)
    angle=0:0.001:2*pi();
    for p=1:numel(unique_pos)
        xcircle_rad=target_radius*cos(angle);
        ycircle_rad=target_radius*sin(angle);
        plot(real(unique_pos(p))+xcircle_rad,imag(unique_pos(p))+ycircle_rad,'--r');
        xcircle_siz=target_size*cos(angle);
        ycircle_siz=target_size*sin(angle);
        plot(real(unique_pos(p))+xcircle_siz,imag(unique_pos(p))+ycircle_siz,'r');
        
    end
    %plot fixation
    xcircle_fix_rad=fixation_radius*cos(angle);
    ycircle_fix_rad=fixation_radius*sin(angle);
    plot(0+xcircle_fix_rad,0+ycircle_fix_rad,'--r');
    xcircle_fix_siz=fixation_size*cos(angle);
    ycircle_fix_siz=fixation_size*sin(angle);
    plot(0+xcircle_fix_siz,0+ycircle_fix_siz,'r');
    
    type_effector=['t_' type '_e_' effector];
    
    
    if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
        for h=1:numel(hands)
            hand=['_' hands{h}];
            plot_raw_internal(groups,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),sac_rea);
        end
    elseif strcmp(sac_rea,'saccades') && str2double(effector)==0
        hand=[];
        plot_raw_internal(groups,type_effector,hand,Plot_settings,Plot_settings.colors.sac),sac_rea;
    end
    
end

set(gca,'xlim',[-30,30],'ylim',[-30,30])
title([type_labels ' ' effector_labels ])
end

function  raw_plotting_xy(input_group1,input_group2,sac_rea,type,effectors,unique_pos,Plot_settings,par,e,s,current_axis,target_radius, target_size)
global GLO
switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
hands={'LH','RH'};
hold on

if strcmp(par,'abort_raw_x')
    groups(1,1).x         =input_group1.abort_raw_x;            groups(2,1).x         =input_group2.abort_raw_x;
    groups(1,1).y         =input_group1.abort_raw_y;            groups(2,1).y         =input_group2.abort_raw_y;
    groups(1,1).states    =input_group1.abort_raw_states;       groups(2,1).states    =input_group2.abort_raw_states;
    groups(1,1).time      =input_group1.abort_raw_time_axis;    groups(2,1).time      =input_group2.abort_raw_time_axis;
    groups(1,1).fix_pos   =input_group1.abort_fix_pos;          groups(2,1).fix_pos   =input_group2.abort_fix_pos;
    groups(1,1).lat       =input_group1.abort_lat;              groups(2,1).lat       =input_group2.abort_lat;
elseif strcmp(par,'success_raw_x')
    groups(1,1).x         =input_group1.success_raw_x;            groups(2,1).x         =input_group2.success_raw_x;
    groups(1,1).y         =input_group1.success_raw_y;            groups(2,1).y         =input_group2.success_raw_y;
    groups(1,1).states    =input_group1.success_raw_states;       groups(2,1).states    =input_group2.success_raw_states;
    groups(1,1).time      =input_group1.success_raw_time_axis;    groups(2,1).time      =input_group2.success_raw_time_axis;
    groups(1,1).fix_pos   =input_group1.success_fix_pos;          groups(2,1).fix_pos   =input_group2.success_fix_pos;
    groups(1,1).lat       =input_group1.success_lat;              groups(2,1).lat       =input_group2.success_lat;
end

idx=-2;
idx=idx+1;
effector=effectors{e};
switch str2double(effector)
    case 0, effector_labels='saccades'; case 1, effector_labels='free gaze reaches'; case 2, effector_labels='joint eye-hand'; case 3, effector_labels='disociated saccades'; case 4, effector_labels='dissociated reaches'; case 6, effector_labels='free gaze with fixation';
end



%     angle=0:0.001:2*pi();
%     for p=1:numel(unique_pos)
%         xcircle_rad=target_radius*cos(angle);
%         ycircle_rad=target_radius*sin(angle);
%         plot(real(unique_pos(p))+xcircle_rad,imag(unique_pos(p))+ycircle_rad,'--r');
%         xcircle_siz=target_size*cos(angle);
%         ycircle_siz=target_size*sin(angle);
%         plot(real(unique_pos(p))+xcircle_siz,imag(unique_pos(p))+ycircle_siz,'r');
%     end


type_effector=['t_' type '_e_' effector];



if strcmp(sac_rea,'reaches') || (strcmp(sac_rea,'saccades') && str2double(effector)~=0)
    for h=1:numel(hands)
        hand=['_' hands{h}];
        plot_raw_internal_xy(groups,type_effector,hand,Plot_settings,Plot_settings.colors.(hands{h}),current_axis,sac_rea);
    end
elseif strcmp(sac_rea,'saccades') && str2double(effector)==0
    hand=[];
    plot_raw_internal_xy(groups,type_effector,hand,Plot_settings,Plot_settings.colors.sac,current_axis,sac_rea);
end




set(gca,'xlim',[-30 30],'ylim',[-30 30])
title([type_labels ' ' effector_labels ])
end

%% INTERNAL LOOP PLOTTING FUNCTIONS

function [idx, e_bar]=plot_internal_per_condition(idx,groups,side,counter_side,type_effector,hand,Plot_settings,col1,par,stat,sig_residuals,precision)
global GLO
if GLO.instructed_only
    conditions={'IN'};
elseif GLO.choice_only
    conditions={'CH'};
else
    conditions={'IN','CH'};
end

for g=1:size(groups,1)
    current_group=groups(g,:);
    
    for p=1:numel(current_group)
        idx=idx+1;
        e_bar{idx}=[];
        for d=1:numel(conditions)
            con=conditions{d};
            col2=Plot_settings.colors.(con);
            condition=[type_effector  '_' side '_' con hand];
            counter_condition=[type_effector  '_' counter_side '_' con hand];
            if ~isfield(current_group(p),condition) || ~isfield(current_group(p),counter_condition)
                continue;
            end
            if GLO.trial_by_trial
                to_plot=[current_group(p).(condition).mean_of_raw];
            else
                to_plot=[current_group(p).(condition).mean_of_mean];
            end
            if all(isnan(to_plot))
                continue;
            else
                isdata=true;
            end
            
            if GLO.trial_by_trial && strcmp(par,'successful')
                temp_mean = current_group(p).(condition).mean_of_raw;
                temp_sem = 0;
            elseif GLO.trial_by_trial && ~strcmp(par,'successful')
                temp_mean = current_group(p).(condition).mean_of_raw;
                temp_sem = current_group(p).(condition).std_of_raw;
                if precision
                    temp_mean = current_group(p).(condition).mean_of_std;
                    temp_sem = current_group(p).(condition).std_of_std;
                end
            else
                temp_mean = current_group(p).(condition).mean_of_mean;
                temp_sem = current_group(p).(condition).sem_of_mean;
                if precision
                    temp_mean = current_group(p).(condition).mean_of_std;
                    temp_sem = current_group(p).(condition).sem_of_std;
                end
                
                %% for precision
                
            end
            
            
            temp_raw = current_group(p).(condition).raw_of_mean;
            
            if precision
                temp_raw = current_group(p).(condition).raw_of_std;
            end
            
            if findstr(par, 'residuals'), temp_sig = sig_residuals(p).(condition).raw_of_mean; end
            if GLO.point_per_batch
                for s=1:numel(temp_raw)
                    if strcmp('lat_r_residuals',par) && temp_sig(s) < 0.05
                        raw_h{idx}.(con) = plot(idx,temp_raw(s),'^','MarkerFaceColor',col2,'MarkerEdgeColor',col1(g,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
                    else
                        raw_h{idx}.(con) = plot(idx,temp_raw(s),'o','MarkerFaceColor',col2,'MarkerEdgeColor',col1(g,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
                    end
                end
                
            else
                e_bar{idx}.(con)=errorbar_constant_barsize_working(idx,temp_mean,temp_sem,temp_sem,0.3,Plot_settings.markers.B);
                set(e_bar{idx}.(con),'MarkerFaceColor',col2,'MarkerEdgeColor',col2,'Color',col1(g,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
            end
            if GLO.text_in_plot
                t=text(idx+(0.3), temp_mean-(temp_mean*0), sprintf('%.3f + %.3f',temp_mean,temp_sem));
                set(t,'Color', [0.3 0.3 0.3], 'FontSize', GLO.fontsize_small, 'FontWeight', 'bold')
            end
            if GLO.hits_in_plot
                n_hits=current_group(p).(condition).num_hits_of_raw;
                t=text(idx+(0.1), temp_mean-(temp_mean*0), sprintf('%.0f',n_hits));
                set(t,'Color', [0.3 0.3 0.3], 'FontSize', GLO.fontsize_small, 'FontWeight', 'bold')
            end
        end
        
    end
    
end
% line between points
if strcmp(type_effector,'t_4_e_0')  && size(groups,1) == 2  && GLO.point_per_batch;
    condition_in = [type_effector  '_' side '_IN'];
    condition_ch = [type_effector  '_' side '_CH'];
    
    for cnd = 1:size(groups(1).(condition_in).raw_of_mean,2)
        plot([idx-1  idx],[groups(1).(condition_in).raw_of_mean(cnd) groups(2).(condition_in).raw_of_mean(cnd)],'color',[0.5 0.5 0.5]);
        plot([idx-1  idx],[groups(1).( condition_ch).raw_of_mean(cnd) groups(2).( condition_ch).raw_of_mean(cnd)],'k');
    end
    
end
end

function [idx, e_bar]=plot_internal_per_condition_consecutive(idx,groups,side,counter_side,type_effector,hand,Plot_settings,col1,par,stat,precision)
global GLO
if GLO.instructed_only
    conditions={'IN'};
elseif GLO.choice_only
    conditions={'CH'};
else
    conditions={'IN','CH'};
end


for d=1:numel(conditions)
    for p=1
        idx=idx+1;
        e_bar{idx}=[];
        for g=1:size(groups,1)
            current_group=groups(g,:);
            con=conditions{d};
            col2=Plot_settings.colors.(con);
            condition=[type_effector  '_' side '_' con hand];
            counter_condition=[type_effector  '_' counter_side '_' con hand];
            if ~isfield(current_group(p),condition) || ~isfield(current_group(p),counter_condition)
                continue;
            end
            if GLO.trial_by_trial
                to_plot=[current_group(p).(condition).mean_of_raw];
            else
                to_plot=[current_group(p).(condition).mean_of_mean];
            end
            if all(isnan(to_plot))
                continue;
            else
                isdata=true;
            end
            
            if GLO.trial_by_trial && strcmp(par,'successful')
                temp_mean = current_group(p).(condition).mean_of_raw;
                temp_sem = 0;
            elseif GLO.trial_by_trial && ~strcmp(par,'successful')
                temp_mean = current_group(p).(condition).mean_of_raw;
                temp_sem = current_group(p).(condition).std_of_raw;
                if precision
                    temp_mean = current_group(p).(condition).raw_of_std;
                    temp_sem = current_group(p).(condition).std_of_std;
                end
            else
                temp_mean = current_group(p).(condition).mean_of_mean;
                temp_sem = current_group(p).(condition).sem_of_mean;
                if precision
                    temp_mean = current_group(p).(condition).raw_of_std;
                    temp_sem = current_group(p).(condition).sem_of_std;
                end
            end
            temp_raw = current_group(p).(condition).raw_of_mean;
            if precision
                temp_raw = current_group(p).(condition).raw_of_std;
            end
            if GLO.point_per_batch
                for s=1:numel(temp_raw)
                    raw_h{idx}.(con) = plot(idx,temp_raw(s),'o','MarkerFaceColor',col2,'MarkerEdgeColor',col1(g,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
                    idx=idx+1;
                    e_bar{idx}=[];
                end
            else
                e_bar{idx}.(con)=errorbar_constant_barsize_working(idx,temp_mean,temp_sem,temp_sem,0.3,Plot_settings.markers.B);
                set(e_bar{idx}.(con),'MarkerFaceColor',col2,'MarkerEdgeColor',col2,'Color',col1(g,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
            end
            if GLO.text_in_plot && ~GLO.point_per_batch
                t=text(idx+(0.1), temp_mean-(temp_mean*0), sprintf('%.5f',temp_mean));
                set(t,'Color', [0.3 0.3 0.3], 'FontSize', GLO.fontsize_small, 'FontWeight', 'bold')
            end
        end
    end
    
end
end

function [idx e_bar]=plot_internal_per_condition_ch_in(idx,groups,side,counter_side,type_effector,hand,Plot_settings,col1,par,stat,sig_residuals)
global GLO


for g=1:size(groups,1)
    current_group=groups(g,:);
    for p=1:numel(current_group)
        idx=idx+1;
        e_bar{idx}=[];
        col2=Plot_settings.colors.IN;
        condition_in=[type_effector  '_' side '_' 'IN' hand];
        condition_ch=[type_effector  '_' side '_' 'CH' hand];
        
        if ~isfield(current_group(p),condition_in) &&  ~isfield(current_group(p),condition_ch)%!!!! CHECK might have to remove the counterpart, see how affect plotting
            continue;
        end
        if GLO.trial_by_trial
            to_plot=[current_group(p).(condition_ch).mean_of_raw]-[current_group(p).(condition_in).mean_of_raw];
        else
            to_plot=[current_group(p).(condition_ch).mean_of_mean]-[current_group(p).(condition_in).mean_of_mean];
        end
        
        if GLO.trial_by_trial && strcmp(par,'successful')
            temp_mean = current_group(p).(condition_ch).mean_of_raw-current_group(p).(condition_in).mean_of_raw;
            temp_sem = 0;
        elseif GLO.trial_by_trial && ~strcmp(par,'successful')
            temp_mean = current_group(p).(condition_ch).mean_of_raw-current_group(p).(condition_in).mean_of_raw;
            temp_sem = current_group(p).(condition_ch).std_of_raw-current_group(p).(condition_in).std_of_raw;
        else
            temp_mean = current_group(p).(condition_ch).mean_of_mean-current_group(p).(condition_in).mean_of_mean;
            temp_sem = current_group(p).(condition_ch).sem_of_mean-current_group(p).(condition_in).sem_of_mean;
        end
        temp_raw = current_group(p).(condition_ch).raw_of_mean-current_group(p).(condition_in).raw_of_mean;
        
        if GLO.point_per_batch
            for s=1:numel(temp_raw)
                raw_h{idx} = plot(idx,temp_raw(s),'o','MarkerFaceColor',col2,'MarkerEdgeColor',col1(g,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
            end
        else
            e_bar{idx}=errorbar_constant_barsize_working(idx,temp_mean,temp_sem,temp_sem,0.3,Plot_settings.markers.B);
            set(e_bar{idx},'MarkerFaceColor',col2,'MarkerEdgeColor',col2,'Color',col1(g,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
        end
        if GLO.text_in_plot
            t=text(idx+(0.1), temp_mean-(temp_mean*0), sprintf('%.5f',temp_mean));
            set(t,'Color', [0.3 0.3 0.3], 'FontSize', GLO.fontsize_small, 'FontWeight', 'bold')
        end
        %         end
    end
    
end
end

function [idx] = plot_accuracy_internal_ellipse(input_group1,input_group2,sac_rea,type,effector,unique_pos,target_radius,target_size,Plot_settings,par,stat, testing,precision)
global GLO;

switch str2double(type)
    case 1, type_labels='fixation';    case 2, type_labels='visually-guided';    case 2.5, type_labels='memory learning';    case 3, type_labels='memory-guided'; case 4, type_labels='delayed visually-guided';
end
switch str2double(effector)
    case 0, effector_labels='saccades'; case 1, effector_labels='free gaze reaches'; case 2, effector_labels='joint eye-hand'; case 3, effector_labels='disociated saccades'; case 4, effector_labels='dissociated reaches'; case 6, effector_labels='free gaze with fixation';
end

input_fieldnames=fieldnames(input_group1);
hold on
angles=[0:pi/100:2*pi];
circle_x=cos(angles);
circle_y=sin(angles);
if GLO.instructed_only
    conditions= {'IN_LH','IN_RH'};
    symbol=['o','o'];
    lin=['-','-'];
elseif GLO.choice_only
    conditions= {'CH_LH','CH_RH'};
    symbol=['s','s'];
    lin=[':',':'];
else
    conditions= {'IN_LH','IN_RH','CH_LH','CH_RH'};
    symbol=['o','o','s','s'];
    lin=['-','-',':',':'];
end

Plot_settings.condition_colors=[0 1 1; 0 1 0; 0 1 1; 0 1 0; 0 0 1; 0 0.5 0; 0 0 1; 0 0.5 0;];
if strcmp(sac_rea,'saccades') && str2double(effector)==0
    hand='';
    if GLO.instructed_only
        conditions={'IN'};
        symbol='o';
        lin='-';
    elseif GLO.choice_only
        conditions={'CH'};
        symbol='s';
        lin=':';
    else
        conditions={'IN','CH'};
        symbol=['o','s'];
        lin=['-',':'];
    end
    Plot_settings.condition_colors=[1 0 0; 1 0 0;0.5 0 0; 0.5 0 0];
end

%%
% t_1 = GetSecs;
temp_struct_1 = create_nan_structure(input_group1);
temp_struct_2 = create_nan_structure(input_group2);
temp_full_struct2 = catstruct(temp_struct_1, input_group2);
temp_full_struct1 = catstruct(temp_struct_2, input_group1);
groups=[temp_full_struct1; temp_full_struct2];
% t_2 = GetSecs - t_1;

for g=1:size(groups,1)
    current_group=groups(g,:);
    idx=0;
    for p=1:numel(current_group)
        %         p=1:numel(current_group)
        idx=idx+1;
        center=[real(unique_pos(p)) imag(unique_pos(p))];
        type_effector=['t_' type '_e_' effector];
        for c=1:numel(conditions)
            condition=conditions{c};
            if g==2 && strcmp(sac_rea,'saccades') && str2double(effector)==0
                current_color=Plot_settings.condition_colors(c+2,:);
            elseif g==2
                current_color=Plot_settings.condition_colors(c+4,:);
            else
                current_color=Plot_settings.condition_colors(c,:);
            end
            type_effector_condition=[type_effector '_' condition];
            if (~ismember(type_effector_condition,input_fieldnames) || ~isfield(current_group(p).(type_effector_condition), 'mean_of_raw') && ~GLO.trial_by_trial) || (~ismember(type_effector_condition,input_fieldnames) || ~isfield(current_group(p).(type_effector_condition), 'mean_of_mean') && GLO.trial_by_trial)
                continue;
            end
            
            if GLO.trial_by_trial
                if isnan(current_group(p).(type_effector_condition).mean_of_raw)
                    continue;
                end
            else
                if isnan(current_group(p).(type_effector_condition).mean_of_mean)
                    continue;
                end
            end
            
            plot(center(1),center(2),'MarkerSize',10,'Marker','+','MarkerEdgeColor','r','MarkerFaceColor','r');
            x_tar_rad = target_radius*circle_x + center(1);
            y_tar_rad = target_radius*circle_y + center(2);
            x_tar_siz = target_size*circle_x + center(1);
            y_tar_siz = target_size*circle_y + center(2);
            plot(x_tar_rad,y_tar_rad,'--r')
            plot(x_tar_siz,y_tar_siz,'r')
            
            if GLO.trial_by_trial
                complex_mean=current_group(p).(type_effector_condition).mean_of_raw;
                complex_raw=current_group(p).(type_effector_condition).raw_of_raw;
                real_std=real(current_group(p).(type_effector_condition).std_of_raw);
                imag_std=imag(current_group(p).(type_effector_condition).std_of_raw);
            else
                complex_mean=current_group(p).(type_effector_condition).mean_of_mean;
                complex_raw=current_group(p).(type_effector_condition).raw_of_mean;
                real_std=mean(real(current_group(p).(type_effector_condition).raw_of_std));
                imag_std=mean(imag(current_group(p).(type_effector_condition).raw_of_std));
            end
            
            complex_all=current_group(p).(type_effector_condition).raw_of_raw;
            
            X_mean=real(complex_mean);
            Y_mean=imag(complex_mean);
            plot(real(complex_mean),imag(complex_mean),'MarkerSize',6,'Marker',symbol(c),'MarkerEdgeColor',current_color,'Color',current_color);
            if GLO.plot_raw_endpoints
                %                 scatter(real(complex_all),imag(complex_all),10,'Marker',symbol(c),'MarkerEdgeColor',current_color,'MarkerFaceColor',current_color);
                scatter(real(complex_all),imag(complex_all),10,'Marker',symbol(c),'MarkerEdgeColor',current_color,'MarkerFaceColor',current_color);
            end
            cl=line((circle_x*real_std)+real(complex_mean),(circle_y*imag_std)+imag(complex_mean));
            
            set(cl,'Color',current_color,'LineWidth',GLO.linewidth,'LineStyle',lin(c));
            
            if GLO.calculate_statististics ==1 && GLO.plot_statististics == 1 && isfield(stat(p).([type_effector '_' condition]),'raw_of_raw')
                if GLO.trial_by_trial && ~strcmp(testing,'patient')
                    px=real(stat(p).([type_effector '_' condition]).raw_of_raw{4});
                    py=imag(stat(p).([type_effector '_' condition]).raw_of_raw{4});
                    px_pre=real(stat(p).([type_effector '_' condition]).raw_of_std{4});
                    py_pre=imag(stat(p).([type_effector '_' condition]).raw_of_std{4});
                else
                    px=real(stat(p).([type_effector '_' condition]).raw_of_mean{4});
                    py=imag(stat(p).([type_effector '_' condition]).raw_of_mean{4});
                    px_pre=real(stat(p).([type_effector '_' condition]).raw_of_std{4});
                    py_pre=imag(stat(p).([type_effector '_' condition]).raw_of_std{4});
                end
                
                n_stars_x= (px<0.001) + (px<0.01) + (px<0.05);
                n_stars_y= (py<0.001) + (py<0.01) + (py<0.05);
                xstars = repmat('x',1,n_stars_x);
                ystars = repmat('y',1,n_stars_y);
                text(center(1),center(2)+2+c,xstars,'fontsize',20,'color',current_color)
                text(center(1),center(2)+6+c,ystars,'fontsize',20,'color',current_color)
                
                n_stars_x_pre= (px_pre<0.001) + (px_pre<0.01) + (px_pre<0.05);
                n_stars_y_pre= (py_pre<0.001) + (py_pre<0.01) + (py_pre<0.05);
                xstars_pre = repmat('-',1,n_stars_x_pre);
                ystars_pre = repmat('|',1,n_stars_y_pre);
                text(center(1),center(2)+2+c,xstars_pre,'fontsize',20,'color',current_color)
                text(center(1),center(2)+6+c,ystars_pre,'fontsize',20,'color',current_color)
            end
        end
    end
end
set(gca,'fontsize', GLO.fontsize_ticks);
xlabel('Horizontal eye position [deg]');
ylabel('Vertical eye position [deg]');
title([type_labels, ' ', effector_labels]);
set(gca,'xtick',[-100:6:100],'ytick',[-100:6:100])
% axis equal
end

function [idx]=plot_par_histograms_internal(idx,groups,side,counter_side,type_effector,hand,Plot_settings,col1,par,precision)
global GLO
if GLO.instructed_only
    conditions={'IN'};
    lin=['-'];
elseif GLO.choice_only
    conditions={'CH'};
    lin=[':'];
else
    conditions={'IN','CH'};
    lin=['-',':'];
end
switch par
    case 'lat', rt_b=0:0.05:0.8;          case 'dur', rt_b=0:0.1:0.5;            case 'n_obs', rt_b=0:1:1000;                case 'velocity', rt_b=0:10:1000;
    case 'lat_r', rt_b=0:0.1:1;             case 'lat_r_residuals', rt_b=0:0.1:1;   case 'dur_r', rt_b=0:0.1:1;                 case 'dur_r_residuals', rt_b=0:0.1:1;
    case 'lat_raw_sac_rea', rt_b=0:0.1:1;   case 'dur_raw_sac_rea', rt_b=0:0.1:1;   otherwise, rt_b=0:0.4:6;
end

for g=1:size(groups,1)
    current_group=groups(g,:);
    for p=1:numel(current_group)
        idx=idx+1;
        for d=1:numel(conditions)
            con=conditions{d};
            condition=[type_effector  '_' side '_' con hand];
            counter_condition=[type_effector  '_' counter_side '_' con hand];
            if ~isfield(current_group(p),condition) || ~isfield(current_group(p),counter_condition)
                continue;
            end
            to_plot=[current_group(p).(condition).mean_of_raw];
            if all(isnan(to_plot))
                continue;
            end
            if ~GLO.CDF
                %                 temp_hist = hist(current_group(p).(condition).raw_of_raw,rt_b)/nansum(current_group(p).(condition).raw_of_raw)*100;
                temp_hist = hist(current_group(p).(condition).raw_of_raw,rt_b);
                if precision
                    %                     temp_hist = hist(current_group(p).(condition).raw_of_std,rt_b)/nansum(current_group(p).(condition).raw_of_std)*100;
                    temp_hist = hist(current_group(p).(condition).raw_of_std,rt_b);
                end
                plot(rt_b,temp_hist,'Color',col1(g,:),'Linewidth',GLO.linewidth,'LineStyle',lin(d))
            else
                
                temp_cdf = current_group(p).(condition).raw_of_raw;
                if precision
                    temp_cdf = current_group(p).(condition).raw_of_std;
                end
                temp_cdf_h = cdfplot(temp_cdf);
                grid off
                set(temp_cdf_h,'Color',col1(g,:),'Linewidth',GLO.linewidth,'LineStyle',lin(d))
            end
        end
    end
end

end

function plot_correlation_internal(groups,side,counter_side,type_effector,hand,Plot_settings,col1,par,h,info)
global GLO

if GLO.instructed_only
    conditions={'IN'};
    mar={'o'};
    lin=['-'];
    col2=[0.5 0.5 0.5];
elseif GLO.choice_only
    conditions={'CH'};
    mar={'s'};
    lin=[':'];
    col2=[0 0 0];
else
    conditions={'IN','CH'};
    mar={'o','s'};
    lin=['-',':'];
    col2=[0.5 0.5 0.5; 0 0 0];
end

id_x=-0.1;
id_y=h*0.03;
col_temp=[];
for g=1:size(groups,1)
    current_group=groups(g,:);
    current_info=info(g,:);
    %     id_x=id_x+0.1;
    
    for d=1:numel(conditions)
        con=conditions{d};
        condition=[type_effector  '_' side '_' con hand];
        counter_condition=[type_effector  '_' counter_side '_' con hand];
        if ~isfield(current_group,condition) || ~isfield(current_group,counter_condition)
            continue;
        end
        to_plot=[current_group.(condition).mean_of_raw];
        if all(isnan(to_plot))
            continue;
        end
        temp_sac = current_group.(condition).raw_of_raw(1,:);
        temp_rea = current_group.(condition).raw_of_raw(2,:);
        
        pre_trial       = current_info.trial.(condition);
        pre_run         = current_info.run.(condition);
        pre_session     = current_info.session.(condition);
        
        info_trial      = pre_trial.raw_of_raw(1,:)';
        info_run        = pre_run.raw_of_raw(1,:)';
        info_session    = pre_session.raw_of_raw(1,:)';
        
        %             [c_R,c_P] = corr(temp_rea',temp_sac','type',GLO.correlation_mode);
        %             [c_R,c_P] = corr(temp_rea',temp_sac','type',GLO.correlation_mode);
        if ~isempty(findstr(type_effector,'e_6')) || ~isempty(findstr(type_effector,'e_1')) || ~isempty(findstr(type_effector,'e_2'))
            [slo, int, c_R, c_P] = beh_myregr_eye_hand(temp_rea',temp_sac',0,GLO.remove_outliers);
            
            
            
            
            xtmp=[temp_rea' ones(length(temp_rea'),1)]; %input matrix for regress function
            ytmp=temp_sac';
            
            %regression coefficients
            [~,~,~,Rint] = regress(ytmp,xtmp);
            
            %check the presence of outliers
            outl=find(ismember(sign(Rint),[-1 1],'rows')==0);
            if ~isempty(outl)
                if GLO.remove_outliers
                    ytmp(outl)=[]; xtmp(outl,:)=[]; info_trial(outl)=[]; info_run(outl,:)=[]; info_session(outl)=[];
                end
            end
            
            xtmp(:,2)=[]; %delete column 2
            %save coefficients value
            temp_rea        = xtmp';
            temp_sac        = ytmp';
            info_trial      = info_trial';
            info_run        = info_run';
            info_session    = info_session';
            
            plot(temp_rea,temp_sac,mar{d},'MarkerFaceColor',col2(d,:),'MarkerEdgeColor',col1(g,:),'Color',col1(g,:));
            % sc = scatter(temp_rea,temp_sac,50,'filled','MarkerEdgeColor',col1(g,:),'MarkerFaceColor',col2(d,:),'linewidth',1)
            
            
            % hold on
            % x_t=get(gca,'XTick');
            % y_t=get(gca,'Ylim');
            %
            % x_t_b=min(x_t):0.025:max(x_t);
            % y_t_b=min(y_t):0.025:max(y_t);
            %
            % hr = hist(temp_rea,x_t_b);
            % hs = hist(temp_sac,y_t_b);
            % hr_n = hr/sum(hr);
            % hs_n = hs/sum(hs);
            % plot(x_t_b,(hr_n/25)+y_t(1),'LineStyle',lin(d,:),'Color',col1(g,:))
            % plot((hs_n/25)+x_t(1),y_t_b,'LineStyle',lin(d,:),'Color',col1(g,:))
            %
            
            
            if GLO.trial_numbers  %&& strcmp(par,'lat_raw_sac_rea')
                for t_idx=1:numel(temp_rea)
                    c_trial(t_idx)=text(temp_rea(t_idx),temp_sac(t_idx), [' ', sprintf('%.0f',info_trial(t_idx)), 'r', sprintf('%.0f',info_run(t_idx)) ]);
                end
                set(c_trial,'FontSize', 10)
            end
            
            c_T=text(max(temp_rea),max(temp_sac), ['r', sprintf('%.2f',c_R), ' p', sprintf('%.2f',c_P), con]);
            %             monkeypsych_analyze_working({'K:\Data\Linus_ina\20160805',2},{'summary',0,'trial_set',[880]})
            %             c_T=text(id_x,id_y, ['r', sprintf('%.2f',c_R), ' p', sprintf('%.2f',c_P),]);
            %             set(gca,'xlim',[-1,1],'ylim',[-1,1])
            set(c_T,'Color', col1(g,:), 'FontSize', 10, 'FontWeight', 'bold')
        end
        
        
    end
    
end
lsline;

end


function [idx, e_bar, Abort_codes]=errors_internal(Abort_codes,idx,groups,side,counter_side,type_effector,hand,Plot_settings,col1,par,stat,sig_residuals)
global GLO
if GLO.instructed_only
    conditions={'IN'};
elseif GLO.choice_only
    conditions={'CH'};
else
    conditions={'IN','CH'};
end
e_bar(idx+1:idx+2*size(Abort_codes.reduced_labels,1))=cell(1,2*size(Abort_codes.reduced_labels,1));

idx_reference=idx-2;
for g=1:size(groups,1)
    current_group=groups(g,:);
    for d=1:numel(conditions)
        idx=idx_reference+g;
        con=conditions{d};
        col2=Plot_settings.colors.(con);
        condition=[type_effector  '_' side '_' con hand];
        counter_condition=[type_effector  '_' counter_side '_' con hand];
        if ~isfield(current_group,condition)
            idx=idx+2*size(Abort_codes.reduced_labels,1);
            continue;
        end
        to_plot=[current_group.(condition)]';
        if all(isempty(to_plot)) && all(isnan(to_plot))
            idx=idx+2*size(Abort_codes.reduced_labels,1);
            continue;
        else
            isdata=true;
        end
        
        h_to_plot=hist(to_plot,1:1:size(Abort_codes.reduced_labels,1));
        h_to_plot=h_to_plot/sum(h_to_plot);
        temp_sem = 0;
        for k=1:numel(h_to_plot)
            temp_mean = h_to_plot(k);
            idx=idx+2;
            e_bar{idx}.(con) = plot(idx,temp_mean,'o','MarkerFaceColor',col2,'MarkerEdgeColor',col1(g,:),'Linewidth',GLO.linewidth,'MarkerSize',Plot_settings.markersize);
            if GLO.text_in_plot
                t1=text(idx+(0.1), temp_mean+0.05, sprintf('%.2f',temp_mean));
                set(t1,'Color', [0.3 0.3 0.3], 'FontSize', GLO.fontsize_small, 'FontWeight', 'bold','Rotation',90)
            end
            
            if g==1
                t2=text(idx+(0.1), -0.3, Abort_codes.reduced_labels{k});
                set(t2,'Color', [0.3 0.3 0.3], 'FontSize', GLO.fontsize_small-1, 'FontWeight', 'bold','Rotation',90)
            end
            
        end
        
        
    end
    
end
end

function plot_raw_internal(groups,type_effector,hand,Plot_settings,col1,sac_rea)
global GLO
if GLO.instructed_only
    conditions={'IN'};
    lin=['-'];
    col2=[0.5 0.5 0.5];
elseif GLO.choice_only
    conditions={'CH'};
    lin=[':'];
    col2=[0 0 0];
else
    conditions={'IN','CH'};
    lin=['-',':'];
    col2=[0.5 0.5 0.5; 0 0 0];
end

% switch str2num(type_effector(3))
%     case 1
%         mov_in_state = [2,3];
%     case 2
%         mov_in_state = [3,4,5];
%     case 3
%         mov_in_state = [3,6,7,9,10,4,5];
%     case 4
%         mov_in_state = [3,6,8,4,5];
% end

mov_in_state = GLO.state_raw_traces;

for g=1:size(groups,1)
    current_group=groups(g,:);
    %     id_x=id_x+0.1;
    for d=1:numel(conditions)
        con=conditions{d};
        condition=[type_effector  '_' con hand];
        hold on
        if ~isfield(current_group.x,condition)
            continue;
        end
        for jj=1:numel(current_group.x.(condition))
            if ~iscell(current_group.states.(condition)(jj))
                continue
            end
            %               x_temp = {};
            %                 y_temp =  {};
            %                 n_bin = [];
            idx_in_state = ismember(current_group.states.(condition){jj},mov_in_state);
            fix_pos=current_group.fix_pos.(condition)(jj);
            %             if strcmp(sac_rea,'reaches')
            % %                 plot(current_group.x.(condition){jj}(idx_in_state)-real(fix_pos),current_group.y.(condition){jj}(idx_in_state)-imag(fix_pos),'Color',col1(g,:),'LineStyle',lin(d),'LineWidth',GLO.linewidth/2)
            %            continue
            %             else
            %             if GLO.trial_by_trial_traces
            plot(current_group.x.(condition){jj}(idx_in_state),current_group.y.(condition){jj}(idx_in_state),'Color',col1(g,:),'LineStyle',lin(d),'LineWidth',GLO.linewidth/2)
            %             else
            %                 x_temp{jj} =  current_group.x.(condition){jj}(idx_in_state);
            %                 y_temp{jj} =  current_group.y.(condition){jj}(idx_in_state);
            %                 n_bin(jj) = length(x_temp{jj});
            %             end
            
        end
        %         if ~GLO.trial_by_trial_traces
        %            n_bin_temp = min(n_bin);
        %            if n_bin_temp == 0
        %                continue
        %            else
        %             disp('ok')
        %            end
        %         end
    end
end

end

function plot_raw_internal_xy(groups,type_effector,hand,Plot_settings,col1,current_axis,sac_rea)
global GLO
if GLO.instructed_only
    conditions={'IN'};
    lin=['-'];
    col2=[0.5 0.5 0.5];
elseif GLO.choice_only
    conditions={'CH'};
    lin=[':'];
    col2=[0 0 0];
else
    conditions={'IN','CH'};
    lin=['-',':'];
    col2=[0.5 0.5 0.5; 0 0 0];
end

switch str2num(type_effector(3))
    case 1
        mov_in_state = [2,3];
    case 2
        mov_in_state = [4,5];
    case 3
        mov_in_state = [9,10,4,5];
    case 4
        mov_in_state = [4,5];
end

% STATE.INI_TRI = 1; % initialize trial
% STATE.FIX_ACQ = 2; % fixation acquisition
% STATE.FIX_HOL = 3; % fixation hold
% STATE.TAR_ACQ = 4; % target acquisition
% STATE.TAR_HOL = 5; % target hold
% STATE.CUE_ON  = 6; % cue on
% STATE.MEM_PER = 7; % memory period
% STATE.DEL_PER = 8; % delay period
% STATE.TAR_ACQ_INV = 9; % target acquisition invisible
% STATE.TAR_HOL_INV = 10; % target hold invisible
% STATE.MAT_ACQ = 11; % target acquisition in sample to match
% STATE.MAT_HOL = 12; % target acquisition in sample to match
% STATE.MAT_ACQ_MSK = 13; % target acquisition in sample to match
% STATE.MAT_HOL_MSK = 14; % target acquisition in sample to match
% STATE.SEN_RET     = 15; % return to sensors for poffenberger
% STATE.ABORT     = 19;
% STATE.SUCCESS   = 20;
% STATE.REWARD    = 21;
% STATE.ITI       = 50;
% STATE.CLOSE     = 99;


for g=1:size(groups,1)
    current_group=groups(g,:);
    %     id_x=id_x+0.1;
    for d=1:numel(conditions)
        con=conditions{d};
        condition=[type_effector  '_' con hand];
        hold on
        if ~isfield(current_group.x,condition)
            continue;
        end
        for jj=1:numel(current_group.x.(condition))
            if ~iscell(current_group.states.(condition)(jj))
                continue
            end
            idx_in_state = ismember(current_group.states.(condition){jj},mov_in_state);
            fix_pos=current_group.fix_pos.(condition)(jj);
            samples_in_states = find(idx_in_state==1);
            
            if strcmp(current_axis,'x')
                if ~isempty(current_group.time.(condition){jj}(idx_in_state))
                    plot(current_group.time.(condition){jj}(idx_in_state)-current_group.time.(condition){jj}(samples_in_states(1)),current_group.x.(condition){jj}(idx_in_state)-real(fix_pos),'Color',col1(g,:),'LineStyle',lin(d),'LineWidth',GLO.linewidth/2)
                end
            else
                if ~isempty(current_group.time.(condition){jj}(idx_in_state))
                    plot(current_group.time.(condition){jj}(idx_in_state)-current_group.time.(condition){jj}(samples_in_states(1)),current_group.y.(condition){jj}(idx_in_state)-imag(fix_pos),'Color',col1(g,:),'LineStyle',lin(d),'LineWidth',GLO.linewidth/2)
                end
            end
            if ((~isempty(current_group.lat.(condition)) && ~all(isnan(current_group.lat.(condition)))) && strcmp(sac_rea,'saccades')) && GLO.saccade_in_raw
                for n_sac=1:numel(current_group.lat.(condition))
                    line([current_group.lat.(condition)(n_sac) current_group.lat.(condition)(n_sac)],[-30 30],'Color',col1(g,:),'LineStyle',lin(d),'LineWidth',GLO.linewidth/4)
                end
            end
            %             plot(current_group.x.(condition){jj}(idx_in_state)-real(fix_pos),current_group.time.(condition){jj}(idx_in_state),'Color',col1(g,:),'LineStyle',lin(d),'LineWidth',GLO.linewidth/2)
            
            %             plot(current_group.x.(condition){jj}(idx_in_state)-real(fix_pos),current_group.y.(condition){jj}(idx_in_state)-imag(fix_pos),'Color',col1(g,:),'LineStyle',lin(d),'LineWidth',GLO.linewidth)
        end
    end
end

end

%% STATISTICS FUNCTIONS

function [par_label_IN, par_label_CH]=stat_sig(idx,sac_rea,groups,side,counter_side,type_effector,hand,Plot_settings,col1,par,stat,precision)
global GLO
par_label_IN = {NaN;NaN};
par_label_CH = {NaN;NaN};
if GLO.instructed_only
    conditions={'IN'};
elseif GLO.choice_only
    conditions={'CH'};
else
    conditions={'IN','CH'};
end
for g=1:size(groups,1)
    current_group=groups(g,:);
    for p=1:numel(current_group)
        idx=idx+1;
        for d=1:numel(conditions)
            con=conditions{d};
            condition=[type_effector  '_' side '_' con hand];
            if ~isfield(current_group(p),condition) %|| ~isfield(current_group(p),counter_condition)
                continue;
            end
            if GLO.trial_by_trial %|| strcmp('successful',par)
                
                if strcmp(con,'IN')
                    par_label_IN{1,:} = {condition};
                    par_label_IN{2,:} = stat.groups.(sac_rea).(par).(condition).raw_of_raw{4};
                    if precision
                        par_label_IN{1,:} = {condition};
                        par_label_IN{2,:} = stat.groups.(sac_rea).(par).(condition).raw_of_std{4};
                    end
                else
                    par_label_CH{1,:} = {condition};
                    par_label_CH{2,:} = stat.groups.(sac_rea).(par).(condition).raw_of_raw{4};
                    if precision
                        par_label_CH{1,:} = {condition};
                        par_label_CH{2,:} = stat.groups.(sac_rea).(par).(condition).raw_of_std{4};
                    end
                end
            else
                if strcmp(con,'IN')
                    par_label_IN{1,:} = {condition};
                    par_label_IN{2,:} = stat.groups.(sac_rea).(par).(condition).raw_of_mean{4};
                    if precision
                        par_label_IN{1,:} = {condition};
                        par_label_IN{2,:} = stat.groups.(sac_rea).(par).(condition).raw_of_std{4};
                    end
                else
                    par_label_CH{1,:} = {condition};
                    par_label_CH{2,:} = stat.groups.(sac_rea).(par).(condition).raw_of_mean{4};
                    if precision
                        par_label_CH{1,:} = {condition};
                        par_label_CH{2,:} = stat.groups.(sac_rea).(par).(condition).raw_of_std{4};
                    end
                end
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

end

%% GENERAL FUNCTIONS

function title_and_save(figure_handle,plot_title)
global GLO
mtit(figure_handle,  plot_title, 'xoff', -0.0, 'yoff', 0.04, 'color', [0 0 0], 'fontsize', GLO.fontsize_titles_big,'Interpreter', 'none');
stampit;
% saveas(gcf,[GLO.folder_to_save plot_title print_out]);
if GLO.create_pdf
    %     wanted_size=[50 30];
    %     set(figure_handle, 'Paperunits','centimeters','PaperSize', wanted_size,'PaperPositionMode', 'manual','PaperPosition', [0 0 wanted_size]);
    set(figure_handle, 'PaperPositionMode', 'auto', 'PaperOrientation', 'portrait');
    switch GLO.append_pdfs
        case 0
            export_fig([GLO.folder_to_save filesep plot_title ], '-pdf','-transparent') % pdf by run
            saveas(gcf,[GLO.folder_to_save filesep plot_title], 'fig')
        case 1
            export_fig([GLO.folder_to_save filesep plot_title, 'appended batches'], '-pdf', '-append','-transparent') % pdf by run
            saveas(gcf,[GLO.folder_to_save filesep plot_title], 'fig')
    end
    close all
end
end

function plot_file_names(batch,summary_figure)

filetable.Control               = vertcat(batch.files_for_input.Control{:}{:});
for idx_batch                   = 1:numel(batch.files_for_input.Control{:})
    id_b.c(idx_batch)           = size(batch.files_for_input.Control{:}{idx_batch},1);
end

background.c=[];
batch_n.c=[];
for idd = 1:numel(id_b.c)
    pair = mod(idd, 2);
    if pair == 0
        background.c = [background.c; repmat([0.5 0.5 0.5],id_b.c(idd),1)];
    else
        background.c = [background.c; repmat([0.7 0.7 0.7],id_b.c(idd),1)];
    end
    batch_n.c= [batch_n.c; repmat(idd,id_b.c(idd),1)];
end
batch_n.c_cell = num2cell(batch_n.c);

if isfield(filetable,'Experimental')
    filetable.Experimental          = vertcat(batch.files_for_input.Experimental{:}{:});
    for idx_batch                   = 1:numel(batch.files_for_input.Experimental{:})
        id_b.e(idx_batch)           = size(batch.files_for_input.Experimental{:}{idx_batch},1);
    end
    
    
    background.e=[];
    batch_n.e=[];
    for idd = 1:numel(id_b.e)
        pair = mod(idd, 2);
        if pair == 0
            background.e = [background.e; repmat([0.5 0.5 0.5],id_b.e(idd),1)];
        else
            background.e = [background.e; repmat([0.7 0.7 0.7],id_b.e(idd),1)];
        end
        batch_n.e= [batch_n.e; repmat(idd,id_b.e(idd),1)];
    end
    batch_n.e_cell = num2cell(batch_n.e);
end

% plot_title                      = [' Summary 6, ' 'filelist'];
% summary_figure                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
% summary_figure1                  = figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);



t1 = uitable('Parent', summary_figure, 'Data', [batch_n.c_cell, filetable.Control],'ColumnName',{'Batch', 'Control sessions','Control runs'},'ColumnWidth',{50 400 50},'RowStriping','off','BackgroundColor', [0.9 0.9 0.9]);
pos = get(subplot(1,2,1),'position');
delete(subplot(1,2,1))
set(t1,'units','normalized')
set(t1,'position',pos)

if isfield(filetable,'Experimental')
    t2 = uitable('Parent', summary_figure, 'Data', [batch_n.e_cell, filetable.Experimental], 'ColumnName', {'Batch', 'Experimental sessions', 'Experimental runs'},'ColumnWidth',{50 400 50},'RowStriping','off','BackgroundColor', [0.8 0.8 0.8]);
    pos = get(subplot(1,2,2),'position');
    delete(subplot(1,2,2))
    set(t2,'units','normalized')
    set(t2,'position',pos)
end
end

function [G, Abort_codes] = error_str_2_num(input_group1, input_group2)
% courtesy of DA
% only input is the OUT_COMP from monkypsychAnalyze code
% first modified: 20140813
% Last modified: 20150430

Abort_codes.original= {
    'NO ABORT', 1;
    'ABORT_USE_INCORRECT_HAND' ,2;
    'ABORT_HND_FIX_ACQ_STATE'  ,3;
    'ABORT_HND_FIX_HOLD_STATE' ,4;
    'ABORT_HND_TAR_ACQ_STATE'  ,5;
    'ABORT_HND_TAR_HOLD_STATE' ,6;
    'ABORT_RELEASE_SENSOR_HOLD_STATE' ,7;
    'ABORT_HND_TAR_ACQ_INV_STATE'  ,8;
    'ABORT_HND_TAR_HOLD_INV_STATE' ,9;
    'ABORT_HND_CUE_ON_STATE' ,10;
    'ABORT_HND_MEM_PER_STATE',11;
    'ABORT_HND_DEL_PER_STATE',12;
    'ABORT_EYE_FIX_ACQ_STATE',13;
    'ABORT_EYE_FIX_HOLD_STATE',14;
    'ABORT_EYE_TAR_ACQ_STATE',15;
    'ABORT_EYE_TAR_HOLD_STATE',16;
    'ABORT_EYE_TAR_ACQ_INV_STATE',17;
    'ABORT_EYE_TAR_HOLD_INV_STATE',18;
    'ABORT_EYE_CUE_ON_STATE',19;
    'ABORT_EYE_MEM_PER_STATE',20;
    'ABORT_EYE_DEL_PER_STATE',21;
    'ABORT_JAW',22;
    'ABORT_BODY',23;
    'UNKNOWN ERROR CODE',24;
    'ABORT_DIRTY_SENSORS', 25;
    };


Abort_codes.labels= {
    'Success', 1;
    'Wrong hnd' ,2;
    'Hnd Fix Acq'  ,3;
    'Hnd Fix Hol' ,4;
    'Hnd Tar Acq'  ,5;
    'Hnd Tar Hol' ,6;
    'Rel Sen Hol' ,7;
    'Hnd Tar Acq Inv'  ,8;
    'Hnd Tar Hol Inv' ,9;
    'Hnd Cue On' ,10;
    'Hnd Mem',11;
    'Hnd Del',12;
    'Eye Fix Acq',13;
    'Eye Fix Hol',14;
    'Eye Tar Acq',15;
    'Eye Tar Hol',16;
    'Eye Tar Acq Inv',17;
    'Eye Tar Hol Inv',18;
    'Eye Cue On',19;
    'Eye Mem',20;
    'Eye Del',21;
    'Jaw moved',22;
    'Body moved',23;
    'Unknown',24;
    'Dirty Sensors', 25;
    };

uniques=[];
groups=[input_group1; input_group2];
for gr= 1:numel(groups)
    fns = fieldnames(groups(gr));
    for fn_idx=1:numel(fns)
        current_group=fns{fn_idx};
        if ~iscell(groups(gr).(current_group)) && isnan(groups(gr).(current_group)), groups(gr).(current_group)={}; end
        [~, G(gr).(current_group)] = ismember([groups(gr).(current_group)]',Abort_codes.original(:,1));
        %unique(ismember([groups(gr).(current_group)]',Abort_codes.original(:,1)))
        [uniques] = [uniques; unique(G(gr).(current_group))];
    end
end
uniques = unique(uniques);

Abort_codes.reduced_original = Abort_codes.original(uniques,1);
Abort_codes.reduced_labels = Abort_codes.labels(uniques,1);

% [~,id] = ismember(G(gr).(current_group),uniques);

for gr= 1:numel(groups)
    fns = fieldnames(groups(gr));
    for fn_idx=1:numel(fns)
        current_group=fns{fn_idx};
        [~, G(gr).(current_group)] = ismember(G(gr).(current_group),uniques);
        unique(G(gr).(current_group));
    end
end

a=1;
end

