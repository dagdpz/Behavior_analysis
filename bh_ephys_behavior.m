function bh_ephys_behavior(project,version,dataset,epoch,basefolder)
% project='Pulv_eye_gaze_position';
% version='paper';
% dataset='_dPulv_PT0_Msac_mov';
project='Pulv_oculomotor';
version='paper';
dataset='_dPulv_PT0_Msac_opt';
epoch={'Thol',5,	0.2,    0.5};
basefolder='C:\Users\lschneider\Desktop';
mainpath=[basefolder filesep project dataset];
load(['Y:\Projects\' project '\ephys\' version '\behaviour_filelist.mat']);

keys.cal.MA_selection                   ={'correct_offset',0,'success',1,'display',0,'keep_raw_data',1,'saccade_definition',4,'reach_1st_pos',1,'correlation_conditions',{}};                        % if you want to run MA with specific settings
keys.colors.fix_offset          =[236 32 38; 16 159 218; 247 148 36]/255;
recalibrated=1;


if recalibrated==1
    keys.path_to_save=['Y:\Projects\' project '\behavior\correctiverecalibrated'];
else
    keys.path_to_save=['Y:\Projects\' project '\behavior\correctivenotrecalibrated'];
end
    keys.plot.export=1;


%% modify filelist to create batches for each session
for monkeys={'Cur','Lin','Both'}
    monkey=monkeys{:};
    if strcmp(monkey,'Both')
        monkey='Lin';
        taskfield=[monkey dataset];
        mainfolder=filelist.(taskfield){1};
        dashstr=strfind(mainfolder,'\');
        mainfolder=mainfolder(1:dashstr(end));
        filelist_in=filelist_formatted.(taskfield);
        filelist_in(:,1)=cellstr(num2str(vertcat(filelist_in{:,1})));
        filelist_in(:,1)=strcat(mainfolder, filelist_in(:,1));
        
        if recalibrated==1
        aa=dir([mainpath filesep monkey filesep]);
        sessions={aa([aa.isdir]).name};
        sessions=sessions(3:end);
        for s=1:numel(sessions)
            filelist_in{s,1}=[mainpath filesep monkey filesep sessions{s}];
        end
        end
        
        monkey='Cur';
        taskfield=[monkey dataset];
        mainfolder=filelist.(taskfield){1};
        dashstr=strfind(mainfolder,'\');
        mainfolder=mainfolder(1:dashstr(end));
        filelist_in2=filelist_formatted.(taskfield);
        filelist_in2(:,1)=cellstr(num2str(vertcat(filelist_in2{:,1})));
        filelist_in2(:,1)=strcat(mainfolder, filelist_in2(:,1));
        
        
        if recalibrated==1
        aa=dir([mainpath filesep monkey filesep]);
        sessions={aa([aa.isdir]).name};
        sessions=sessions(3:end);
        for s=1:numel(sessions)
            filelist_in2{s,1}=[mainpath filesep monkey filesep sessions{s}];
        end
        end
        
        filelist_in=vertcat(filelist_in, filelist_in2);
        monkey='Both';
    else
        taskfield=[monkey dataset];
        mainfolder=filelist.(taskfield){1};
        dashstr=strfind(mainfolder,'\');
        mainfolder=mainfolder(1:dashstr(end));
        filelist_in=filelist_formatted.(taskfield);
        filelist_in(:,1)=cellstr(num2str(vertcat(filelist_in{:,1})));
        filelist_in(:,1)=strcat(mainfolder, filelist_in(:,1));
        
        
        if recalibrated==1
        aa=dir([mainpath filesep monkey filesep]);
        sessions={aa([aa.isdir]).name};
        sessions=sessions(3:end);
        for s=1:numel(sessions)
            filelist_in{s,1}=[mainpath filesep monkey filesep sessions{s}];
        end
        end
    end
    MA_input={};
    for c=1:size(filelist_in,1)
        MA_input=[MA_input, {{filelist_in{c,1}, filelist_in{c,2}}},{keys.cal.MA_selection}];
    end
    data=monkeypsych_analyze_working(MA_input{:});
    
    all_saccades=[data{:}];
    all_saccades=vertcat(all_saccades.saccades);
    positions=unique(round([all_saccades.tar_pos] - [all_saccades.fix_pos]));
    positions(isnan(positions))=[];
    fixations=unique(real([all_saccades.fix_pos])); % will need to use precision here
    distances=abs(positions);
    unique_distances=unique(distances);
    
    clear per_pos mean_RTs mean_dur mean_amp mean_vel
    for s=1:numel(data)
        saccades=[data{s}.saccades];
        states=[data{s}.states];
        raw=[data{s}.raw];
        mean_RTs(s)=nanmean([saccades.lat])*1000;
        mean_dur(s)=nanmean([saccades.dur])*1000;
        mean_amp(s)=nanmean(abs([saccades.endpos]-[saccades.startpos]));
        mean_vel(s)=nanmean([saccades.velocity]);
        mean_prec(s)=nanmean([saccades.endpos]-[saccades.fix_pos]);
        
        fix_pos=[saccades.fix_pos];
        for f=1:numel(fixations)
            
            for p=1:numel(positions)
                tr_idx=(round([saccades.tar_pos] - fix_pos))==positions(p) & real([saccades.fix_pos])==fixations(f);
                
                
                tr_idx_int=find(tr_idx);
                av_pos=[];
                acq_dur=[];
                for t=1:numel(tr_idx_int)
                    idx=tr_idx_int(t);
                    on=states(idx).MP_states_onset([states(idx).MP_states]==epoch{2});
                    time=raw(idx).time_axis >= on+epoch{3}  & raw(idx).time_axis <= on+epoch{4}; %>=?
                    av_pos(t)=median(raw(idx).x_eye(time))+1i*median(raw(idx).y_eye(time));
                    
                    acq_dur(t)=states(idx).MP_states_onset([states(idx).MP_states]==10)-...
                        states(idx).MP_states_onset([states(idx).MP_states]==9);
                end
                
                per_pos(f,p).sac_end_all{s}=[saccades(tr_idx).endpos]-fix_pos(tr_idx)+fixations(f);
                per_pos(f,p).eye_in_thol_all{s}=av_pos-fix_pos(tr_idx)+fixations(f);
                
                per_pos(f,p).sac_off(s)=nanmean([saccades(tr_idx).accuracy_xy]);
                per_pos(f,p).RTs(s)=nanmean([saccades(tr_idx).lat])*1000;
                per_pos(f,p).dur(s)=nanmean([saccades(tr_idx).dur])*1000;
                per_pos(f,p).vel(s)=nanmean([saccades(tr_idx).velocity]);
                per_pos(f,p).acq_dur(s)=nanmean(acq_dur-[saccades(tr_idx).lat])*1000;
                per_pos(f,p).tar_rad(s)=nanmean([saccades(tr_idx).tar_rad]);
                per_pos(f,p).sac_end(s)=nanmean([saccades(tr_idx).endpos]-fix_pos(tr_idx)+fixations(f));
                per_pos(f,p).sac_prec(s)=nanstd(real([saccades(tr_idx).endpos]-fix_pos(tr_idx)));
                per_pos(f,p).eye_in_thol(s)=nanmean(av_pos-fix_pos(tr_idx)+fixations(f));
                per_pos(f,p).amp(s)=nanmean(abs([saccades(tr_idx).endpos]-[saccades(tr_idx).startpos]));
                per_pos(f,p).finalamp(s)=nanmean(abs(av_pos-fix_pos(tr_idx)));
                per_pos(f,p).final_prec(s)=nanstd(av_pos-fix_pos(tr_idx));
                
                
                per_pos(f,p).meanRTs= nanmean([per_pos(f,p).RTs]);
                per_pos(f,p).meandur= nanmean([per_pos(f,p).dur]);
                per_pos(f,p).meanvel= nanmean([per_pos(f,p).vel]);
                per_pos(f,p).meanamp= nanmean([per_pos(f,p).amp]);
                per_pos(f,p).meanacq_dur= nanmean([per_pos(f,p).acq_dur]);
                per_pos(f,p).meanfinalamp= nanmean([per_pos(f,p).finalamp]);
                
                per_pos(f,p).fixation(s)=f;
                per_pos(f,p).distance(s)=abs(positions(p));
                per_pos(f,p).position(s)=p;
            end
        end
    end
    
    %% plotting
    
    
    %% accuracy
    plot_title=[monkey ' accuracy'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    Colors=keys.colors.fix_offset;
    
    angles=0:pi/180:2*pi;
    
    hold on
    for f=1:numel(fixations)
        
        for p=1:numel(positions)
            
            r=5;
            x0=positions(p)+fixations(f);
            scatter(real(x0),imag(x0),'+','markeredgecolor','k');
            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
            
            plot([(nanmean([per_pos(f,p).sac_end])) , nanmean([per_pos(f,p).eye_in_thol])],'color',Colors(f,:));
            scatter(real(nanmean([per_pos(f,p).sac_end])),imag(nanmean([per_pos(f,p).sac_end])),30,Colors(f,:),'filled');
            %nanmean(per_pos(f,p).dur);
        end
    end
    axis equal
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    
    %% accuracy per session
    plot_title=[monkey ' accuracy per session'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    Colors=keys.colors.fix_offset;
    
    angles=0:pi/180:2*pi;
    
    hold on
    for f=1:numel(fixations)
        
        for p=1:numel(positions)
            
            r=5;
            x0=positions(p)+fixations(f);
            scatter(real(x0),imag(x0),'+','markeredgecolor','k');
            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
            
            %plot([(nanmean(per_pos(f,p).sac_end)) , nanmean(per_pos(f,p).eye_in_thol)],'color',Colors(f,:));
            scatter(real([per_pos(f,p).sac_end]),imag([per_pos(f,p).sac_end]),30,Colors(f,:),'filled');
            nanmean(per_pos(f,p).dur);
        end
    end
    axis equal
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% accuracy per session
    plot_title=[monkey ' final accuracy per session'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    Colors=keys.colors.fix_offset;
    
    angles=0:pi/180:2*pi;
    
    hold on
    for f=1:numel(fixations)
        
        for p=1:numel(positions)
            
            r=5;
            x0=positions(p)+fixations(f);
            scatter(real(x0),imag(x0),'+','markeredgecolor','k');
            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
            
            %plot([(nanmean(per_pos(f,p).sac_end)) , nanmean(per_pos(f,p).eye_in_thol)],'color',Colors(f,:));
            plot([[per_pos(f,p).sac_end]' , [per_pos(f,p).eye_in_thol]']','color',Colors(f,:));
            scatter(real([per_pos(f,p).eye_in_thol]),imag([per_pos(f,p).eye_in_thol]),15,Colors(f,:),'filled');
            %nanmean(per_pos(f,p).dur);
        end
    end
    axis equal
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    
    %% accuracy per session
    Colors=keys.colors.fix_offset;
    
    angles=0:pi/180:2*pi;
    
    for s=1:numel(data)
        plot_title=[monkey ' final accuracy, session ' num2str(s)];
        filename=plot_title;
        figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
        hold on
        for f=1:numel(fixations)
            
            for p=1:numel(positions)
                
                r=5;
                x0=positions(p)+fixations(f);
                scatter(real(x0),imag(x0),'+','markeredgecolor','k');
                plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
                
                %plot([(nanmean(per_pos(f,p).sac_end)) , nanmean(per_pos(f,p).eye_in_thol)],'color',Colors(f,:));
                plot([[per_pos(f,p).sac_end_all{s}]' , [per_pos(f,p).eye_in_thol_all{s}]']','color',Colors(f,:));
                scatter(real([per_pos(f,p).eye_in_thol_all{s}]),imag([per_pos(f,p).eye_in_thol_all{s}]),15,Colors(f,:),'filled');
                %nanmean(per_pos(f,p).dur);
            end
        end
        axis equal
        ph_title_and_save(figure_handle,filename,plot_title,keys)
    end
    
    %% precision per session
    plot_title=[monkey ' precision per session'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    Colors=keys.colors.fix_offset;
    
    angles=0:pi/180:2*pi;
    
    hold on
    for f=1:numel(fixations)
        
        for p=1:numel(positions)
            
            r=5;
            x0=positions(p)+fixations(f);
            scatter(real(x0),imag(x0),'+','markeredgecolor','k');
            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
            
            x0=nanmean([per_pos(f,p).sac_end]);
            r=nanmean([per_pos(f,p).sac_prec]);
            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'color',Colors(f,:));
            drawring([real(x0) imag(x0)],r-sterr([per_pos(f,p).sac_prec]), r+sterr([per_pos(f,p).sac_prec]),angles,Colors(f,:));
            %scatter(real([per_pos(f,p).sac_prec]),imag([per_pos(f,p).sac_prec]),30,Colors(f,:),'filled');
            %nanmean(per_pos(f,p).dur);
        end
    end
    axis equal
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    
    %% precision per session
    plot_title=[monkey ' precision of session means'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    Colors=keys.colors.fix_offset;
    
    angles=0:pi/180:2*pi;
    
    hold on
    for f=1:numel(fixations)
        
        for p=1:numel(positions)
            
            r=5;
            x0=positions(p)+fixations(f);
            scatter(real(x0),imag(x0),'+','markeredgecolor','k');
            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
            
            x0=nanmean([per_pos(f,p).sac_end]);
            r=nanstd([per_pos(f,p).sac_end]);
            plot(real(x0)+r*sin(angles) , imag(x0)+r*cos(angles),'color',Colors(f,:));
            %scatter(real([per_pos(f,p).sac_prec]),imag([per_pos(f,p).sac_prec]),30,Colors(f,:),'filled');
            %nanmean(per_pos(f,p).dur);
        end
    end
    axis equal
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% final precision per session
    plot_title=[monkey ' final precision per session'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    Colors=keys.colors.fix_offset;
    
    angles=0:pi/180:2*pi;
    
    hold on
    for f=1:numel(fixations)
        
        for p=1:numel(positions)
            
            r=5;
            x0=positions(p)+fixations(f);
            scatter(real(x0),imag(x0),'+','markeredgecolor','k');
            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
            
            x0=nanmean([per_pos(f,p).eye_in_thol]);
            r=nanmean([per_pos(f,p).final_prec]);
            plot(real(x0)+r*sin(angles) , imag(x0)+r*cos(angles),'color',Colors(f,:));
            drawring([real(x0) imag(x0)],r-sterr([per_pos(f,p).final_prec]), r+sterr([per_pos(f,p).final_prec]),angles,Colors(f,:));
            
            %shadedErrorBar(real(x0)+r*sin(angles), imag(x0)+r*cos(angles), repmat(sterr([per_pos(f,p).final_prec]),size(angles)) ,{'color',Colors(f,:)});
            %scatter(real([per_pos(f,p).sac_prec]),imag([per_pos(f,p).sac_prec]),30,Colors(f,:),'filled');
            %nanmean(per_pos(f,p).dur)
        end
    end
    axis equal
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    
    %% final precision per session
    plot_title=[monkey ' final precision of session means'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    Colors=keys.colors.fix_offset;
    
    angles=0:pi/180:2*pi;
    
    hold on
    for f=1:numel(fixations)
        
        for p=1:numel(positions)
            
            r=5;
            x0=positions(p)+fixations(f);
            scatter(real(x0),imag(x0),'+','markeredgecolor','k');
            plot(real(x0)+r*sin(angles), imag(x0)+r*cos(angles),'k');
            
            x0=nanmean([per_pos(f,p).eye_in_thol]);
            r=nanstd([per_pos(f,p).eye_in_thol]);
            plot(real(x0)+r*sin(angles) , imag(x0)+r*cos(angles),'color',Colors(f,:));
            %scatter(real([per_pos(f,p).sac_prec]),imag([per_pos(f,p).sac_prec]),30,Colors(f,:),'filled');
            %nanmean(per_pos(f,p).dur);
        end
    end
    axis equal
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% duration
    
    plot_title=[monkey ' duration'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    bins=min(floor([per_pos.meandur])):1:max(ceil([per_pos.meandur]));
    clear histo
    for d=1:numel(unique_distances)
        histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meandur],bins);
    end
    bar(bins,histo,'stacked');
    x_lim=get(gca,'xlim');
    y_lim=get(gca,'ylim');
    
    for f=1:numel(fixations)
        for d=1:numel(unique_distances)
%             means(d)=nanmean([per_pos(:,distances==unique_distances(d)).meandur]);
%             stds(d)=nanstd([per_pos(:,distances==unique_distances(d)).meandur]);
%             SEs(d)=sterr([per_pos(:,distances==unique_distances(d)).meandur]);
%             text_to_plot=['M: ' num2str(means(d)) ', std: ' num2str(stds(d)) ', SE: ' num2str(SEs(d))];
%             text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
            persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).dur),1);
            text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
            text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
        end
    end
    
    pd = anovan([per_pos.dur]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
    pp = anovan([per_pos.dur]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
    
    text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
    text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
    
    xlabel('saccade duration [ms]');
    ylabel('N conditions');
    title(['across targets: ' num2str(mean([per_pos(:).meandur])) ' + ' num2str(sterr([per_pos(:).meandur]))...
        'across sessions: ' num2str(mean(mean_dur)) ' + ' num2str(sterr(mean_dur)) ...
        'across targets, end as tihol: ' num2str(mean([per_pos(:).meanacq_dur])) ' + ' num2str(sterr([per_pos(:).meanacq_dur])) ]);
    
    legend(cellstr(num2str(unique_distances')));
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% latencies
    
    plot_title=[monkey ' latency'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    bins=min(floor([per_pos.meanRTs])):1:max(ceil([per_pos.meanRTs]));
    clear histo
    for d=1:numel(unique_distances)
        histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanRTs],bins);
    end
    bar(bins,histo,'stacked');
    x_lim=get(gca,'xlim');
    y_lim=get(gca,'ylim');
%     
%     for d=1:numel(unique_distances)
%         means(d)=nanmean([per_pos(:,distances==unique_distances(d)).meanRTs]);
%         stds(d)=nanstd([per_pos(:,distances==unique_distances(d)).meanRTs]);
%         SEs(d)=sterr([per_pos(:,distances==unique_distances(d)).meanRTs]);
%         text_to_plot=['M: ' num2str(means(d)) ', std: ' num2str(stds(d)) ', SE: ' num2str(SEs(d))];
%         text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
%     end

    for f=1:numel(fixations)
        for d=1:numel(unique_distances)
            persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).RTs),1);
            text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
            text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
        end
    end
    
    pd = anovan([per_pos.RTs]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
    pp = anovan([per_pos.RTs]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
    
    text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
    text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
    
    xlabel('saccade reaction time [ms]');
    ylabel('N conditions');
    title(['across targets: ' num2str(mean([per_pos(:).meanRTs])) ' + ' num2str(sterr([per_pos(:).meanRTs]))...
        'across sessions: ' num2str(mean(mean_RTs)) ' + ' num2str(sterr(mean_RTs))]);
    
    legend(cellstr(num2str(unique_distances')));
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% amplitude
    
    plot_title=[monkey ' amplitude'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    bins=min(floor([per_pos.meanamp])):1:max(ceil([per_pos.meanamp]));
    clear histo
    for d=1:numel(unique_distances)
        histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanamp],bins);
    end
    bar(bins,histo,'stacked');
    x_lim=get(gca,'xlim');
    y_lim=get(gca,'ylim');
%     
%     for d=1:numel(unique_distances)
%         means(d)=nanmean([per_pos(:,distances==unique_distances(d)).meanamp]);
%         stds(d)=nanstd([per_pos(:,distances==unique_distances(d)).meanamp]);
%         SEs(d)=sterr([per_pos(:,distances==unique_distances(d)).meanamp]);
%         text_to_plot=['M: ' num2str(means(d)) ', std: ' num2str(stds(d)) ', SE: ' num2str(SEs(d))];
%         text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
%     end
    
    
    for f=1:numel(fixations)
        for d=1:numel(unique_distances)
            persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).amp),1);
            text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
            text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
        end
    end
    
    pd = anovan([per_pos.amp]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
    pp = anovan([per_pos.amp]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
    
    text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
    text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
    
    xlabel('saccade amplitude [�]');
    ylabel('N conditions');
    title(['across targets: ' num2str(mean([per_pos(:).meanamp])) ' + ' num2str(sterr([per_pos(:).meanamp]))...
        'across sessions: ' num2str(mean(mean_amp)) ' + ' num2str(sterr(mean_amp))  ]);
    
    legend(cellstr(num2str(unique_distances')));
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    %% final amplitude
    
    plot_title=[monkey ' final amplitude'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    bins=min(floor([per_pos.meanfinalamp])):1:max(ceil([per_pos.meanfinalamp]));
    clear histo
    for d=1:numel(unique_distances)
        histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanfinalamp],bins);
    end
    bar(bins,histo,'stacked');
    x_lim=get(gca,'xlim');
    y_lim=get(gca,'ylim');
%     
%     for d=1:numel(unique_distances)
%         means(d)=nanmean([per_pos(:,distances==unique_distances(d)).meanfinalamp]);
%         stds(d)=nanstd([per_pos(:,distances==unique_distances(d)).meanfinalamp]);
%         SEs(d)=sterr([per_pos(:,distances==unique_distances(d)).meanfinalamp]);
%         text_to_plot=['M: ' num2str(means(d)) ', std: ' num2str(stds(d)) ', SE: ' num2str(SEs(d))];
%         text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
%     end
%     
    
    
    for f=1:numel(fixations)
        for d=1:numel(unique_distances)
            persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).finalamp),1);
            text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
            text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
        end
    end
    
    pd = anovan([per_pos.finalamp]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
    pp = anovan([per_pos.finalamp]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
    
    text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
    text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
    
    xlabel('final distance [�]');
    ylabel('N conditions');
    title(['across targets: ' num2str(mean([per_pos(:).meanfinalamp])) ' + ' num2str(sterr([per_pos(:).meanfinalamp])) ]);
    
    legend(cellstr(num2str(unique_distances')));
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    
    
    %% velocity
    
    plot_title=[monkey ' velocity'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    bins=min(floor([per_pos.meanvel]/10)*10):10:max(ceil([per_pos.meanvel]/10)*10);
    clear histo
    for d=1:numel(unique_distances)
        histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanvel],bins);
    end
    bar(bins,histo,'stacked');
    x_lim=get(gca,'xlim');
    y_lim=get(gca,'ylim');
%     
%     for d=1:numel(unique_distances)
%         means(d)=nanmean([per_pos(:,distances==unique_distances(d)).meanvel]);
%         stds(d)=nanstd([per_pos(:,distances==unique_distances(d)).meanvel]);
%         SEs(d)=sterr([per_pos(:,distances==unique_distances(d)).meanvel]);
%         text_to_plot=['M: ' num2str(means(d)) ', std: ' num2str(stds(d)) ', SE: ' num2str(SEs(d))];
%         text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
%     end
%     
    
    for f=1:numel(fixations)
        for d=1:numel(unique_distances)
            persession=nanmean(vertcat(per_pos(f,distances==unique_distances(d)).vel),1);
            text_to_plot=[num2str(nanmean(persession)) '+' num2str(sterr(persession))];
            text(x_lim(1)+diff(x_lim)/5*f,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
        end
    end
    
    pd = anovan([per_pos.vel]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
    pp = anovan([per_pos.vel]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
    
    text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
    text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
    
    xlabel('velocity [�/s]');
    ylabel('N conditions');
    title(['across targets: ' num2str(mean([per_pos(:).meanvel])) ' + ' num2str(sterr([per_pos(:).meanvel]))...
        'across sessions: ' num2str(mean(mean_vel)) ' + ' num2str(sterr(mean_vel))  ]);
    
    legend(cellstr(num2str(unique_distances')));
    ph_title_and_save(figure_handle,filename,plot_title,keys)
    
    
    %% precision
    
    plot_title=[monkey ' precision'];
    filename=plot_title;
    figure_handle= figure('units','normalized','outerposition',[0 0 1 1],'name',plot_title);
    
    bins=min(floor([per_pos.meanvel]/10)*10):10:max(ceil([per_pos.meanvel]/10)*10);
    clear histo
    for d=1:numel(unique_distances)
        histo(:,d)=hist([per_pos(:,distances==unique_distances(d)).meanvel],bins);
    end
    bar(bins,histo,'stacked');
    x_lim=get(gca,'xlim');
    y_lim=get(gca,'ylim');
    
%     for d=1:numel(unique_distances)
%         means(d)=nanmean([per_pos(:,distances==unique_distances(d)).meanvel]);
%         stds(d)=nanstd([per_pos(:,distances==unique_distances(d)).meanvel]);
%         SEs(d)=sterr([per_pos(:,distances==unique_distances(d)).meanvel]);
%         text_to_plot=['M: ' num2str(means(d)) ', std: ' num2str(stds(d)) ', SE: ' num2str(SEs(d))];
%         text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*d,text_to_plot);
%     end
    
    pd = anovan([per_pos.sac_prec]',[[per_pos.fixation]' [per_pos.distance]'],'model','full','display','off');
    pp = anovan([per_pos.sac_prec]',[[per_pos.fixation]' [per_pos.position]'],'model','full','display','off');
    
    text_to_plot=['p anova FxD: ' num2str(pd(1)) ', ' num2str(pd(2)) ', ' num2str(pd(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+1),text_to_plot);
    text_to_plot=['p anova FxP: ' num2str(pp(1)) ', ' num2str(pp(2)) ', ' num2str(pp(3)) ];
    text(x_lim(1)+diff(x_lim)/2,y_lim(2)-diff(y_lim)/20*(d+2),text_to_plot);
    
    xlabel('precision [�]');
    ylabel('N conditions');
    title(['across targets: ' num2str(mean([per_pos(:).sac_prec])) ' + ' num2str(sterr([per_pos(:).sac_prec]))]);
    
    legend(cellstr(num2str(unique_distances')));
    ph_title_and_save(figure_handle,filename,plot_title,keys)
end
end


function hArrow = drawArrow(p0,p1,color)
% drawArrow(p0,p1)
% Draws a simple arrow in 2D, from p0 to p1.
% from: https://de.mathworks.com/matlabcentral/fileexchange/55181-drawarrow by Matthew Kelly
%
% INPUTS:
%   p0 = [x0; y0] = position of the tail
%   p1 = [x1; y1] = position of the tip
%   color = arrow color. Optional: default is black
%       --> can be 'r','g','b','c','m','y','w', 'k' or a 1x3 color vector
%
% OUTPUTS:
%   hArrow = handle to the patch object representing the arrow
%
% Defaults:
if nargin == 2
    color = 'k';
end
% Parameters:
W1 = 0.08;   % half width of the arrow head, normalized by length of arrow
W2 = 0.014;  % half width of the arrow shaft
L1 = 0.18;   % Length of the arrow head, normalized by length of arrow
L2 = 0.13;  % Length of the arrow inset
% Unpack the tail and tip of the arrow
x0 = p0(1);
y0 = p0(2);
x1 = p1(1);
y1 = p1(2);
% Start by drawing an arrow from 0 to 1 on the x-axis
P = [...
    0, (1-L2), (1-L1), 1, (1-L1), (1-L2), 0;
    W2,    W2,     W1, 0,    -W1,    -W2, -W2];
% Scale,rotate, shift and plot:
dx = x1-x0;
dy = y1-y0;
Length = sqrt(dx*dx + dy*dy);
Angle = atan2(-dy,dx);
P = P.*repmat([Length;1],1,7);   %Scale Length*
P = [cos(Angle), sin(Angle); -sin(Angle), cos(Angle)]*P;  %Rotate
P = p0(:)*ones(1,7) + P;  %Shift
% Plot!
hArrow = patch(P(1,:), P(2,:),color);  axis equal;
set(hArrow,'EdgeColor',color);
end

function drawring(center,rin, rout,angles,color)
xin = center(1) + rin*cos(angles);
xout = center(1) + rout*cos(angles);
yin = center(2) + rin*sin(angles);
yout = center(2) + rout*sin(angles);
% Make patch
hp = patch([xout,xin],[yout,yin],color,'linestyle','none');
end