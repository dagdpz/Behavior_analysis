%% Linus ephys [20160203 20160506] datasets

clc
clear batch,
close all
warning('off','all')
% EH_compare_groups(batch,testing)

global GLO
GLO.accuracy_as_absolute            =   1; %% 1 meaning it computes the averages of x and y first and then creates the euclidean
GLO.delete_last                     =   0; %1 delete the last trial of each run
GLO.fontsize_titles_big             =   16;
GLO.fontsize_small                  =   8;
GLO.fontsize_ticks                  =   12;
GLO.fontsize_labels                 =   12;
GLO.linewidth                       =   2;
GLO.plot_raw_endpoints              =   1; %1 means 1 point per trial, 0 means average across trials
GLO.calculate_statististics         =   1;
GLO.parametric_testing              =   1;
GLO.plot_statististics              =   1;
GLO.plot_it                         =   1;
GLO.create_pdf                      =   1;
GLO.append_pdfs                     =   0;
GLO.parent_folder                   =   '';
GLO.folder_to_save                  =   'Y:\Projects\Simultaneous_dPul_PPC_recordings\behavior\beh_analysis\Linus\LIP_dPul_control';
GLO.type_of_free_gaze               =   '6';
GLO.one_subject                     =   0;
GLO.trial_by_trial                  =   0; % for statistics, 0 means calculate statistics bases on average per run
GLO.CDF                             =   0; % 1 plot cumulative distribution function
GLO.text_in_plot                    =   1; % plot mean and std in text
GLO.same_day                        =   0;
GLO.testing_patient                 =   0;
GLO.instructed_only                 =   0;
GLO.choice_only                     =   0;
GLO.only_significant                =   1; % for sigstar
GLO.only_success_for_accuracy       =   0;
GLO.only_between_group_comparissons =   0;
GLO.point_per_batch                 =   0; %0 average across session , 1 % 1 point per run
GLO.summary                         =   [1 2 5 10 11];  %1 2 5 10 11 which plots (vector of number of -1 for all plots) [1 2 3 5 10]
GLO.target_locations_in_raw         =   0; %not used anymore
GLO.saccade_in_raw                  =   0; %only for one plot, see beh_compare_groups, not used anymore
GLO.modify_positions                =   0; % used in reallocate_positions_from_mpa ?
GLO.euclideans_reach                =   [-15, 15]; % used in reallocate_positions_from_mpa ?
GLO.trial_numbers                   =   0; %in correlation plot
GLO.keep_raw_output                 =   1; % 1 = save raw data in the output structure and plot raw traces in some plots 
GLO.hits_in_plot                    =   1; % plot the number of hits per condition on plots
GLO.min_hits                        =   0; %or 1 for 50 hits min
GLO.only_successful_side_selection  =   1; %0 takes in account aborted trial from state_inf, 1 only successful trials
%next 3 for plotting
GLO.saccades_effectors              = {'0'};
GLO.reaches_effectors               = {'4'};
GLO.types_to_plot                   = {'4'};
GLO.saccades.effectors_raw_xy       = {'0'};
GLO.reaches.effectors_raw_xy        = {'4'};
GLO.state_raw_traces                = [4 5];



% steady.passfilter                   =   {'saccades','lat',0.08, 0.5;'reaches','lat',0.1, 10};

%steady.passfilter                   =   {'saccades','lat',0, 3;'reaches','lat',0.2, 10;}; % add
steady.passfilter                   =   {'saccades','lat',0.01, 3;'reaches','lat',0, 10;}; % add

% steady.passfilter                   =   {'saccades','lat',0.01, 10;'reaches','lat',0.01, 10;}; %% !!!!!! remove  


GLO.remove_outliers                 =   1;
% 1 ERROR BARS
% 2 HISTOGRAMS
% 3 ACCURACY
% 4 CORRELATIONS
% 5 HAND RELEASE
% 6 RAW XY 2D 
% 7X X VS TIME
% 7Y Y VS TIME
% 8 FILELIST
% 9 ERRORS
% 10 CH - IN

% % 1: saccade, that ended up closest to the target and was big enough
% % 2: biggest saccade that ended up close enough
% % 3: last saccade in the state
% % 4: first saccade in the state
% % 5: the first that is bigger than 'sac_min_amp'
steady.reach_1st_pos                =   1; %%new used for all but
steady.reach_1st_pos_in             =   0; %%new
steady.reach_pos_at_state_change    =   0; %%new

% 'nsacc_max'                       maximum number of saccades calculated for one state             integer values >=1
% 'sac_ini_t'                       initiation velocity threshold for saccades (in degrees/sec)     defines saccade detection, be careful
% 'sac_end_t'                       end velocity threshold for saccades (in degrees/sec)            defines saccade detection, be careful
% 'sac_min_dur'                     Minimum saccade duration (in sec) for saccade detection         defines saccade detection, be careful
% 'min_sac_amplitude'               Minimum saccade amplitude (in deg) for saccade detection        defines the minimum amplitude for selected saccades in saccade_definition for closest saccades (1)
% 'max_sac_dist'                    maximum distance from target for saccades                       [max_dist_x max_dist_y] defines if a miss is still taken as intended to go to that target position
%                                                                                                   which is crucial for choice trials to define the desired target location
%                                                                                                   also defines maximum distance to the target for selected saccades in saccade_definition for biggest saccades (2)


steady.remove_outliers              =   GLO.remove_outliers;
steady.sac_ini_t                  = 200;
steady.sac_end_t                  = 50;
steady.sac_min_dur                   = 0.01;
% steady.max_sac_dist               = 7;
steady.eyetracker_sample_rate     = 220; % Hertz
steady.correct_offset             = 1;
steady.saccade_definition         = 1;
steady.smoothing_samples          = 15; %downsampling does a diff to find relevant changing samples, then it inerpolates to 1 ms between such samples and then it smoothens, the smoothing is set such that it grabs one sample before and one after plus a bit more but not reaching 2 samples 4.54 ms each sample of the eyetracker
steady.sac_min_amp                = 2;
steady.keep_raw_data              = 1;
steady.correlation_mode           = 'pearson';
% steady.display                    = 0;
steady.downsampling               = 1;
        steady.tar_range_x                  =   [NaN;NaN];
        steady.tar_range_y                  =   [NaN;NaN];

%  load('Y:\Projects\Simultaneous_dPul_PPC_recordings\ephys\dPul_control_LIP_Lin_8s_\behaviour_filelist.mat');
 
%  filelist_formatted_control=filelist_formatted.Lin_LIP_R_PT0_Dsac_han;
%  filelist_formatted_inactivation=filelist_formatted.Lin_LIP_R_PT1_Dsac_han;

subject_ID{1}='Control';
group{1}                        = repmat({'Linus'},1,8);
dates_subject_in{1}             = {20210623,20210729,20210910,20211013,20211028,20211029,20211203,20211210};
batching{1}.runs                = {2;4;2;2;2;2;3;2};  % either empty or specific runs specified
batching{1}.inactivation_sites   = {'R','R','R','R','R','R','R','R'};
batching_type{1}                = 1; % 1 run by run, 2 session by session, 3 group by group
batching{1}.range_of_dates      = 0;

% subject_ID{2}='Experimental';
% group{2}                        = repmat({'Linus'},size(filelist_formatted_inactivation,1),1);
% dates_subject_in{2}                = filelist_formatted_inactivation(:,1);
% batching{2}.runs                = filelist_formatted_inactivation(:,2);  % either empty or specific runs specified
% batching{2}.inactivation_sites   = {'R','R','R','R','R','R','R','R','R','R','R','R','R','L','L','L','L','L','L'};
% batching_type{2}                = 2; % 1 run by run, 2 session by session, 3 group by group
% batching{2}.range_of_dates      = 0;% 

subject_ID{2}='Experimental';
group{2}                        = repmat({'Linus'},1,8);
dates_subject_in{2}             = {20210623,20210729,20210910,20211013,20211028,20211029,20211203,20211210};
batching{2}.runs                = {3;5;3;3;3;3;4;3};  % either empty or specific runs specified
batching{2}.inactivation_sites   = {'R','R','R','R','R','R','R','R'};
batching_type{2}                = 1; % 1 run by run, 2 session by session, 3 group by group
batching{2}.range_of_dates      = 0;% 

% subject_ID{2}='Experimental';
% group{2}                        = repmat({'Bacchus'},1,2);
% dates_subject_in{2}             = repmat({20210225},1,2);;
% batching{2}.runs                = {3;4};  % either empty or specific runs specified
% batching_type{2}                = 1; % 1 run by run, 2 session by session, 3 group by group
% batching{2}.range_of_dates      = 0;
% 


%% over all

run beh_run_analysis




