function [files_arranged, file_as_i_want_it_cell] = bh_arrange_trials_eye_hand(folder_with_session_days, dates, batching)
dir_folder_with_session_days=dir(folder_with_session_days); % dir

if batching.range_of_dates == 1
    if dates==0
        all_files_in_base_path={dir_folder_with_session_days.name}; % all files from the main monkey folder
    else
        starting_date= dates(1);
        end_date=dates(2);
        folder_with_session_days_1= [folder_with_session_days filesep num2str(dates(1))];
        folder_with_session_days_2= [folder_with_session_days filesep num2str(dates(2))];
        all_files_in_base_path=[];
        ctr=1;
        for k=1: length(dir_folder_with_session_days)
            X=str2num(dir_folder_with_session_days(k).name);
            if X==dates(1) |  ( X<=  dates(2) & X >  dates(1))
                all_files_in_base_path{ctr}= dir_folder_with_session_days(k).name;
                ctr=ctr+1;
            end
        end
    end
else
    dates_temp=num2cell(dates);
    dates_temp_cell = cellfun(@num2str, dates_temp,'UniformOutput',0);  
    all_files_in_base_path=[];
    ctr=1;
    
    for k=1: length(dir_folder_with_session_days)
        for i=1:numel(dates)
            X=str2num(dir_folder_with_session_days(k).name);
            if dates(i) == X
                all_files_in_base_path{ctr}= dir_folder_with_session_days(k).name;
                ctr=ctr+1;
            end
        end
    end
end

if isempty(all_files_in_base_path)
disp(['no valid dates ' num2str(dates) ' inside folder ', folder_with_session_days])
files_arranged=[];
file_as_i_want_it_cell={};
    return
end

all_files_in_base_path = cellfun(@num2str,num2cell(unique(cellfun(@str2num, all_files_in_base_path))),'UniformOutput',0);

c=0;
files_arranged=[];
i_run=1;
for in_folders = 1:length(all_files_in_base_path)
    if length(all_files_in_base_path{in_folders})==8 % take all the folders named as yyyymmdd from the main monkey folder
        c=c+1;
        folder_to_look_at{c} = all_files_in_base_path{in_folders}; % start looping through session folders
        individual_day_folder = [folder_with_session_days filesep folder_to_look_at{c}]; % session of interest
        d_individual_day_folder=dir(individual_day_folder); % dir
        files_inside_session_folder={d_individual_day_folder.name}'; % files inside session folders
        for number_of_files = 1:length(files_inside_session_folder) % start looping within the session folder
            if length(files_inside_session_folder{number_of_files}) > 4 && strcmp(files_inside_session_folder{number_of_files}(end-3:end),'.mat')
                %load([individual_day_folder filesep files_inside_session_folder{number_of_files}])
                day_file.task_name(i_run,:)=[individual_day_folder filesep files_inside_session_folder{number_of_files}];
                file_as_i_want_it_cell(i_run,:)= {day_file.task_name(i_run,1:end-21) str2double(day_file.task_name(i_run,end-5:end-4))};
                i_run=i_run+1;
            end
        end
    end
    clc;
end

