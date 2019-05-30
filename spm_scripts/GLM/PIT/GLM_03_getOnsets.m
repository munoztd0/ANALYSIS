function GLM_03_getOnsets()

% intended for REWOD PIT run

% get onsets for first control model (reward vs neutral)
% durations
% Simplified model on ONSETs 3*CS with modulator and grips as control
% last modified on APRIL 2019

%% define paths

homedir = '/home/REWOD/';
%homedir = '/home/cisa/CISA/REWOD';

mdldir        = fullfile (homedir, '/DATA/STUDY/MODELS/SPM');
sourcefiles   = fullfile(homedir, '/DATA/STUDY/CLEAN');
addpath (genpath(fullfile(homedir,'/ANALYSIS/my_tools')));

ana_name      = 'GLM-03';
%session       = {'second'};
task          = {'PIT'};
subj          = subID; %{'01';'02';'03';'04';'05';'06';'07';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'20';'21';'22';'23';'24';'25';'26'}; %doing it with 19 & 01?


%% create folder  
mkdir (fullfile (mdldir, char(task), ana_name)); % this is only because we have one run per task

%% extract and save data
for j = 1:length(task)
    
    taskX      = char(task(1));
    %sessionX  = char(session(j));
    
    for  i=1:length(subj)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load participants data
        subjX=[char(subj(i))];

        subjdir=fullfile(mdldir, char(task), ana_name,  ['sub-' subjX],'timing');
        mkdir (subjdir)

        cd (fullfile(sourcefiles,['sub-' subjX], 'func')); 
        behavfile = ['sub-' num2str(subjX) '_ses-second' '_task-' taskX '_run-01_events.mat'];
        fprintf('participant number: %s task: %s \n', subj{i}, task{1})
        disp(['file ' num2str(i) ' ' behavfile]);
        load (behavfile);
        
        %% FOR SPM
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets and durations for CS FOR RIM
        onsets.Reminder.trial          = RIM.ONSETS.trialstart;
        durations.Reminder.trial       = RIM.DURATIONS.trialstart;

        %replaced grip_frq by mob_effort
        modulators.Reminder.trial      = RIM.BEHAVIOR.mobilized_effort;
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets grips %%
        onsets.Reminder.grips           = RIM.ONSETS.grips;
        durations.Reminder.grips       = zeros (length(onsets.Reminder.grips),1);
        modulators.Reminder.grips      = ones  (length(onsets.Reminder.grips),1);
        
              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets and durations for CS FOR PE
        onsets.PartialExtinction.trial          = PE.ONSETS.trialstart;
        durations.PartialExtinction.trial       = PE.DURATIONS.trialstart;

        %replaced grip_frq by mob_effort
        modulators.PartialExtinction.trial      = PE.BEHAVIOR.mobilized_effort;

 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets grips %%
        onsets.PartialExtinction.grips           = PE.ONSETS.grips;
        durations.PartialExtinction.grips       = zeros (length(onsets.PartialExtinction.grips),1);
        modulators.PartialExtinction.grips      = ones  (length(onsets.PartialExtinction.grips),1);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets and durations for CS FOR PIT
        onsets.CS.CSp          = PIT.ONSETS.trialstart(strcmp ('CSplus', PIT.CONDITIONS));
        onsets.CS.CSm          = PIT.ONSETS.trialstart(strcmp ('CSminus', PIT.CONDITIONS));
        onsets.CS.Baseline     = PIT.ONSETS.trialstart(strcmp ('Baseline', PIT.CONDITIONS));
        
        durations.CS.CSp       = PIT.DURATIONS.trialstart(strcmp ('CSplus', PIT.CONDITIONS));
        durations.CS.CSm       = PIT.DURATIONS.trialstart(strcmp ('CSminus', PIT.CONDITIONS));
        durations.CS.Baseline  = PIT.DURATIONS.trialstart(strcmp ('Baseline', PIT.CONDITIONS));
        
        %replaced grip_frq by mob_effort
        modulators.CS.CSp      = BEHAVIOR.mobilized_effort(strcmp ('CSplus', PIT.CONDITIONS));
        modulators.CS.CSm      = BEHAVIOR.mobilized_effort(strcmp ('CSminus', PIT.CONDITIONS));
        modulators.CS.Baseline = BEHAVIOR.mobilized_effort(strcmp ('Baseline', PIT.CONDITIONS));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets grips %%?
        onsets.grips           = PIT.ONSETS.grips;
        durations.grips       = zeros (length(onsets.grips),1);
        modulators.grips      = ones  (length(onsets.grips),1);

       
        %% FOR FSL
        onsets.CS.CSp          = PIT.ONSETS.trialstart(strcmp ('CSplus', PIT.CONDITIONS));
        onsets.CS.CSm          = PIT.ONSETS.trialstart(strcmp ('CSminus', PIT.CONDITIONS));
        onsets.CS.Baseline     = PIT.ONSETS.trialstart(strcmp ('Baseline', PIT.CONDITIONS));
        
        % go in the directory where data will be saved
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
        cd (subjdir) %save all info in the participant directory
        
        % create text file with 3 colons: onsets, durations, parametric modulators
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        name = {'CS'; 'grips'; 'Reminder'; 'PartialExtinction'};
        
        for ii = 1:length(name)
            
            nameX = char(name(ii));
            
            if strcmp (nameX, 'CS')  % for structure that contains substuctures
                substr = {'CSp'; 'CSm'; 'Baseline'};% specify the substructures names
                
                for iii = 1:length(substr)
                    substrX = char(substr(iii));
                    nameXX  = [nameX '_' substrX]; % name that combines the structure and the substructures
                    % database with three rows of interest
                    database.(nameXX) = [num2cell(onsets.(nameX).(substrX)), num2cell(durations.(nameX).(substrX)), num2cell(modulators.(nameX).(substrX))];
                    % save the database in a txt file
                    fid = fopen ([ana_name '_task-' taskX '_' nameX '_' substrX '.txt'],'wt');
                    formatSpec = '%f\t%f\t%f\n';
                    [nrows,~] = size(database.(nameXX));
                    for row = 1:nrows
                        fprintf(fid,formatSpec,database.(nameXX){row,:});
                    end
                    fclose(fid);
                end
                
            elseif strcmp (nameX, 'PartialExtinction')  % for structure that contains substuctures
                     substr = {'trial'; 'grips'};% specify the substructures names

                for iii = 1:length(substr)
                    substrX = char(substr(iii));
                    nameXX  = [nameX '_' substrX]; % name that combines the structure and the substructures
                    % database with three rows of interest
                    database.(nameXX) = [num2cell(onsets.(nameX).(substrX)), num2cell(durations.(nameX).(substrX)), num2cell(modulators.(nameX).(substrX))];
                    % save the database in a txt file
                    fid = fopen ([ana_name '_task-' nameX '_' substrX '.txt'],'wt');
                    formatSpec = '%f\t%f\t%f\n';
                    [nrows,~] = size(database.(nameXX));
                    for row = 1:nrows
                        fprintf(fid,formatSpec,database.(nameXX){row,:});
                    end
                    fclose(fid);
                end
                
            elseif strcmp (nameX, 'Reminder')  % for structure that contains substuctures
                     substr = {'trial'; 'grips'};% specify the substructures names

                for iii = 1:length(substr)
                    substrX = char(substr(iii));
                    nameXX  = [nameX '_' substrX]; % name that combines the structure and the substructures
                    % database with three rows of interest
                    database.(nameXX) = [num2cell(onsets.(nameX).(substrX)), num2cell(durations.(nameX).(substrX)), num2cell(modulators.(nameX).(substrX))];
                    % save the database in a txt file
                    fid = fopen ([ana_name '_task-' nameX '_' substrX '.txt'],'wt');
                    formatSpec = '%f\t%f\t%f\n';
                    [nrows,~] = size(database.(nameXX));
                    for row = 1:nrows
                        fprintf(fid,formatSpec,database.(nameXX){row,:});
                    end
                    fclose(fid);
                end

          else
                % database with three rows of interest 
                database.(nameX) = [num2cell(onsets.(nameX)), num2cell(durations.(nameX)), num2cell(modulators.(nameX))];
                % save the database in a txt file
                fid = fopen ([ana_name '_task-' taskX '_' nameX '.txt'],'wt');
                formatSpec = '%f\t%f\t%f\n';
                [nrows,~] = size(database.(nameX));
                for row = 1:nrows
                    fprintf(fid,formatSpec,database.(nameX){row,:});
                end
                fclose(fid);
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save data
        mat_name = [ana_name '_task-' taskX '_onsets'];
        save (mat_name, 'onsets', 'durations', 'modulators')
  
        end
               
        

        
    end
    
end