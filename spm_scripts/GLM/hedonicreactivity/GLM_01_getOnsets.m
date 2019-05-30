function GLM_01_getOnsets()
%watchout
% intended for REWOD hedonic reactivity run

% get onsets for first control model (reward vs neutral)
% Stick functions
% Simplified model on ONSETs 7 (STARTTRIAL, 2*odor with modulator (liking
% ratings) 3*questions 1 RINSE)
% last modified on MARCH 2018

%% define paths

homedir = '/home/REWOD/';
%homedir = '/home/cisa/CISA/REWOD';

mdldir        = fullfile (homedir, '/DATA/STUDY/MODELS/SPM');
sourcefiles   = fullfile(homedir, '/DATA/STUDY/CLEAN');
addpath (genpath(fullfile(homedir,'/ANALYSIS/my_tools')));

ana_name      = 'GLM-01';
%session       = {'second'};
task          = {'hedonic'};
subj          = subID; %{'01';'02';'03';'04';'05';'06';'07';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'20';'21';'22';'23';'24';'25';'26'}; %doing it with 19 & 01?

%% create folder
mkdir (fullfile (mdldir, char(task), ana_name)); % this is only because we have one run per task

%% extract and save data
%for j = 1:length(task)

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
        % Get onsets and durations for start
        onsets.trialstart    = ONSETS.trialstart;
        durations.trialstart = zeros (length(onsets.trialstart),1); %??
        %durations.trialstart    = DURATIONS.trialstart;
        modulators.trialstart = ones (length(onsets.trialstart),1); %??

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets and durations for odor valveopen?
        %onsets.taste.reward      = ONSETS.liquid(strcmp ('MilkShake', CONDITIONS));
        onsets.odor.reward      = ONSETS.sniffSignalOnset(strcmp ('chocolate', CONDITIONS));
        onsets.odor.neutral     = ONSETS.sniffSignalOnset(strcmp ('neutral', CONDITIONS));
        onsets.odor.control     = ONSETS.sniffSignalOnset(strcmp ('empty', CONDITIONS));

        
        %get zero values for durations (stick functions)
        durations.odor.reward   = zeros (length(onsets.odor.reward),1);
        durations.odor.neutral  = zeros (length(onsets.odor.neutral),1);
        durations.odor.control  = zeros (length(onsets.odor.control),1);

       
        modulators.odor.reward  = BEHAVIOR.liking (strcmp ('chocolate', CONDITIONS));
        modulators.odor.neutral = BEHAVIOR.liking (strcmp ('neutral', CONDITIONS));
        modulators.odor.control = BEHAVIOR.liking (strcmp ('empty', CONDITIONS));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets and duration questions
        onsets.liking            = ONSETS.liking;
        durations.liking         = DURATIONS.liking;
        modulators.liking        = ones (length(onsets.liking),1);

        onsets.intensity         = ONSETS.intensity;
        durations.intensity      = DURATIONS.intensity;
        modulators.intensity     = ones (length(onsets.intensity),1);



        %% FOR FSL

        % go in the directory where data will be saved
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
        cd (subjdir) % let's save all info in the participant directory

        % create text file with 3 colons: onsets, durations, paretric modulators
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        name = {'trialstart'; 'odor'; 'liking'; 'intensity'}; % 'familiarity'; 'rinse'};

        for ii = 1:length(name)

            nameX = char(name(ii));

            if strcmp (nameX, 'odor')  % for structure that contains substuctures
                substr = {'reward'; 'control'; 'neutral'};% specify the substructures names

                for iii = 1:length(substr)
                    substrX = char(substr(iii));
                    nameXX  = [nameX '_' substrX]; % name that combines the structure and the substructures
                    % database with three rows of interest
                    database.(nameXX) = [num2cell(onsets.(nameX).(substrX)), num2cell(durations.(nameX).(substrX)), num2cell(modulators.(nameX).(substrX))];
                    % save the database in a txt file
                    fid = fopen ([ana_name '_task-' taskX '_' nameX '_' substrX '.txt'],'wt');
                    formatSpec = '%f\t%f\t%d\n';
                    [nrows,~] = size(database.(nameXX));
                    for row = 1:nrows
                        fprintf(fid,formatSpec,database.(nameXX){row,:});
                    end
                    fclose(fid);
                end

          else
                % database with three rows of interest %%%% ADD MODULATORS
                database.(nameX) = [num2cell(onsets.(nameX)), num2cell(durations.(nameX)), num2cell(modulators.(nameX))];
                % save the database in a txt file
                fid = fopen ([ana_name '_task-' taskX '_' nameX '.txt'],'wt');
                formatSpec = '%f\t%f\t%d\n';
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


