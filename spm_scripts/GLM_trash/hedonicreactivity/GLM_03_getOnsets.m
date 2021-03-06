function GLM_03_getOnsets()
%watchout
% intended for REWOD HEDONIC reactivity 

% get onsets for first control model (2nd level modulators)
% Stick functions
% Simplified model on ONSETs 7 (STARTTRIAL, 2*odor with modulator (liking
% ratings) 3*questions)
% last modified on MAY 2019

%% define paths

homedir = '/home/REWOD/';
%homedir = '/home/cisa/CISA/REWOD';

mdldir        = fullfile (homedir, '/DATA/STUDY/MODELS/SPM');
sourcefiles   = fullfile(homedir, '/DATA/STUDY/CLEAN');
addpath (genpath(fullfile(homedir,'/ANALYSIS/my_tools')));

ana_name      = 'GLM-03';
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
        durations.trialstart    = DURATIONS.trialstart;
        modulators.trialstart = ones (length(onsets.trialstart),1); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets and durations for odor valveopen?
        onsets.odor.reward      = ONSETS.sniffSignalOnset(strcmp ('chocolate', CONDITIONS));
        onsets.odor.neutral     = ONSETS.sniffSignalOnset(strcmp ('neutral', CONDITIONS));
        onsets.odor.control     = ONSETS.sniffSignalOnset(strcmp ('empty', CONDITIONS));

        
        %get durations
        durations.odor.reward   = DURATIONS.trialstart(strcmp ('chocolate', CONDITIONS));
        durations.odor.neutral   = DURATIONS.trialstart(strcmp ('neutral', CONDITIONS));
        durations.odor.control   = DURATIONS.trialstart(strcmp ('empty', CONDITIONS));
       
        
        %for liking
        modulators.odor.reward.lik  = BEHAVIOR.liking (strcmp ('chocolate', CONDITIONS));
        modulators.odor.neutral.lik = BEHAVIOR.liking (strcmp ('neutral', CONDITIONS));
        modulators.odor.control.lik = BEHAVIOR.liking (strcmp ('empty', CONDITIONS));

        
        %for intensity
        modulators.odor.reward.int  = BEHAVIOR.intensity (strcmp ('chocolate', CONDITIONS));
        modulators.odor.neutral.int = BEHAVIOR.intensity (strcmp ('neutral', CONDITIONS));
        modulators.odor.control.int = BEHAVIOR.intensity (strcmp ('empty', CONDITIONS));
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get onsets and duration questions
        onsets.liking            = ONSETS.liking;
        durations.liking         = DURATIONS.liking;
        modulators.liking        = ones (length(onsets.liking),1);
        
        onsets.intensity         = ONSETS.intensity;
        durations.intensity      = DURATIONS.intensity;
        modulators.intensity     = ones (length(onsets.liking),1);
        

        


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
                substr = {'reward'; 'control'; 'neutral'};% specify the substructures names ?neutral?
                subsubstr = {'lik'; 'int'};
                for iii = 1:length(substr)
                    substrX = char(substr(iii));
                    for iiii =  1:length(subsubstr)
                        subsubstrX = char(subsubstr(iiii));
                        nameXX  = [nameX '_' substrX '_' subsubstrX]; % name that combines the structure and the substructures
                        % database with three rows of interest
                        database.(nameXX) = [num2cell(onsets.(nameX).(substrX)), num2cell(durations.(nameX).(substrX)), num2cell(modulators.(nameX).(substrX).(subsubstrX))];
                        % save the database in a txt file
                        fid = fopen ([ana_name '_task-' taskX '_' nameX '_' substrX '_' subsubstrX '.txt'],'wt');
                        formatSpec = '%f\t%f\t%f\n';
                        [nrows,~] = size(database.(nameXX));
                        for row = 1:nrows
                            fprintf(fid,formatSpec,database.(nameXX){row,:});
                        end
                        fclose(fid);
                    end
                end
             

          else
                % database with three rows of interest %%%% ADD MODULATORS
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

              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save data
        mat_name = [ana_name '_task-' taskX '_onsets'];
        save (mat_name, 'onsets', 'durations', 'modulators')
        end



    end
end

