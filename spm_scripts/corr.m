%sub-26_ses-second_task-hedonic_run-01_events.mat'


% intended for REWOD HEDONIC

% get liking and intensity in .txt

%% define paths

%homedir = '/home/REWOD/';
homedir = '/home/cisa/CISA/REWOD';

%mdldir        = fullfile (homedir, '/DATA/STUDY/MODELS/SPM');
sourcefiles   = fullfile(homedir, '/DATA/STUDY/CLEAN');
addpath (genpath(fullfile(homedir,'/ANALYSIS/my_tools')));

%ana_name      = 'GLM-02';
%session       = {'second'};
task          = {'hedonic'};
subj          = {'01';'02';'03';'04';'05';'06';'07';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'20';'21';'22';'23';'24';'25';'26'}; %doing it with 19 & 01?

%% create folder
%mkdir (fullfile (mdldir, char(task), ana_name)); % this is only because we have one run per task

%% extract and save data
for j = 1:length(task) %useless in this case

    taskX      = char(task(1));
    %sessionX  = char(session(j));

    for  i=1:length(subj)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load participants data
        subjX=[char(subj(i))];

        %subjdir=fullfile(mdldir, char(task), ana_name,  ['sub-' subjX],'timing');
        %mkdir (subjdir)

        cd (fullfile(sourcefiles,['sub-' subjX], 'func')); 
        behavfile = ['sub-' num2str(subjX) '_ses-second' '_task-' taskX '_run-01_events.mat'];
        fprintf('participant number: %s task: %s \n', subj{i}, task{1})
        disp(['file ' num2str(i) ' ' behavfile]);
        load (behavfile);
        
        
        
         database.corr = [num2cell(BEHAVIOR.intensity), num2cell(BEHAVIOR.liking)];
         % save the database in a txt file
         fid = fopen ('corr_task-hedonic.txt','wt');
         formatSpec = '%f\t%f\n';
         [nrows,~] = size(database.corr);
         for row = 1:nrows
             fprintf(fid,formatSpec,database.corr{row,:});
         end
         fclose(fid);
     end

end