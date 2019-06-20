%function GLM_04_ndLevel
% with covariate not zscored
% PIT

% does t-test and full_factorial
do_ttest = 1;
do_covariate = 1;
remove = 1; %watcha
removesub = {'sub-24'}; %'sub-10';
removedsub = '24'; 
%removesub = 'sub-24'; % % which sub do we want to remove

%% define path

%homedir = '/home/REWOD';
homedir = '/home/cisa/CISA/REWOD';

mdldir   = fullfile(homedir, '/DATA/STUDY/MODELS/SPM/PIT');% mdl directory (timing and outputs of the analysis)
funcdir  = fullfile(homedir, '/DATA/STUDY/CLEAN');% directory with  post processed functional scans
covdir   = fullfile (homedir, '/DATA/STUDY/MODELS/SPM/PIT/GLM-04/group_covariates'); % director with the extracted second level covariates

name_ana = 'GLM-04'; % output folder for this analysis
groupdir = fullfile (mdldir,name_ana, 'group/');


%% specify spm param
%addpath('/usr/local/external_toolboxes/spm12/');
addpath /usr/local/MATLAB/R2018a/spm12 ;
addpath ([homedir '/ANALYSIS/spm_scripts/GLM/dependencies']);
spm('Defaults','fMRI');
spm_jobman('initcfg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define constrasts and constrasts names
if do_covariate
    
    % covariate of interest name become folder 
    covariateNames = {'CSp-CSm_effort' %1
        'CSp-Baseline_effort' %2
        'CSm-Baseline_effort'%3
        'CSp-CSm&Baseline_effort'}; %4

    % These contrast names become sub-folders
    contrastNames = {'CSp-CSm' %1
        'CSp-Baseline'% 2
        'CSp-CSm&Baseline'%3
        'grips'%4
        'CSm-Baseline'}; %5
    
    conImages = {
        'con-0001'
        'con-0002'
        'con-0003'
        'con-0004'
        'con-0005'};
    
    %% prepare batch for each contrasts
    
    for c = 1:length(covariateNames)
        
        covariateX = covariateNames{c};
        
        filename = fullfile(covdir, [covariateX '.txt']);
        delimiterIn = '\t';
        headerlinesIn = 1;
        C = importdata(filename,delimiterIn,headerlinesIn); %importdata
        
        cov.ID   = C.data(:,1);
        cov.data = C.data(:,2);
        
   

        if remove
            for i = 1:length(removesub)
                idx            = str2double(removesub{i}(5:end));
                torm           = find(cov.ID==idx);

                cov.ID(torm)   = [];
                cov.data(torm) = [];
            end
        end
        
        for n = 1:length(contrastNames)
            
            clear matlabbatch
            
            conImageX = conImages{n};
            contrastX = contrastNames{n};
            
            if remove 
                contrastFolder = fullfile (groupdir, 'covariate', covariateX, ['removing-' removedsub], contrastX);
            else
                contrastFolder = fullfile (groupdir, 'covariate', covariateX, 'all', contrastX);
            end
            
            mkdir(contrastFolder);
            
            % create the group level spm file
            matlabbatch{1}.spm.stats.factorial_design.dir = {contrastFolder}; % directory
            
            % select contrasts only for participants that have the behavioral covariate
            for s = 1:length(cov.ID)
                cov.IDX      = cov.ID(s);
           
                Scue = deblank(['sub-' sprintf('%02d ', cov.IDX)]);
                conAll (s,:) = spm_select('List',groupdir,['^' Scue '.*' conImageX '.nii']); % select constrasts
            end
            
            for j =1:size(conAll,1)
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{j,1} = [groupdir conAll(j,:) ',1'];
            end
            
            if remove % remove subject from analysis
                disp(['removing subject: ' removedsub]);
                allsub = matlabbatch{1}.spm.stats.factorial_design.des.t1.scans; % let's put this in a smaller variable
                for ii = 1:length(removesub)
                    idx = (regexp(allsub,removesub{ii})); % find string containing the sub id
                    idxtoRemove = find(~cellfun(@isempty,idx)); % get the index of that string
                    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans(idxtoRemove) = []; % remove the string from the scans selected for the analysis
                    allsub = matlabbatch{1}.spm.stats.factorial_design.des.t1.scans;
                end
            end
            
            matlabbatch{1}.spm.stats.factorial_design.cov.c      = cov.data;
            matlabbatch{1}.spm.stats.factorial_design.cov.cname  = covariateX;
            matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;
            
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1; %%??
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            % extimate design matrix
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {[contrastFolder  '/SPM.mat']};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            
            % specify one sample tconstrast
            matlabbatch{3}.spm.stats.con.spmmat(1)                = {[contrastFolder  '/SPM.mat']};
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name     = contrastX (1:end);
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights  = [1 0]; % in the covariate the second colon is the one of interest
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep  = 'none';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name     = ['Neg ' contrastX(1:end)];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights  = [-1 0]; % in the covariate the second colon is the one of interest
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep  = 'none';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.name     = covariateX (1:end);
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights  = [0 1]; % in the covariate the second colon is the one of interest
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep  = 'none';
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.name     = ['Neg ' covariateX(1:end)];
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights  = [0 -1];
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep  = 'none';
            
            disp ('***************************************************************') 
            disp (['running batch for: '  contrastX ': ' covariateX] ) 
            disp ('***************************************************************') 
               
            spm_jobman('run',matlabbatch)
            
        end
    end
end

clear all


%end