%function GLM_03_ndLevel()

%HEDONIC

%does t-test and full_factorial
do_ttest = 1;
remove = 1;
removesub = {'sub-23'} ; %because 23 neutral_int zero var
removedsub = '23';

%% define path

%homedir = '/home/REWOD';
homedir = '/home/cisa/REWOD';
%homedir = '/Users/davidmunoz/REWOD';
mdldir   = fullfile (homedir, '/DATA/STUDY/MODELS/SPM/hedonic');% mdl directory (timing and outputs of the analysis)
name_ana = 'GLM-03c'; 
groupdir = fullfile (mdldir,name_ana, 'group/');

%% specify spm param
%addpath /usr/local/external_toolboxes/spm12/ ;
addpath /usr/local/MATLAB/R2018a/spm12 ;
addpath ([homedir '/ANALYSIS/spm_scripts/GLM/dependencies']);

spm('Defaults','fMRI');
spm_jobman('initcfg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define constrasts and constrasts names
if do_ttest
    
    % These contrast names become folders
    contrastNames = {'reward-control'%1
        'reward-neutral'%2
        'Odor-NoOdor'%3
        'odor_presence'%4
        'reward_lik-control'%5
        'reward_int-control'%6
        'reward_lik-neutral_lik'%7
        'reward_int-neutral_int'%8
        'Odor_lik-NoOdor'%9
        'Odor_int-NoOdor'%10
        'odor_lik_presence'%11
        'odor_int_presence' %12
        'reward_lik-control_lik'%13
        'reward_int-control_int'}; %14

   
    
    conImages = {'con_0001'
        'con_0002'
        'con_0003'
        'con_0004'
        'con_0005'
        'con_0006'
        'con_0007'
        'con_0008'
        'con_0009'
        'con_0010'
        'con_0011'
        'con_0012'
        'con_0013'
        'con_0014'};
    
    
      %% prepare batch for each contrasts
    
    for n = 1:length(contrastNames)
        
        clear matlabbatch
        
        conImageX = conImages{n};
        contrastX = contrastNames{n};
        
        if remove
            contrastFolder = fullfile (groupdir, 'ttests', ['removing-' removedsub], contrastX);
        else
            contrastFolder = fullfile (groupdir, 'ttests', 'all', contrastX);
        end
        
        mkdir(contrastFolder);
        
        % create the group level spm file
        matlabbatch{1}.spm.stats.factorial_design.dir = {contrastFolder}; % directory
        
        %  FORMAT [dirs] = spm_select('List',direc,'dir',filt)
        conAll     = spm_select('List',groupdir,['^'  '.*' conImageX '.nii']); % select constrasts 
     
        for j =1:length(conAll)
            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{j,1} = [groupdir conAll(j,:) ',1'];
        end
        
        if remove % remove subject from analysis
            disp(['removing subject: ' removedsub]);
            allsub = matlabbatch{1}.spm.stats.factorial_design.des.t1.scans; % let's put this in a smaller variable
            for i = 1:length(removesub)
                idx = (regexp(allsub,removesub{i})); % find string containing the sub id
                idxtoRemove = find(~cellfun(@isempty,idx)); % get the index of that string
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans(idxtoRemove) = []; % remove the string from the scans selected for the analysis
                allsub = matlabbatch{1}.spm.stats.factorial_design.des.t1.scans;
            end
               
        end
        
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        % estimate design matrix
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {[contrastFolder  '/SPM.mat']};
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        
        % specify one sample tconstrast
        matlabbatch{3}.spm.stats.con.spmmat(1)               = {[contrastFolder  '/SPM.mat']};
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name     = contrastX (1:end);
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights  = [1];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep  = 'none';
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name     = ['Neg ' contrastX(1:end)];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights  = [-1];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep  = 'none';
        
        spm_jobman('run',matlabbatch)
        
    end
end
%end
