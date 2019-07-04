function GLM_MF_13_ndLevel_covariate()

% BY THE END THIS SCRIPT NEEDS TO BE COMBINED IN A SINGLE SCRIPT


% does t-test and full_factorial
do_ttest = 0;
do_covariate = 1;
removesub = 'sub-28'; %'sub-24'; % which sub do we want to remove

%% define path

homedir = '/home/eva/PAVMOD/';

%homedir = '/Users/evapool/mountpoint/';

funcdir  = fullfile (homedir, '/DATA/brain/cleanBIDS');% directory with  post processed functional scans
mdldir   = fullfile (homedir, '/DATA/brain/MODELS/SPM');% mdl directory (timing and outputs of the analysis)
covdir   = fullfile (homedir, '/DATA/group_covariates'); % director with the extracted second level covariates

name_ana = 'GLM-MF-13'; % output folder for this analysis
groupdir = fullfile (mdldir,name_ana, 'group/');


%% specify spm param
addpath('/usr/local/matlab/R2014a/toolbox/spm12b');
addpath ([homedir '/ANALYSIS/spm_script/GLM/dependencies']);
spm('Defaults','fMRI');
spm_jobman('initcfg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define constrasts and constrasts names
if do_covariate
    
    % covariate of interest name become folder
    covariateNames = {'RT_devaluation_efficiency'
        'CS_devaluation_efficiency'
        'US_devaluation_efficiency'
        'PU_devaluation_efficiency'};

    % These contrast names become sub-folders
    contrastNames = {
        'RUN123.CS.INT1' %4
        'RUN12.CS.INT1'% 5
        'RUN12.ANT.INT1'%6
        'RUN123.ANT.INT1'%7
        'RUN12.ALL.INT1'%8
        'RUN123.ALL.INT1'%9
        'RUN12.CS.INT2'%10
        'RUN123.CS.INT2'%11
        'RUN12.ANT.INT2'%12
        'RUN123.ANT.INT2'%13
        'RUN12.ALL.INT2'%14
        'RUN123.ALL.INT2'%15
        'RUN12.CS.INT3'%16
        'RUN12.ANT.INT3'%17
        'RUN12.ALL.INT3'};%18;
    
    conImages = {
        'con-0004'
        'con-0005'
        'con-0006'
        'con-0007'
        'con-0008'
        'con-0009'
        'con-0010'
        'con-0011'
        'con-0012'
        'con-0013'
        'con-0014'
        'con-0015'
        'con-0016'
        'con-0017'
        'con-0018'};
    
    %% prepare batch for each contrasts
    
    for c = 1:length(covariateNames)
        
        covariateX = covariateNames{c};
        
        filename = fullfile(covdir,[covariateX '.txt']);
        delimiterIn = ' ';
        headerlinesIn = 1;
        C = importdata(filename,delimiterIn,headerlinesIn);
        
        cov.ID   = C.data(:,1);
        cov.data = C.data(:,2);
       
        if removesub
            idx            = str2double(removesub(5:end));
            torm           = find(cov.ID==idx);
            
            cov.ID(torm)   = [];
            cov.data(torm) = [];
        end
        
        for n = 1:length(contrastNames)
            
            clear matlabbatch
            
            conImageX = conImages{n};
            contrastX = contrastNames{n};
            
            if removesub %
                contrastFolder = fullfile (groupdir, 'covariate', covariateX, removesub, contrastX);
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
            
            if removesub % remove subject from analysis
                
                disp(['removing subject: ' removesub]);
                allsub = matlabbatch{1}.spm.stats.factorial_design.des.t1.scans; % let's put this in a smaller variable
                idx = (regexp(allsub,removesub)); % find string containing the sub id
                idxtoRemove = find(~cellfun(@isempty,idx)); % get the index of that string
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans(idxtoRemove) = []; % remove the string from the scans selected for the analysis
                
            end
            
            matlabbatch{1}.spm.stats.factorial_design.cov.c      = cov.data;
            matlabbatch{1}.spm.stats.factorial_design.cov.cname  = covariateX;
            matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;
            
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
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



end