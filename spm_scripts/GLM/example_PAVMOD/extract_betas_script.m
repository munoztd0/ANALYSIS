clear
clc

% which model?
ana_name = 'GLM-04o';

% path
dir_base  = '/Users/evapool/mountpoint/'; %ANALYSIS/fsl_script/ROI/GLM-03o'; %'home/evapool/PAVMOD/ANALYSIS/fsl_script/ROI/GLM-03p';
dir_data  = fullfile(dir_base, 'DATA','brain', 'MODELS', 'SPM', ana_name, 'group');
dir_output= fullfile(dir_base, 'DATA','brain', 'MODELS', 'SPM', ana_name, 'group');

% intialize spm 
spm('defaults','fmri');
spm_jobman('initcfg');

cd(dir_data)

sub_list = ['01';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';...
    '15';'16';'17';'18';'19';'20';'21';'22';'23';'24';'25';'26';'27';'28';'29';...
    '30'];

%list ROIs/clusters from which to extract betas
roi_dir = fullfile(dir_base, 'ANALYSIS','fsl_script', 'ROI', 'GLM-03o');
%roi_list = char(spm_select('FPList', roi_dir, '^.*\.img$'), spm_select('FPList', roi_dir, '^.*\.nii$'));
roi_list = char(spm_select('FPList', roi_dir, ['^'  '.*' 'nii']));


%list contrasts of interest
con_names = {'interValue.run2-run3'
    'run2.value-devalue'
    'run3.value-devalue'};

con_list = {'con-0008.nii,1' %interaction value run2-run3
    'con-0006.nii,1' % value-devalue run2
    'con-0007.nii,1'};% value-devalue run3


%loop across ROIs first
for r=1:size(roi_list,1)
    result = con_names';
    maskName = roi_list(r,:);
    v_mask = spm_vol(maskName); %extract voxels from ROIs as mask
    outputFile = [v_mask.fname(1:end-4) '_betas_' ana_name '.mat'];
    
    if exist(outputFile,'file')==0 %extract betas only if file doesn't exist    
        Y = spm_read_vols(v_mask);
        roi_volume_mm = length(find(Y > 0))*abs(det(v_mask.mat));
        clear Y;
        
        %loop across contrasts
        for c = 1:length(con_list)
            conName = con_names(c,:);
            % List of files to extract data from
            for s0 = 1:length(sub_list)  % select con image from every subject
                conDir = fullfile( dir_data, ['sub-' sub_list(s0,:) '_' char(con_list(c,:))]); %add Model directory if necessary
                images(s0,:) = conDir;
            end
            numImages = size(images,1);

            fprintf(1,[char(conName) ': Memory mapping images...\n']);
            for i=1:numImages
                v{i} = spm_vol(images(i,:));
                % Verify images are in same space
                if ~isequal(v{i}.dim(1:3),v{1}.dim(1:3))
                error('Images must have same dimensions.')
                end
                % Verify orientation/position are the same
                if ~isequal(v{i}.mat,v{1}.mat)
                error('Images must have same orientation/position.')
                end
            end
            fprintf(1,'done.\n');

            [Y, XYZmm] = spm_read_vols(v{1});
            clear Y

            XYZmask = inv(v_mask.mat)*([XYZmm; ones(1,size(XYZmm,2))]);
            ind = find(spm_sample_vol(v_mask, XYZmask(1,:), XYZmask(2,:), XYZmask(3,:),0) > 0);

            ind = ind(:)';
            vals = [];
            for j=1:length(v)
                Y = spm_read_vols(v{j});
                % rows:  images.  cols:  voxels
                vals = [vals; Y(ind)];
            end

            % Consider only values which are finite for all images
            vals = vals(:, find(all(isfinite(vals))));
            intersection_volume_mm = size(vals,2)*abs(det(v{1}.mat));
            if intersection_volume_mm==0
                error('No voxels in ROI are in-brain for all images to be sampled.');
            end
            m = mean(vals, 2);

            result(2:length(v)+1,c) = num2cell(m); 

            save(outputFile,'result');

        end
 
    end
         
end