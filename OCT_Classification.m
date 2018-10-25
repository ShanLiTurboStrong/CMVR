clear all;clc;
% set path for training samples and test samples directory
addpath('large_scale_svm');
addpath('sift');
addpath(genpath('sparse_coding'));
% directory setup for training samples
img_dir = 'image';% directory for dataset images
flattened_img_dir = 'Flattened_image';
data_dir = 'data';                  % directory to save the sift features of the chosen dataset
dataSet = 'OCTT';


%% parameter setting
% sift descriptor extraction
skip_cal_sift =true;                % if 'skip_cal_sift' is true, skip extracting sift feature patches for training
gridSpacing = 8;
patchSize = 16;
maxImSize = 1000;                    % The Caltech101 is 300,remember to change it
nrml_threshold = 1;                  % low contrast region normalization threshold (descriptor length)
% dictionary training for sparse coding
skip_dic_training = true;           % if 'skip_dic_training' is true, skip dictionary training
nBases = 1024;
nsmp =60000;                      % the total patches
beta = 1e-5;                        % a small regularization for stablizing sparse coding
num_iters = 1;                     % the iteration times. why 50? it takes too much time.
% feature pooling parameters
pyramid = [1, 2, 4];                % spatial block number on each level of the pyramid
gamma = 0.15;
knn = 200;                          % find the k-nearest neighbors for approximate sparse coding
% ScSPM for training sets.
skip_cal_sparse= true;
skip_svm=true;
%%============================================ Codebook=====================================================================%

rt_img_dir = fullfile(img_dir, dataSet);
rt_data_dir = fullfile(data_dir, dataSet);
rt_flattened_img_dir = fullfile(flattened_img_dir, dataSet);
%%=================================================================================================================%



%==================calculate sift features or retrieve the database directory==============%
if skip_cal_sift
    database = retr_database_dir(rt_data_dir);
else
    [database, lenStat, patchNumberToatal, patchNumberPerImage] = CalculateSiftDescriptor(rt_flattened_img_dir, rt_data_dir, gridSpacing, patchSize, maxImSize, nrml_threshold);
    save('log.mat', 'patchNumberToatal', 'patchNumberPerImage');
    %[database,lenStat] = CalculateLBPSiftDescriptor(rt_img_dir, rt_data_dir, gridSpacing, patchSize, maxImSize, mapping,nrml_threshold);
end
%=====================================================================================%

%===========================% for exclude three persons' images out for cross-validation========================%
NUM_PATIENT = 3;
% randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n inclusive.
Test(1,:)=randperm(NUM_PATIENT,NUM_PATIENT);%AMD
Test(2,:)=randperm(NUM_PATIENT,NUM_PATIENT);%DME
Test(3,:)=randperm(NUM_PATIENT,NUM_PATIENT);%NOR

Tim=database.imnum(1,Test(1,NUM_PATIENT))+database.imnum(2,Test(2,NUM_PATIENT))+database.imnum(3,Test(3,NUM_PATIENT));
num_img=sum(sum(database.imnum))-Tim;
%==============================train the k-dimension(B) base for sparse coding=======================================%
Bpath = ['dictionary/dict_' dataSet '_' num2str(nBases) '.mat'];
if skip_dic_training,
    load(Bpath);
else
    X = rand_sampling(database,nsmp,Test, num_img);
    [B, S, stat] = reg_sparse_coding(X, nBases, eye(nBases), beta, gamma, num_iters);
    save(Bpath, 'B', 'S', 'stat');
end
nBases = size(B, 2);% size of the dictionary
%==================================================================================================%




%======================================== calculate the sparse coding feature==========================================%
Spath=['dictionary/Sparse_' dataSet '_' num2str(nBases) '.mat'];
if skip_cal_sparse
    load(Spath);
else
    dimFea = sum(nBases*pyramid.^2);
    sc_fea = zeros(dimFea, num_img);
    sc_label = zeros(num_img, 1);

    disp('==================================================');
    fprintf('Calculating the sparse coding feature...\n');
    fprintf('Regularization parameter: %f\n', gamma);
    disp('==================================================');
    [l,n]=size(Test);
    iter1=0;
    for i=1:l
        for j=1:n-1
            for m=1:database.imnum(i,Test(i,j))
                iter1 = iter1+1;
                if ~mod(iter1, 50),
                    fprintf('.\n');
                else
                    fprintf('.');
                end;
                fpath = database.path{i,Test(i,j),m};
                load(fpath);
                if knn,
                    sc_fea(:, iter1) = sc_approx_pooling(feaSet, B, pyramid, gamma, knn);
                else
                    sc_fea(:, iter1) = sc_pooling(feaSet, B, pyramid, gamma);
                end
                sc_label(iter1) = feaSet.label;
            end
        end
    end
    save(Spath, 'sc_label', 'sc_fea');
end
%===============================================================================================================%




%============================= train SVM using sparse code and label of each image=============================================%

Wpath=['dictionary/svm.mat'];
if skip_svm
    load(Wpath);
else
    lambda = 0.1;                       % regularization parameter for w
    TSTART=tic;
    [w, b, class_name] = li2nsvm_multiclass_lbfgs(sc_fea', sc_label, lambda);
    save(Wpath,'w','b','class_name');
end
%===============================================================================================================%





%===================================validation for the 3 persons' imageSet picked up earlier==================================================%
for l=1:database.nclass
    P{l}=database.path{l,Test(l,NUM_PATIENT),1};
    index=strfind(P{l},'/');
    P{l}=P{l}(index(2)+1:index(4)-1); % AMD/AMD1, ...
end

for ii=1:database.nclass
    subname = P{ii};
    Results=zeros(3,1);
    fprintf('Test Sample: %s\n',subname);
    if ~strcmp(subname, '.DS_Store') & ~strcmp(subname, '.') & ~strcmp(subname, '..')
        frames = dir(fullfile(rt_flattened_img_dir, subname,'*.tif'));
        c_num = length(frames);
        fprintf('length of frame= %d\n', c_num)
        for jj = 1:c_num
            imgpath = fullfile(rt_flattened_img_dir, subname, frames(jj).name);
            fprintf(frames(jj).name)
            I = imread(imgpath);%I(I>=255)=0;I = im2double(I);
            fprintf('Test Sample: %s_%d_%s\n',subname,jj,frames(jj).name);
            %        fprintf('.');
            % calculate sift features
            [feaSet_test,lenStat_test]=CalculateSiftDescriptor_Test(I,gridSpacing, patchSize, maxImSize, nrml_threshold);
            %[feaSet_test]=CalculateLBPSiftDescriptor_Test(I, gridSpacing, patchSize,mapping,nrml_threshold);
            % calculate the sparse coding feature
            if knn,
                sc_fea_test= sc_approx_pooling(feaSet_test, B, pyramid, gamma, knn);
            else
                sc_fea_test= sc_pooling(feaSet_test, B, pyramid, gamma);
            end


            [C, Y] = li2nsvm_multiclass_fwd(sc_fea_test', w, b, class_name);
            if C==1
                result='AMD';
                Results(1)=Results(1)+1;
            else
                if C==2
                    result='DME';
                    Results(2)=Results(2)+1;
                else
                    result='NOR';
                    Results(3)=Results(3)+1;
                end
            end
            fprintf('Result: %s\n',result);
        end
    end
    fprintf('%s result:AMD:%d, DME:%d, NOR:%d \n',subname,Results(1),Results(2),Results(3));
    ResultA{ii}=sprintf('%s result:AMD:%d, DME:%d, NOR:%d \n',subname,Results(1),Results(2),Results(3));
end
%============================================================================================================================%
save('allWorkspaceResult_Hospital(time-logged).mat');
save('Result.mat','ResultA');
