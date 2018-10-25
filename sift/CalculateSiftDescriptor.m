function [database, lenStat, patchNumberToatal, patchNumberPerImage] = CalculateSiftDescriptor(rt_img_dir, rt_data_dir, gridSpacing, patchSize, maxImSize, nrml_threshold)
%==========================================================================

disp('Extracting SIFT features...');
subfolders = dir(rt_img_dir);%OCTT

siftLens = [];

database = [];

database.imnum = 0; % total image number of the database
database.cname = {}; % name of each class
database.path = {}; % contain the pathes for each image of each class
database.nclass = 0;
c_num_sum=0;
patchNumberToatal=0;
patchNumberPerImage=0;
for ii = 1:length(subfolders),
    subname = subfolders(ii).name;%AMD/DME/NOR
    %     disp(subfolders(ii).name);
    if ~strcmp(subname, '.DS_Store')& ~strcmp(subname, '.') & ~strcmp(subname, '..'),
        database.nclass = database.nclass + 1;
        database.cname{database.nclass} = subname;
        subfolders1 = dir([rt_img_dir,'/',subname]);%AMD1~15//DME1~15//NOR1~15
        count=0;
        for kk=1:length(subfolders1),
            subname1 = subfolders1(kk).name;
            
            if ~strcmp(subname1, '.DS_Store')& ~strcmp(subname1, '.') & ~strcmp(subname1, '..'),
                count=count+1;
                frames = dir(fullfile(rt_img_dir, subname,subname1, '*.tif'));
                c_num = length(frames);
                c_num_sum=c_num_sum+c_num;
                database.imnum(database.nclass,count) = c_num;
                
                %                     disp(num2str(database.nclass));
                %                     disp(num2str(count));
                %                     disp(num2str(fullfile(rt_img_dir, subname,subname1, '*.tif')));
                %                     disp(num2str(database.imnum(database.nclass,count)));
                %
                
                
                
                siftpath = fullfile(rt_data_dir, subname,subname1);
                if ~isdir(siftpath),
                    mkdir(siftpath);
                end;
                for jj = 1:c_num,
                    imgpath = fullfile(rt_img_dir, subname,subname1, frames(jj).name);
                    I = imread(imgpath);
                    if ndims(I) == 3,
                        I = im2double(rgb2gray(I));
                    else
                        I(I>=255)=0;
                        I = im2double(I);
                    end;
                    %                         %   BM3D       Read a grayscale image and scale its intensities in range [0,1]
                    %                         y = I;
                    %                         % Generate the same seed used in the experimental results of [1]
                    %                         randn('seed', 0);
                    %                         %rng(0,'v4')
                    %                         % Standard deviation of the noise --- corresponding to intensity
                    %                         %  range [0,255], despite that the input was scaled in [0,1]
                    %                         sigma = 25;%initial is 25, lasttime is 35
                    %                         % Add the AWGN with zero mean and standard deviation 'sigma'
                    %                         z = y + (sigma/255)*randn(size(y));
                    %                         % Denoise 'z'. The denoised image is 'y_est', and 'NA = 1' because
                    %                         %  the true image was not provided
                    %                         [NA, y_est] = BM3D(1, z, sigma);
                    %                         I=y_est;
                    
                    
                    
                    %                         if ndims(I) == 3,
                    %                             I = im2double(rgb2gray(I));
                    %                         else
                    %                             I(I>=255)=0;
                    %                             I = im2double(I);
                    %                         end;
                    % %                         [im_h,im_w] = size(I);
                    % %                         w=fspecial('gaussian');
                    % %                         I=imfilter(I,w);
                    % %                         gridnum=128;
                    % %                         num=floor(im_w/gridnum);
                    % %                         right=im_h*ones(1,gridnum);
                    % %                         for ii=1:gridnum
                    % %                             fea=sum(I(:,(ii-1)*num+1:ii*num),2);
                    % %                             [p]=max(fea);
                    % %                             for i=im_h:-1:1
                    % %                                 if fea(i)>p/2
                    % %                                     right(ii)=i;
                    % %                                     break;
                    % %                                 end
                    % %                             end
                    % %                         end
                    % %                         xq=1:gridnum;
                    % %                         f1=1/20*[1 2 3 4 -20 4 3 2 1];
                    % %                         index2=abs(conv(padarray(right,[0,4],'symmetric'),f1,'valid')<10)&(right~=im_h);
                    % %                         p=polyfit((xq(index2)-1)*num+1,right(index2),2);
                    % %                         xq1=1:im_w;
                    % %                         right=polyval(p,xq1);
                    % %                         I=flat(I,right);
                    % %                         I=I(end-200:end,:);
                    %
                    %                         [im_h, im_w] = size(I);
                    %                         w=fspecial('gaussian',[35,35],10);
                    %                         I1=imfilter(I,w);
                    %                         for ii=1:im_w
                    %                             p=max(I1(:,ii));
                    %                             I1(I1(:,ii)<p/4,ii)=0;
                    %                         end
                    %                         I=I1.*I;
                    %                         %pause(0.01);imshow(I);
                    %                         a=sum(I,2);
                    %                         left=1;right=length(a);
                    %                         count1=0;count2=0;
                    %                         [~,peak]=max(a);
                    %                         ii=1;
                    %                         while ii<length(a);
                    %                             if a(ii)<15
                    %                                 if left~=1
                    %                                 count2=count2+1;
                    %                                 if count2>5
                    %                                     right=ii-count2;
                    %                                     break;
                    %                                 end
                    %                                 else
                    %                                     count1=0;
                    %                                 end
                    %                             else
                    %                                 if left~=1
                    %                                    count2=0;
                    %                                 else
                    %                                     count1=count1+1;
                    %                                     if count1>5
                    %                                        left=ii-count1;
                    %                                        ii=peak;
                    %                                     end
                    %                                 end
                    %
                    %                             end
                    %                             ii=ii+1;
                    %                         end
                    %                         I=I(left:right,:);
                    [im_h, im_w] = size(I);
                    % make grid sampling SIFT descriptors
                    remX = mod(im_w-patchSize,gridSpacing);
                    offsetX = floor(remX/2)+1;
                    remY = mod(im_h-patchSize,gridSpacing);
                    offsetY = floor(remY/2)+1;
                    
                    [gridX,gridY] = meshgrid(offsetX:gridSpacing:im_w-patchSize+1, offsetY:gridSpacing:im_h-patchSize+1);
                    
                    fprintf('Processing %s: wid %d, hgt %d, grid size: %d x %d, %d patches\n', ...
                        frames(jj).name, im_w, im_h, size(gridX, 2), size(gridX, 1), numel(gridX));
                     patchNumberToatal=patchNumberToatal+numel(gridX);
                    % find SIFT descriptors
                    siftArr = sp_find_sift_grid(I, gridX, gridY, patchSize, 0.8);
                    [siftArr, siftlen] = sp_normalize_sift(siftArr, nrml_threshold);
                    
                    siftLens = [siftLens; siftlen];
                    
                    feaSet.feaArr = siftArr';
                    feaSet.x = gridX(:) + patchSize/2 - 0.5;
                    feaSet.y = gridY(:) + patchSize/2 - 0.5;
                    feaSet.label=database.nclass;
                    feaSet.width = im_w;
                    feaSet.height = im_h;
                    [pdir, fname] = fileparts(frames(jj).name);
                    fpath = fullfile(rt_data_dir, subname,subname1,[fname, '.mat']);
                    save(fpath, 'feaSet');
                    database.path{database.nclass,count,jj} = fpath;
                    %                         disp(num2str(database.nclass));
                    %                         disp(num2str(count));
                    %                         disp(num2str(jj));
                    %                         disp(num2str(database.path{database.nclass,count,jj}));
                end;
            end
        end
    end;
end;
patchNumberPerImage=floor(patchNumberToatal/c_num_sum);
lenStat = hist(siftLens, 100);
