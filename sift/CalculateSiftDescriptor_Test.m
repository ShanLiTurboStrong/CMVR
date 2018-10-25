function [feaSet,lenStat] = CalculateSiftDescriptor_Test(I, gridSpacing, patchSize, maxImSize, nrml_threshold)
if ndims(I) == 3,
    I = im2double(rgb2gray(I));
else
    I(I>=255)=0;
    I = im2double(I);
end;          






% [im_h,im_w] = size(I);
% w=fspecial('gaussian');
% I=imfilter(I,w);
% gridnum=128;
% num=floor(im_w/gridnum);
% right=im_h*ones(1,gridnum);
% for ii=1:gridnum
%     fea=sum(I(:,(ii-1)*num+1:ii*num),2);
%     [p]=max(fea);
%     for i=im_h:-1:1
%         if fea(i)>p/2
%             right(ii)=i;
%             break;
%         end
%     end
% end
% xq=1:gridnum;
% f1=1/20*[1 2 3 4 -20 4 3 2 1];
% index2=abs(conv(padarray(right,[0,4],'symmetric'),f1,'valid')<10)&(right~=im_h);            
% p=polyfit((xq(index2)-1)*num+1,right(index2),2);
% xq1=1:im_w;
% right=polyval(p,xq1);
% I=flat(I,right);
% I=I(end-200:end,:);



% %   BM3D       Read a grayscale image and scale its intensities in range [0,1]
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






[im_h, im_w] = size(I);




% w=fspecial('gaussian',[35,35],10);
% I1=imfilter(I,w);
% for ii=1:im_w
%     p=max(I1(:,ii));
%     I1(I1(:,ii)<p/4,ii)=0;
% end
% I=I1.*I;
% a=sum(I,2);
% left=1;right=length(a);
% count1=0;count2=0;
% [~,peak]=max(a);
% ii=1;
% while ii<length(a);
%     if a(ii)<15 
%         if left~=1
%         count2=count2+1;
%         if count2>5
%             right=ii-count2;
%             break;
%         end
%         else
%             count1=0;
%         end
%     else
%         if left~=1
%            count2=0;                    
%         else
%             count1=count1+1;
%             if count1>5
%                left=ii-count1;
%                ii=peak;
%             end
%         end
% 
%     end
%     ii=ii+1;
% end
% I=I(left:right,:);


[im_h, im_w] = size(I);

% make grid sampling SIFT descriptors
remX = mod(im_w-patchSize,gridSpacing);
offsetX = floor(remX/2)+1;
remY = mod(im_h-patchSize,gridSpacing);
offsetY = floor(remY/2)+1;
[gridX,gridY] = meshgrid(offsetX:gridSpacing:im_w-patchSize+1, offsetY:gridSpacing:im_h-patchSize+1);
% find SIFT descriptors
siftArr = sp_find_sift_grid(I, gridX, gridY, patchSize, 0.8);
[siftArr, siftlen] = sp_normalize_sift(siftArr, nrml_threshold);
feaSet.feaArr = siftArr';
feaSet.x = gridX(:) + patchSize/2 - 0.5;
feaSet.y = gridY(:) + patchSize/2 - 0.5;
feaSet.width = im_w;
feaSet.height = im_h;
lenStat = hist(siftlen, 100);
end
