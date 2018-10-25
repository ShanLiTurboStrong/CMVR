function  [averageImgFlattentime,averageWidth,averageHeight,heightVector]=FlattenImage(rt_img_dir,rt_data_dir)
disp('Flatten Input Image...');
subfolders = dir(rt_img_dir);
addpath('BM3D');
skip_threshold=70;
c_num_sum=0;
flatten_time_sum=0;
flatten_height_sum=0;
flatten_width_sum=0;
heightVector=zeros(1,465);
for nclass = 1:length(subfolders)
    subname = subfolders(nclass).name;%AMD,DME,NORMAL
    if ~strcmp(subname, '.DS_Store')&~strcmp(subname, '.') & ~strcmp(subname, '..')
        
        subfolders1 = dir([rt_img_dir,'/',subname]);%(AMD1~15)/(DME1~15)/(NORMAL1~15)
        for kk=1:length(subfolders1)
            subname1 = subfolders1(kk).name;%AMD/DME/NORMAL for specific person
            if ~strcmp(subname1, '.DS_Store')&~strcmp(subname1, '.') & ~strcmp(subname1, '..'),
                frames = dir(fullfile(rt_img_dir, subname,subname1, '*.tif'));
                c_num = length(frames);
                flattenPath = fullfile(rt_data_dir, subname,subname1);
                if ~isdir(flattenPath),
                    mkdir(flattenPath);
                end;
                for jj = 1:c_num
                    imgpath = fullfile(rt_img_dir, subname,subname1, frames(jj).name);
                    originalImage=imread(imgpath);%read the Grayscale Picture(uint8 image)
                    disp(imgpath);
                    TSTART=tic;
                    if ndims(originalImage) == 3,
                        I = im2double(rgb2gray(originalImage));
                    else
                        originalImage(originalImage>=255)=0;%fill up white margin
                        I = im2double(originalImage);
                    end;
                    [m,n]=size(I);%obtain the original image row&&column
                    I_Container=zeros(m+200,n);
                    I_Container(101:100+m,:)=I;
                    [m,n]=size(I_Container);
                    cropceil=1;%for cropping use
                    cropfloor=m;%for cropping use
                    uncroppedFinalImage=zeros(m,n);%define out-put image
                    
                    
                    I_Container=im2double(I_Container);%indispensable, because must transfer uint8 to double for functions in matlab tools
                    %          Read a grayscale image and scale its intensities in range [0,1]
                    y = I_Container;
                    % Generate the same seed used in the experimental results of [1]
                    randn('seed', 0);
                    %rng(0,'v4')
                    % Standard deviation of the noise --- corresponding to intensity
                    %  range [0,255], despite that the input was scaled in [0,1]
                    sigma = 35;%initial is 25, lasttime is 45
                    % Add the AWGN with zero mean and standard deviation 'sigma'
                    z = y + (sigma/255)*randn(size(y));
                    % Denoise 'z'. The denoised image is 'y_est', and 'NA = 1' because
                    %  the true image was not provided
                    [NA, y_est] = BM3D(1, z, sigma);
                    [level,EM] = graythresh(y_est); %using Ostu's algorithm to determine the optimal threshold
                    bw = im2bw(y_est,level);%transform to binary image
                    mf= medfilt2(bw,[35,35]);%median filter: 35*35?original is 6*6?last time is 20*20?lasttime is 35---
                    se1 = strel('disk',40);%Create morphological structuring element.lasttime is 40
                    closedI=imclose(mf,se1);%morphological closing
                    se2 = strel('disk',5);%lasttime is 5,modified for hospital-amd13005
                    openedI=imopen(closedI,se2);%morphological opening
                    %figure(figureNumber);
                    %subplot(221);
                    %imshow(I0);%show gray image
                    %subplot(222);
                    %imshow(bw);%show binary image
                    %subplot(223);
                    %imshow(mf);%show median filter image
                    %subplot(224);
                    %imshow(openedI);%show closing-opening image
                    
                    
                    %                     hf=figure(figureNumber+1);
                    %                     haxes1=axes('parent',hf,'tag','axes1','units','normalized','position',[0.1 0.1 0.8 0.8]);
                    %                     axes(haxes1);
                    %                     imshow(openedI);
                    %                     hold on;
                    %                     axis on;%!!!!!the original point of pic is left-top, but of poly is left-bottom,
                    %               the final origin-point is determined by the first thing showed on the
                    %               figure
                    
                    
                    %----------------------------find the ceil and floor of openedI to get acurate dataset--------------------------------%
                    flag=zeros(1,n);
                    cropceil=1;
                    cropfloor=m;
                    for R=1:m
                        for C=1:n
                            if openedI(R,C)==1 && R>flag(C)
                                ll=0;
                                while ll<skip_threshold
                                    ll=ll+1;
                                    if openedI(R+ll,C)==0
                                        flag(C)=R+ll;
                                        break;
                                    end
                                    if ll==skip_threshold
                                        cropceil=R;
                                        R=m;
                                        C=n;
                                    end
                                end
                            end
                            if C==n
                                break;
                            end
                        end
                        if R==m
                            break;
                        end
                    end
                    R=m;
                    for R=m:-1:1
                        for C=1:n
                            if openedI(R,C)==1
                                cropfloor=R;
                                R=1;%notice: can not use the inner loop to access the outer loop iterate varible
                                break;
                            end
                        end
                        if R==1
                            break;
                        end
                    end
                    %--------------------------------------------------------------------------------------------------------------%
                    openedI(1:cropceil,:)=0;
                    openedI(cropfloor:m,:)=0;
                    
                    %----------remove the columns within which the number of 0-pixels less than X in openedI----------------------%
                    sumedOpenedI=sum(openedI,1);
                    half_maxOpendedI=max(sumedOpenedI)*(1/3);
                    row=ones(1,n);
                    for left=1:n
                        if sumedOpenedI(left)<=skip_threshold || sumedOpenedI(left)<=half_maxOpendedI%X pixels is based on experience
                            row(left)=0;
                        end
                    end
                    %--------------------------------------------------------------------------------------------------------------%
                    
                    
                    datasetXsum=sum(row,2);
                    datasetMiddle=zeros(1,datasetXsum,'double');%calculate the least square curve fitting Y-dataset
                    datasetBottom=zeros(1,datasetXsum,'double');
                    datasetXcounter=0;
                    for C=1:n
                        if row(C)==1
                            datasetXcounter=datasetXcounter+1;
                        else
                            continue;
                        end
                        count=0;
                        number=0;
                        for R = cropceil:cropfloor
                            if openedI(R,C)==1
                                count=count+R;
                                number=number+1;
                            end
                        end
                        datasetMiddle(1,datasetXcounter)=count./number;%use"m-count./number"because the axsis
                        %in plot and image is opposite
                    end
                    datasetXcounter=0;
                    for R = 1:n
                        if row(R)==1
                            datasetXcounter=datasetXcounter+1;
                        else
                            continue;
                        end
                        count=cropfloor;
                        while openedI(count,R)~=1
                            count=count-1;
                        end
                        datasetBottom(1,datasetXcounter)=count;
                    end
                    
                    
                    dataset1=datasetMiddle;
                    dataset2=datasetBottom;
                    x=find(row);%the least square curve fitting X-dataset
                    
                    %datasetMiddle
                    a1 = polyfit(x,dataset1,1);%linear fit
                    a2 = polyfit(x,dataset1,2);%parabalo fit
                    %datasetBottom
                    a3 = polyfit(x,dataset2,1);%linear fit
                    a4 = polyfit(x,dataset2,2);%parabalo fit
                    if a2(1,1)>0%!!!!!!judge if middle-parabalo's a is
                        %negative(=>a2(1,1)>0, because the original
                        %is on left-top),then use the bottom dataset
                        a1=a3;
                        a2=a4;
                        dataset=dataset2;
                        %                         disp([num2str(ii),'->bottom']);
                    else%use middle dataset
                        %                         disp([num2str(ii),'->middle']);
                        dataset=dataset1;
                    end
                    x0=(-a2(1,2))./(2*a2(1,1));%calculate parabalo-axsis
                    z1=polyval(a1,x);
                    z2=polyval(a2,x);
                    coef1=corrcoef(z1,dataset);
                    coef2=corrcoef(z2,dataset);
                    %                         disp(coef1);
                    %                         disp(coef2);
                    
                    %now if (-b)/2a lies beyond the x-dimension of original picture
                    %or correlation coefficients :linear>=parabalo
                    %then use linear fit
                    if coef1(1,2)>=coef2(1,2) || a2(1,1)>0%*********here are linear codes********
                        %disp('linear-fitting');
                        %plot(x,z1);%draw the linear line onto the image
                        angle=atand(a1(1,1));
                        uncroppedFinalImage=imrotate(I_Container,angle,'bilinear','crop');
                        %crop the image
                        openedI1=imrotate(openedI,angle,'bilinear','crop');%for cropping use
                    else%*********here are parabalo codes********
                        %disp('parabalo-fitting');
                        %                         plot(x,z2);%draw the parabola onto the image,!!!!!the showed plot
                        %                         %       may up-towards,but the "a" is negative,cause the original point is
                        %                         %       left-top,which is determined by the pic
                        P1=4*a2(1,1)*a2(1,3);
                        P2=a2(1,2)*a2(1,2);
                        P3=4*a2(1,1);
                        if x0<=1
                            vertex=a2(1,1)+a2(1,2)+a2(1,3);
                        elseif  x0>=n
                            temp1=a2(1,1)*n;
                            temp1=temp1*n;
                            temp2=a2(1,2)*n;
                            vertex=temp1+temp2+a2(1,3);
                        else
                            vertex=(P1-P2)./P3;%calculate the vertex of parabola
                        end
                        subvector=zeros(1,n);
                        openedI1=zeros(m,n);
                        %jude the parabalo a(1,1)'s value
                        if a2(1,1) >=0%if the parabalo's up-towards
                            for C = 1:n
                                temp1=a2(1,1)*C;
                                temp1=temp1*C;
                                temp2=a2(1,2)*C;
                                subvector(1,C)=round(temp1+temp2+a2(1,3)-vertex);%calculate
                                %the bias to horizon at each x-cordinate of parabola
                                %we should gurantee that subvector is always positive
                                %however if the parabalo fit failed(like a(1,1)<0)
                                %then there will be negative number thus "I0(round(R-subvector(1,C)),C)"
                                %may come up error like "index out of bounds"
                            end
                            
                            for R = 1:m
                                for C = 1:n
                                    if(R+subvector(1,C)<=m)
                                        uncroppedFinalImage(R,C)=I_Container(R+subvector(1,C),C);%substract each column of picture
                                        openedI1(R,C)=openedI(R+subvector(1,C),C);%for cropping use
                                    end
                                end
                            end
                        else%if the parabalo's down-towards
                            for C = 1:n
                                temp1=a2(1,1)*C;
                                temp1=temp1*C;
                                temp2=a2(1,2)*C;
                                subvector(1,C)=round(vertex-temp1-temp2-a2(1,3));%calculate
                                %the bias to horizon at each x-cordinate of parabola
                            end
                            for R = m:-1:1
                                for C = 1:n
                                    if(R-subvector(1,C)>=1)
                                        uncroppedFinalImage(R,C)=I_Container(R-subvector(1,C),C);%substract each column of picture
                                        openedI1(R,C)=openedI(R-subvector(1,C),C);%for cropping use
                                    end
                                end
                            end
                        end
                    end
                    %for crop use
                    flag=zeros(1,n);
                    for R=1:m
                        for C=1:n
                            if openedI1(R,C)==1 && R>flag(C)
                                ll=0;
                                while ll<skip_threshold
                                    ll=ll+1;
                                    if openedI1(R+ll,C)==0
                                        flag(C)=R+ll;
                                        break;
                                    end
                                    if ll==skip_threshold
                                        cropceil=R;
                                        R=m;
                                        C=n;
                                    end
                                end
                            end
                            if C==n
                                break;
                            end
                        end
                        if R==m
                            break;
                        end
                    end
                    R=m;
                    for R=m:-1:1
                        for C=1:n
                            if openedI1(R,C)==1
                                cropfloor=R;
                                R=1;%notice: can not use the inner loop to access the outer loop iterate varible
                                break;
                            end
                        end
                        if R==1
                            break;
                        end
                    end
                    finalImage=imcrop(uncroppedFinalImage,[1 cropceil (n-1) (cropfloor-cropceil)]);
                    [withTmp,hightTmp]=size(finalImage);
                    if c_num==168
                        heightVector(1,jj)=hightTmp;
                    else
                        heightVector(1,jj+168)=hightTmp;
                    end
                    flatten_height_sum=flatten_height_sum+withTmp;
                    flatten_width_sum=flatten_width_sum+hightTmp;
                    flatten_time_sum=flatten_time_sum+toc(TSTART);
                    imwrite(finalImage,[flattenPath,'/','aligned_',frames(jj).name]);
                    
                end
            end
        end
    end
end
averageImgFlattentime=flatten_time_sum/c_num_sum;
averageHeight=flatten_height_sum/c_num_sum;
averageWidth=flatten_width_sum/c_num_sum;
disp('Flatten finished!');
end
