function [X] = rand_sampling(database, num_smp,Test, num_img)
% sample local features for unsupervised codebook training
num_per_img=round(num_smp/num_img);
num_smp = num_per_img*num_img;
load(database.path{1,1,1});
dimFea = size(feaSet.feaArr, 1);
X = zeros(dimFea, num_smp);
cnt = 0;
[l,n]=size(Test);
c=0;
for i=1:l
    for j=1:n-1
        for m=1:database.imnum(i,Test(i,j))
            c=c+1;
            fpath = database.path{i,Test(i,j),m};
            load(fpath);
            num_fea = size(feaSet.feaArr, 2);
            rndidx = randperm(num_fea);
            X(:, cnt+1:cnt+num_per_img) = feaSet.feaArr(:, rndidx(1:num_per_img));
            cnt = cnt+num_per_img;
        end
    end
end
