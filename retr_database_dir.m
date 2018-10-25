function [database] = retr_database_dir(rt_data_dir)
%=========================================================================
% inputs
% rt_data_dir   -the rootpath for the database. e.g. '../data/caltech101'
% outputs
% database      -a tructure of the dir
%                   .path   pathes for each image file
%                   .label  label for each image file
% written by Jianchao Yang
% Mar. 2009, IFP, UIUC
%=========================================================================

fprintf('dir the database...');
subfolders = dir(rt_data_dir);

database = [];

database.imnum = 0; % total image number of the database
database.cname = {}; % name of each class
database.path = {}; % contain the pathes for each image of each class
database.nclass = 0;

for ii = 1:length(subfolders),
    subname = subfolders(ii).name;
    
    if ~strcmp(subname, '.DS_Store')&~strcmp(subname, '.') & ~strcmp(subname, '..'),
        database.nclass = database.nclass + 1;
        database.cname{database.nclass} = subname;       
        subfolders1 = dir([rt_data_dir,'/',subname]);
        count=0;
        for kk=1:length(subfolders1),
                subname1 = subfolders1(kk).name;               
                if ~strcmp(subname1, '.DS_Store')&~strcmp(subname1, '.') & ~strcmp(subname1, '..'),
                    count=count+1;
                    frames = dir(fullfile(rt_data_dir, subname,subname1, '*.mat'));        
                    c_num = length(frames);           
                    database.imnum(database.nclass,count) = c_num;                       
                    for jj = 1:c_num,
                        [pdir, fname] = fileparts(frames(jj).name);                        
                        fpath = fullfile(rt_data_dir, subname,subname1,[fname, '.mat']);            
                        database.path{database.nclass,count,jj} = fpath;
                    end;
                end
        end
    end;
end;
disp('done!');