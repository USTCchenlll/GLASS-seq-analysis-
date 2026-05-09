% 通过flash拼接R1和R2，输出完整序列的fastq
% 增加了去除接头的功能 20251009
clc
clear

addpath(genpath(pwd))
raw_path = 'D:\CL\AFISH1\afish_data_20251009\raw_seq_data\';
allFolders = dir(raw_path);  % 获取路径下的所有文件和文件夹
allFolders = allFolders([allFolders.isdir]);  % 筛选出文件夹
allFolders = allFolders(~ismember({allFolders.name}, {'.', '..'}));  % 排除"."和".."这两个特殊文件夹
allFolders = struct2cell(allFolders).';

% 显示文件夹列表
for i1 = 12:length(allFolders)
    if i1==6||i1==9||i1==11||i1==12;else;continue;end
    disp(i1)
    k1 = strfind(allFolders{i1,1},'_');
    allFolders{i1,7}=allFolders{i1,1}(k1+1:end);
    k2 = strfind(allFolders{i1,1},'-');
    allFolders{i1,8}=allFolders{i1,1}(k2(2)+1:end);
    
    r1 = [raw_path '/' allFolders{i1,1} '/' allFolders{i1,7} '_combined_R1.fastq.gz'];
    r2 = [raw_path '/' allFolders{i1,1} '/' allFolders{i1,7} '_combined_R2.fastq.gz'];
    r1_out = [raw_path '/' allFolders{i1,1} '/' allFolders{i1,7} '_cutadapt_R1.fastq.gz' ];
    r2_out = [raw_path '/' allFolders{i1,1} '/' allFolders{i1,7} '_cutadapt_R2.fastq.gz' ]; 
    
    cmd = ['call C:\ProgramData\anaconda3\condabin\conda.bat activate afish_seq && cutadapt ', ...
           ' -j 16 -a AGATCGGAAGAGC -A AGATCGGAAGAGC ', ...
           '-o ' r1_out ' -p ' r2_out ' ', r1 ' ' r2 '' ];
    status = system(cmd);
    if status ~= 0
        error('cutadapt 运行失败，请确认是否正确安装并配置环境变量');
    end

    system(['flash ' raw_path '/' allFolders{i1,1} '/' allFolders{i1,7} '_cutadapt_R1.fastq.gz' ...
        ' ' raw_path '/' allFolders{i1,1} '/' allFolders{i1,7} '_cutadapt_R2.fastq.gz -m 40 -M 200 ' ...
        '-o ' allFolders{i1,8}])
end


