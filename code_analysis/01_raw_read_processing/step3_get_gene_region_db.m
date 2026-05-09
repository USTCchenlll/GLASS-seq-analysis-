% 读取片段序列，分割出坐标，UMI和目标片段
% 重新存为数据，进行blast

clc
clear
% load data
root_path = 'D:\CL\AFISH1\afish_data_20251009\';
addpath(genpath(pwd))
sample_list = {'1a'};

% 读取基因序列信息
raw_region_gene_list = importdata("D:\CL\AFISH1\code_20251009\input/slices_1535g_info.xlsx");
raw_region_list = raw_region_gene_list.textdata.Sheet1;
gene_list = raw_region_gene_list.textdata.gene;

for i1 = 1:length(raw_region_list)
    r = strsplit(raw_region_list{i1,1},'_');
    % 将整数转换为三位的字符串
    rr = [sprintf('%03d', find(strcmp(gene_list,r{2})==1)) '_' sprintf('%03d', str2double(r{3})) '_'];
    region_list(2*i1-1,:) = [ raw_region_list{i1,2}  {[rr '1']}];
    region_list(2*i1,:) = [ raw_region_list{i1,3}  {[rr '2']}];
end
clear raw_region_list raw_region_gene_list
sequence_list = region2seq(region_list,'1535_region_list.txt');
