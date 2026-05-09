% 读取match数据，加上sample-id，找到对应的x和y，找到对应的probe，合并umi
clc
clear
% load data
root_path = 'D:\CL\AFISH1\afish_data_20251009\';
addpath(genpath(pwd))
sample_list = {'1a','1b','1s','2a','2b','2s','3a','3b','3s'};

% load stage umi
stage = readcell('D:\CL\AFISH1\code_20251009\input\stage_umi.xlsx','Sheet','Sheet2');
stage_x = regexprep(stage, {'A','T','C','G'}, {'1','2','3','4'});
stage_x = str2double(stage_x(:,1));
stage_y = cellfun(@(s) seqrcomplement(s), stage, 'UniformOutput', false);
stage_y = regexprep(stage_y, {'A','T','C','G'}, {'1','2','3','4'});
stage_y = str2double(stage_y(:,1));
% load probe info
probe_info = readcell('D:\CL\AFISH1\code_20251009\input\slices_1535g_info.xlsx','Sheet','Sheet2');
probe_info = cell2mat(probe_info(:,4:5));

for i1 = 1:length(sample_list)
    disp(i1)
    inputname = ['D:/CL/AFISH1/afish_data_20251009/match_data/' sample_list{i1} '_part_full_info_simple_exact.txt'];
    t_data = readmatrix( inputname ,'Delimiter','_');
    t_data = [i1*ones(length(t_data(:,1)),1) t_data];    
    % 1-6[sample f/r x y umi1 umi2 ]
    % 7-10[gene region left/right wrong ] 
    % 11-14[gene region left/right wrong ]    
    % get stage umi
    tx = t_data(:,3:5);
    for i2 = 1:length(stage)
        tx(tx==stage_x(i2)) = i2;     
    end
    ty = t_data(:,6:8);
    for i2 = 1:length(stage)
        ty(ty==stage_y(i2)) = i2;     
    end   
    
    % tx: N x 3 ; 情况1：a 跟 b 或 a 跟 c 相等 → 取 a; 情况2：前面没命中，并且 b == c → 取 b
    a = tx(:,1);b = tx(:,2);c = tx(:,3);
    tx1 = zeros(size(tx,1),1);
    m1 = (a == b) | (a == c);
    tx1(m1) = a(m1);
    m2 = ~m1 & (b == c);
    tx1(m2) = b(m2);
    
    % ty: N x 3
    a = ty(:,1);b = ty(:,2);c = ty(:,3);
    ty1 = zeros(size(ty,1),1);
    m1 = (a == b) | (a == c);
    ty1(m1) = a(m1);
    m2 = ~m1 & (b == c);
    ty1(m2) = b(m2);

    % 你后面的组装
    t_data = [t_data(:,1:2), tx1, ty1, t_data(:,9:16)];
    t_data(:,5) = t_data(:,5)*1000000+t_data(:,6);
    
    % 假设这两列都是正整数
    colA = t_data(:,11);   % 对应 probe_info(i2,1)
    colB = t_data(:,10);   % 对应 probe_info(i2,2)
    maxB = max(colB);      % 用这列来做进位
    % 给 t_data 每一行生成唯一key
    t_key = colA * (maxB+1) + colB;   % N×1
    % 给 probe_info 每一行也生成同样的key
    p_key = probe_info(:,1) * (maxB+1) + probe_info(:,2);   % M×1
    % 我们要一个 size = max(t_key) 的查表向量
    lookup = zeros(max(t_key), 1);    % 注意：如果 key 很稀疏、最大值很大，下面有别的写法
    % 把 probe_info 的 key 映射成索引 i2
    % probe_info 的行号就是你想写进去的 i2
    lookup(p_key) = 1:numel(p_key);
    % 最后：直接用查表的方式给 t_data 第6列赋值
    t_data(:,6) = lookup(t_key);
  
    h5create([root_path 'mtx_data\' sample_list{i1} '_mtx_simple_data.h5'], '/dataset1', size(t_data));   % 创建数据集
    h5write([root_path 'mtx_data\' sample_list{i1} '_mtx_simple_data.h5'], '/dataset1', t_data);          % 写入数据
end

%% filter data
clc
clear
root_path = 'D:\CL\AFISH1\afish_data_20251009\';
addpath(genpath(pwd))
sample_list = {'1a','1b','1s','2a','2b','2s','3a','3b','3s'};

% make data simple
data = zeros(0,8);
for i1 = 1:length(sample_list)
    disp(i1)
    t_data = h5read([root_path 'mtx_data\' sample_list{i1} '_mtx_simple_data.h5'],'/dataset1');    
    t_data = t_data((t_data(:,7)==t_data(:,10))&(t_data(:,8)==t_data(:,11))&(t_data(:,9)==2)&(t_data(:,12)==1)==1,:);   
    t_data = t_data(t_data(:,3)<9&t_data(:,4)<13,:);
    t_data = t_data(:,1:8);
    data = [data; t_data];
end

h5create([root_path 'mtx_data\full_filtered_mtx_data.h5'], '/dataset1', size(data));   % 创建数据集
h5write([root_path 'mtx_data\full_filtered_mtx_data.h5'], '/dataset1', data);          % 写入数据

%%
clc
clear
root_path = 'D:\CL\AFISH1\afish_data_20251009\';
addpath(genpath(root_path))
fdata = h5read([root_path 'mtx_data\full_filtered_mtx_data.h5'],'/dataset1');
x = fdata(:,1);  % 提取第二列
x_new = fdata(:,3);
x_new(x >= 4 & x <= 6) = x_new(x >= 4 & x <= 6) + 8;
x_new(x >= 7 & x <= 9) = x_new(x >= 7 & x <= 9) + 16;
fdata(:,3) = x_new;

% 假设 fdata 的列意义： [~, x, y, umi, grp]，尺寸 N×5
x   = fdata(:,3);
y   = fdata(:,4);
umi = fdata(:,5);
grp = fdata(:,6);              % i1 分组（假定取值为 1..G）

Xmax = 24;  Ymax = 12;  G = 5476;   % 你的网格与分组上限

% 过滤越界/非法（可选）
ok = x>=1 & x<=Xmax & y>=1 & y<=Ymax & grp>=1 & grp<=G & ~isnan(umi);
x = x(ok);  y = y(ok);  grp = grp(ok);  umi = umi(ok);
% 1) 每个 (x,y,grp) 的条数
map_probe_full = accumarray([x y grp], 1, [Xmax Ymax G], @sum, 0);
% 2) 每个 (x,y,grp) 的 UMI 去重计数
map_probe_umi  = accumarray([x y grp], umi, [Xmax Ymax G], ...
                      @(v) numel(unique(v(~isnan(v)))), 0);
save('map_probe_umi.mat','map_probe_umi')
save('map_probe_full.mat','map_probe_full')                  
 %% 将探针地图转换为基因地图                 

clc
clear
root_path = 'D:\CL\AFISH1\afish_data_20251009\';
addpath(genpath(root_path))
% load probe info
probe_info = readcell('D:\CL\AFISH1\code_20251009\input\slices_1535g_info.xlsx','Sheet','Sheet2');
probe_info = cell2mat(probe_info(:,4:5));
gene_info = readcell('D:\CL\AFISH1\code_20251009\input\slices_1535g_info.xlsx','Sheet','gene');

load map_probe_umi
load map_probe_full

g = max(probe_info(:,2));

for i1 = 1:g
    t = map_probe_umi(:,:,probe_info(:,2)==i1); 
    map_gene_umi(:,:,i1) = mean(t,3); 
end

map_gene_umi_have_counts = map_gene_umi(:,:,squeeze(sum(sum(map_gene_umi)))>0);
gene_have_counts = gene_info(squeeze(sum(sum(map_gene_umi)))>0,:);

save('full_info_with_gapdh.mat')









