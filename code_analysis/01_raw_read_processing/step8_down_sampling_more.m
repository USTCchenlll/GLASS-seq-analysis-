%% filter data
clc
clear
root_path = 'E:\CL\SLICES\AFISH1\afish_data_20251009\';
addpath(genpath(pwd))
sample_list = {'1a','1b','1s','2a','2b','2s','3a','3b','3s'};
raw_data = h5read([root_path 'mtx_data\full_mtx_data.h5'],'/dataset1');    

% load probe info
probe_info = readcell('E:\CL\SLICES\AFISH1\code_20251009\input\slices_1535g_info.xlsx','Sheet','Sheet2');
probe_info = cell2mat(probe_info(:,4:5));
gene_info = readcell('E:\CL\SLICES\AFISH1\code_20251009\input\slices_1535g_info.xlsx','Sheet','gene');

for i00 = 1:1%20

for i1 = 1:1%10:10:180
    disp(i1)
    fdata = zeros(0,0);    
    for i2 = 1:9
        k = raw_data(raw_data(:,1)==i2,:);
        cy = round(length(k(:,1))*i1/180);      
        p = randperm(length(k(:,1)));
        selected_data = k(p(1:cy), :);
        fdata = [fdata;selected_data];
    end
    
    fdata = fdata((fdata(:,7)==fdata(:,10))&(fdata(:,8)==fdata(:,11))&(fdata(:,9)==2)&(fdata(:,12)==1)==1,:);   
    fdata = fdata(fdata(:,3)<9&fdata(:,4)<13,:);
    fdata = fdata(:,1:8);
    
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
    save(['down_sampling_more\' num2str(i00) '_'  num2str(i1) '_map_probe_umi.mat'],'map_probe_umi')
    save(['down_sampling_more\' num2str(i00) '_'  num2str(i1) '_map_probe_full.mat'],'map_probe_full')       

    g = max(probe_info(:,2));
    for i2 = 1:g
        t = map_probe_umi(:,:,probe_info(:,2)==i2); 
        map_gene_umi(:,:,i2) = mean(t,3); 
    end
    save([ 'down_sampling_more\' num2str(i00) '_' num2str(i1) '_map_gene_umi.mat'],'map_gene_umi')

end


end
















