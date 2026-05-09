clc
clear
addpath(genpath(pwd))
% load data
root_path = 'D:\CL\AFISH1\afish_data_20250916\mtx_data\';
x = [1,2];
y = [1,2,3,4,5,6,7,8];

ly = length(y);lx = length(x);   

region_list = {'CB-1','CB-2','CT-1','CT-2','STR-1','STR-2','HIP','HB'};

gene_list = readcell('slices_1535g_info.xlsx','Sheet','gene');
lg = length(gene_list);

for i2 = 1:lx
    clear mtx_data
    mtx_data = zeros(ly,lg);
    norm_mtx_data = zeros(ly,lg);
    
    umi_data = zeros(ly,lg);
    norm_umi_data = zeros(ly,lg); 
    sum_k = zeros(0,5);
    for i1 = 1:ly
        disp(['sample_' num2str(i2) '_region_' num2str(i1)])
        a = h5read([root_path 'x_' num2str(x(i2)) '_y_' num2str(y(i1)) '_mtx_data.h5' ],'/dataset1');
        a = a(a(:,7)==a(:,11)&a(:,8)==a(:,12)&a(:,9)==2&a(:,13)==1&a(:,10)<3&a(:,14)<3,:);
        
        for i3 = 1:lg
            t = a(a(:,7)==i3,:);
            if isempty(t);continue;end
            k = sortrows(tabulate(t(:,8)), 2);
            k(:,5)=i3;
            mtx_data(i1,i3)= mean(k(:,2));         
            for i4 = 1:length(k(:,1))
                tu = t(t(:,8)==k(i4,1),:);
                if isempty(tu);continue;end               
                % 假设 M 是 10000×10 的 matrix，Col5 和 Col6 是 UMI
                umi1 = tu(:,5);
                umi2 = tu(:,6);
                % 拼接成组合 UMI，转为字符串表示
                joint_umi = strcat(string(umi1), "_", string(umi2));
                % 找出唯一组合及其所在行（首次出现的行）
                [~, ia] = unique(joint_umi, 'stable');  % 保留第一次出现的顺序
                k(i4,4)=length(ia);               
            end
            sum_k=[sum_k;k];
            umi_data(i1,i3)=mean(k(:,4));           
        end
        norm_mtx_data(i1,:) = mtx_data(i1,:)/sum(mtx_data(i1,:));
        norm_umi_data(i1,:) = umi_data(i1,:)/sum(umi_data(i1,:));
               
    end
    
    mtx = [[{''},gene_list.'];[region_list.',num2cell(mtx_data)]];
    writecell(mtx,['ad_' num2str(i2) '_e.xlsx'])    
    mtx = [[{''},gene_list.'];[region_list.',num2cell(norm_mtx_data)]];
    writecell(mtx,['ad_' num2str(i2) '_ne.xlsx'])    
    mtx = [[{''},gene_list.'];[region_list.',num2cell(umi_data)]];
    writecell(mtx,['ad_' num2str(i2) '_u.xlsx'])    
    mtx = [[{''},gene_list.'];[region_list.',num2cell(norm_umi_data)]];
    writecell(mtx,['ad_' num2str(i2) '_nu.xlsx'])      
    writematrix(sum_k,['ad_' num2str(i2) '_probe_umi.xlsx'])         
end







