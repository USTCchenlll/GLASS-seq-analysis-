function encoding_sequence_list = region2seq(es_generated,file_name)
%将潜在探针序列，整合成fasta格式的txt文档
%输出结果举例：
% >|125|ENSMUST00000085374.7|ENSMUSG00000070570.6|OTTMUSG00000058649.1|OTTMUST00000144033.1|Slc17a7-201|Slc17a7|2915|protein_coding|
% CACAGCCACCATGGAGTTCCGGCAGGAGGAGT

%2022.11.22 CL
% update 20240415 CL
encoding_sequence_list = fopen(file_name,'a+');
for i = 1:length(es_generated(:,1))
    fprintf(encoding_sequence_list,'>');
    fprintf(encoding_sequence_list,es_generated{i,2});
    fprintf(encoding_sequence_list,'\n');
    fprintf(encoding_sequence_list,es_generated{i,1});
    fprintf(encoding_sequence_list,'\n');
end
    fclose all;
    disp('sequence已导出'); 