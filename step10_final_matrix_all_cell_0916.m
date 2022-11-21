clear
clc

stage_num = 949;
gene_code = importdata('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\0_code_data\code_gene_20220910.xlsx');
gene_code_name = gene_code.textdata.sum_list;
gene_code = gene_code.data.sum_list;
l_gene_code = length(gene_code_name); 


load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\2_stage_data\all_slice_file_list.mat');

load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\3_cell_data\cell_list1008_dilate.mat');
cell_list(cell_list(:,6)==1,:) = [];

load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\all_spot_1008.mat');
all_spot = all_spot_1008; clear all_spot_1008
all_spot(all_spot(:,7)<all_spot(:,6),:) = [];

for slice_i = 1:length(all_slice_file_list(1,1,:))
    tic
    temp_slice_stage = all_slice_file_list(:,:,slice_i);
    stage_list = sort(unique( temp_slice_stage));
    stage_list(1,:) = [];
    
    temp_slice_spot  = zeros(0,9);
    temp_cell_list = cell_list(cell_list(:,7)==slice_i,:);
    all_check = zeros(length(all_spot(:,1)),1);
    for i1 = 1:length(stage_list(:,1))
        stage_i = stage_list(i1,1);
        all_check  = all_check+(all_spot(:,8)==stage_i);
    end
    temp_spot_list = all_spot(logical(all_check),:);
    
    % hongguan check
   hongguan(slice_i,1) = length(temp_spot_list);
   hongguan(slice_i,2) = length(temp_cell_list); 
    
    slice_gene_panel = zeros(length(temp_cell_list(:,1)),l_gene_code+3);
    slice_gene_panel(:,l_gene_code+1) = temp_cell_list(:,3); %x
    slice_gene_panel(:,l_gene_code+2) = temp_cell_list(:,4); %y
    slice_gene_panel(:,l_gene_code+3) = temp_cell_list(:,8); %area
    
    
for cell_i = 1:length(temp_cell_list)
    if mod(cell_i,100)==0  disp([ num2str(cell_i) '/' num2str(length(temp_cell_list))  ]); end
    temp_stage =  temp_spot_list(:,8) == temp_cell_list( cell_i,2); %stage
    temp_cell =  temp_spot_list(:,9) == temp_cell_list( cell_i,1); %cell

    if sum(temp_stage.*temp_cell)==0 continue; end
    temp_spot = temp_spot_list(logical((temp_stage.*temp_cell)),:);
    
    
	for code_i = 1:length(gene_code)
    code_1 =  temp_spot(:,6) == gene_code(code_i,1);
    code_2 =  temp_spot(:,7) == gene_code(code_i,2);    
    slice_gene_panel(cell_i,code_i) = sum(code_1.*code_2);    
	end    
     
end


save(['D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\4_final_data\gene_panel_no1_slice' num2str(slice_i) '_1008_dilate.mat'],'slice_gene_panel');
xlswrite(['D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\4_final_data\gene_panel_no1_slice' num2str(slice_i) '_1008_dilate.xlsx'],slice_gene_panel)

xlswrite(['D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\4_final_data\cell_list_no1_slice' num2str(slice_i) '_1008_dilate.xlsx'],temp_cell_list)

toc
   
end

%%


% l_gene_code = length(gene_code_name); 

% gene_panel = zeros(length(cell_list(:,1)),l_gene_code+3);
% gene_panel(:,l_gene_code+1) = cell_list(:,3); %x
% gene_panel(:,l_gene_code+2) = cell_list(:,4); %y
% gene_panel(:,l_gene_code+3) = cell_list(:,8); %area
% tic
% for cell_i = 1:5000
%     if mod(cell_i,100)==0  disp([ num2str(cell_i) '/' num2str(length(cell_list))  ]); end
%     temp_stage =  all_spot(:,8) == cell_list( cell_i,2); %stage
%     temp_cell =  all_spot(:,9) == cell_list( cell_i,1); %cell
% 
%     if sum(temp_stage.*temp_cell)==0 continue; end
%     temp_spot = all_spot(logical((temp_stage.*temp_cell)),:);
%     
%     
% 	for code_i = 1:length(gene_code)
%     code_1 =  temp_spot(:,6) == gene_code(code_i,1);
%     code_2 =  temp_spot(:,7) == gene_code(code_i,2);    
%     gene_panel(cell_i,code_i) = sum(code_1.*code_2);    
% 	end    
%      
% end
% toc
% gene_panel_filterd = gene_panel(cell_list(:,6)==0,:);
% save('D:\CL\jisui\jisui_image_analysis\20220830_img_flow\intermediate_data\4_final_data\gene_panel_no1_0907.mat','gene_panel');
% save('D:\CL\jisui\jisui_image_analysis\20220830_img_flow\intermediate_data\4_final_data\gene_panel_filterd_no1_0907.mat','gene_panel_filterd');
% 
% a = cell_list(cell_list(:,6)==0,[1,5,7]);% 相对序号 绝对序号 stage序号
% xlswrite('D:\CL\jisui\jisui_image_analysis\20220830_img_flow\intermediate_data\final_data\gene_panel_filterd_no1_0907.xlsx',gene_panel_filterd)
% xlswrite('D:\CL\jisui\jisui_image_analysis\20220830_img_flow\intermediate_data\final_data\cell_list_no1_0907.xlsx',a)
% 








