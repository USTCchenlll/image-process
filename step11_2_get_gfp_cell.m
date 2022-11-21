%生成gfp细胞空间矩阵

clc
clear
%gfp_stage_list = load('D:\CL\jisui\jisui_image_analysis\20220806_img_flow\intermediate_data\gfp_imformation\have_gfp_stage.txt');
load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\2_stage_data\coordinate_stage1003.mat');
load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\2_stage_data\all_slice_file_list.mat');
dim = [2048,2048];


all_gfp_spot = zeros(0,13);
stage_num  = 949;
num_gfp = 0;
load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\all_spot_unfilter1008.mat');
all_spot = all_spot_unfilter1008; clear all_spot_unfilter1008
all_spot(all_spot(:,7)<all_spot(:,6),:) = [];
all_spot(all_spot(:,6)<all_spot(:,5),:) = [];


%%
% 读取视野坐标 stage位置
% x y zstack focus round code1 code2 stage cell cellx celly slice area

for stage_i = 1:stage_num
    disp(stage_i)
    clear xy_location temp_spot
    xyz = find(all_slice_file_list == stage_i);
    [xs,ys,zs] = ind2sub([8,10,16],xyz);
    % 读取gfp 图像
    raw_gfp_img = double(imread(['Z:\CL\20221003_jisui_no18\z_projection\cell_img\6_ex_gfp_cell\gfp_cell_s' num2str(stage_i) '_img.tif']));
    % 处理gfp图像，确定gfp细胞区域和坐标
    if isempty(find(raw_gfp_img==9999, 1)) continue; end
    gfp_img = zeros(dim);
    gfp_img(raw_gfp_img==9999) = 1; 
    % imshow(raw_gfp_img,[1,9999])   
    L = bwlabel(gfp_img,8);      
    L(L>0) = L(L>0) + num_gfp;
    num_gfp = num_gfp + length(unique(L)) - 1;
    tbl = tabulate(L(:)); 
    num_L = unique(L);
    num_L(1,:) = [];
    
    for i1 = 1:length(num_L)        
       [y,x] = find(L==num_L(i1,1));
       xy_location(i1,1) = num_L(i1,1); xy_location(i1,2) = mean(x); xy_location(i1,3) = mean(y);
    end
    % imshow(gfp_img)
	% 读取过滤后的spot(分别过滤)
	temp_spot = all_spot(all_spot(:,8)==stage_i,:);    
    for i1 = 1:length(temp_spot(:,1))
    temp_spot(i1,9) = L(round(temp_spot(i1,2)),round(temp_spot(i1,1)));
    if L(round(temp_spot(i1,2)),round(temp_spot(i1,1)))>0
    temp_spot(i1,13) = tbl(tbl(:,1)==L(round(temp_spot(i1,2)),round(temp_spot(i1,1))),2);
    end
    if temp_spot(i1,9)>0
        temp_spot(i1,10) = xy_location(xy_location(:,1)==temp_spot(i1,9),2)+coordinate_stage{all_slice_file_list == stage_i}(2);
        temp_spot(i1,11) = xy_location(xy_location(:,1)==temp_spot(i1,9),3)+coordinate_stage{all_slice_file_list == stage_i}(1);
    end
    
    end
    temp_spot(temp_spot(:,9)==0,:) = [];
    temp_spot(:,12) = zs;
    
    all_gfp_spot  = [all_gfp_spot;temp_spot];
end
save('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\5_gfp_data\all_gfp_spot1008_dilate.mat','all_gfp_spot')
save('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\5_gfp_data\all_gfp_spot1008_dilate.txt', 'all_gfp_spot', '-ascii', '-tabs')

%%

% 生成单细胞矩阵 
gene_code = importdata('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\0_code_data\code_gene_20220910.xlsx');
gene_code_name = gene_code.textdata.sum_list;
gene_code = gene_code.data.sum_list;

gene_panel_gfp = zeros(max(all_gfp_spot(:,9)),length(gene_code_name)+3);

for i1 = 1:max(all_gfp_spot(:,9))
    temp_gene_spot = all_gfp_spot(all_gfp_spot(:,9)==i1,:);
    if isempty(temp_gene_spot) continue; end
    gene_panel_gfp(i1,length(gene_code_name)+1) = temp_gene_spot(1,10); %x
    gene_panel_gfp(i1,length(gene_code_name)+2) = temp_gene_spot(1,11); %x
    gene_panel_gfp(i1,length(gene_code_name)+3) = temp_gene_spot(1,12); %x
    gene_panel_gfp(i1,length(gene_code_name)+4) = temp_gene_spot(1,13); %x
    for gene_i = 1:length(gene_code_name)
	code_1 =  find(temp_gene_spot(:,6) == gene_code(gene_i,1));
    code_2 =  find(temp_gene_spot(:,7) == gene_code(gene_i,2));
    if ~isempty(intersect(code_1,code_2))
    gene_panel_gfp(i1,gene_i) = length(intersect(code_1,code_2));   
    end    
    
        
    end

end

% clear gene_panel_gfp
xlswrite('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\5_gfp_data\gene_panel_gfp1008.xlsx',gene_panel_gfp);














































