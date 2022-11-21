%生成gfp细胞空间矩阵

clc
clear

dim = [2048,2048];
gfp_intensity = 12000;
small_size = 5;
SE1 = strel('disk',small_size);

all_gfp_spot = zeros(0,13);
stage_num  = 949;
num_gfp = 0;


%%
% 读取视野坐标 stage位置
% x y zstack focus round code1 code2 stage cell cellx celly slice area

for stage_i = 1:stage_num
    disp(stage_i)
    clear xy_location temp_spot
    % 读取gfp 图像
    raw_gfp_img = double(imread(['Z:\CL\20221003_jisui_no18\z_projection\round1\max_20221003_jisui18_round11_GFP-only-3_s' num2str(stage_i) '.ome.tif' ]));
    raw_cell_img = double(imread(['Z:\CL\20221003_jisui_no18\z_projection\cell_img\3_drift_z_cell\max_20221003_jisui18_no1_round11_DAPI_s' num2str(stage_i) '.ome.tif_stardist.tiff'])); 
    if sum(sum(raw_gfp_img>gfp_intensity))<500 
        imwrite(uint16(raw_cell_img),['Z:\CL\20221003_jisui_no18\z_projection\cell_img\5_raw_gfp_cell\gfp_cell_s' num2str(stage_i) '_img.tif'])
        continue; 
%     else
%         continue;
    end
    
    raw_cell_img = imtranslate(raw_cell_img,[-16 1]);
    % 处理gfp图像，确定gfp细胞区域和坐标
    temp_cell_list = sort(unique(raw_cell_img));temp_cell_list(1,:)=[];
    temp_cell_num = length(unique(raw_cell_img))-1;
    for cell_i = 1:temp_cell_num
        if length(intersect(find(raw_cell_img==temp_cell_list(cell_i,1)),find(raw_gfp_img>gfp_intensity)))/length(find(raw_cell_img==temp_cell_list(cell_i,1)))>0.5
        raw_cell_img(raw_cell_img==temp_cell_list(cell_i,1)) = 9999;   
        end       
    end
%     imshow(raw_cell_img,[0,2000])
%     hold on 
%     imshow(raw_gfp_img,[15000,20000])
    gfp_img = zeros(dim);
    gfp_img(raw_gfp_img>gfp_intensity) = 1;
    gfp_img = imerode(gfp_img,SE1) ;    
    gfp_img = imdilate(gfp_img ,SE1);
    raw_cell_img(gfp_img==1) = 9999;   
    % imshow(raw_cell_img,[1,2000])
    imwrite(uint16(raw_cell_img),['Z:\CL\20221003_jisui_no18\z_projection\cell_img\5_raw_gfp_cell\gfp_cell_s' num2str(stage_i) '_img.tif'])

end
