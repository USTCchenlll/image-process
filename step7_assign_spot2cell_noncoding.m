% 对识别出的点信号进行细胞的分配 for noncoding 561 r12345

clc
clear
%%
% x y zstack unfocus round code1 code2 stage cell_i

stage_num = 949;

loc_all = zeros(0,8);

for stage_i = 1:stage_num
    disp(num2str(stage_i))
	
    loc = load(['D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\noncoding\noncoding_raw_localizations_s' num2str(stage_i) '.txt'],'-ascii') ;

    if ~isempty(loc)
    loc(:,6) = loc(:,5)+10;
    loc(:,7) = loc(:,5)+10;
    loc_filterd =  loc;
    loc_filterd(:,8) = stage_i;
     cell_img = double(imread(['Z:\CL\20221003_jisui_no18\z_projection\cell_img\4_ex_z_cell\max_20221003_jisui18_no1_round11_DAPI_s' num2str(stage_i)  '.ome.tif_stardist.tiff'   ]));

    for i1 = 1:length(loc_filterd(:,1))
    loc_filterd(i1,9) = cell_img(round(loc_filterd(i1,2)),round(loc_filterd(i1,1)));
    end
    else 
        loc_filterd = loc;
    end
    loc_all = [loc_all ;loc_filterd];
    
end

save( 'D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\noncoding_spot_assignment.txt','loc_all','-ascii','-tabs');
    

% a = double(imread('Z:\CL\20220806coding_jisui\z_projection_0806\z_round1\max_jisui0806_round11_Cy7_s1.ome.tif'));
% figure
% imshow(cell_img,[0,1]);
% hold on 
% scatter(loc_filterd(:,1),loc_filterd(:,2),5,'green')    