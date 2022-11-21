% 对识别出的点信号进行细胞的分配 for coding

clc
clear
%%

stage_num = 949;
loc_all = zeros(0,8);

for stage_i = 1:stage_num
    disp(num2str(stage_i))
	
    loc = load(['D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\coding\pair_coding_s' num2str(stage_i) '_decoding.txt'],'-ascii') ;
    choose_spot = loc(:,6)-loc(:,5);
    loc(choose_spot<0,:) = [];
    if ~isempty(loc)
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

save( 'D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\pair_coding_spot_assignment.txt','loc_all','-ascii','-tabs');
    

% a = double(imread('Z:\CL\20220806coding_jisui\z_projection_0806\z_round1\max_jisui0806_round11_Cy7_s1.ome.tif'));
% figure
% imshow(cell_img,[0,1]);
% hold on 
% imwrite(uint16(cell_img2),'cell2.tif')
% imwrite(uint16(dapi_img),'dapi.tif')
% scatter(loc_filterd(:,1),loc_filterd(:,2),5,'green')    



