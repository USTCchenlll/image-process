%因为第二轮dapi更亮

clc
clear

stage_num = 949;

for stage_i = 1:stage_num
    disp([num2str(stage_i) '/949'])
    cell_img = double(imread(['Z:\CL\20221003_jisui_no18\z_projection\cell_img\2_raw_z_cell\max_20221003_jisui18_round11_DAPI_s' num2str(stage_i) '.ome.tif_stardist.tiff'  ]));
    cell_img = imtranslate(cell_img,[-16 1]);
    imwrite(uint16(cell_img),['Z:\CL\20221003_jisui_no18\z_projection\cell_img\3_drift_z_cell\max_20221003_jisui18_no1_round11_DAPI_s' num2str(stage_i) '.ome.tif_stardist.tiff'])
end
     