% get all cell location and ronghe them and get off the overlap cell in near fov

clc
clear 

stage_num = 949;

% coding_spot = load('D:\CL\jisui\jisui_image_analysis\20220830_img_flow\intermediate_data\spots_data\coding_spot_assignment.txt');
% noncoding_spot = load('D:\CL\jisui\jisui_image_analysis\20220830_img_flow\intermediate_data\spots_data\noncoding_spot_assignment.txt');

load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\2_stage_data\coordinate_stage1003.mat');
load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\2_stage_data\all_slice_file_list.mat');

slice_num = 16;
image_size = [2048,2048];
overlap_size = [205,205]; 

% cell_list :  1(code for this cell in this stage) 2(stage_i) 3(x) 4(y)
% 5(unique_id) 6(if overlap 1=overlap) 7(slice_i) 8(area) 
cell_list = zeros(0,8);
num_cell = 0;

%%
for slice_i = 1:slice_num
    tic
    disp(['slice_' num2str(slice_i)])
    clear temp_boundary temp_cell_list
    
   
    temp_coordinate_stage = coordinate_stage(:,:,slice_i);
    temp_slice = all_slice_file_list(:,:,slice_i);
    stage_list = sort(unique(temp_slice));
    stage_list(1,:) = [];
    
	% stage_boundary
    for stage_i = 1:length(stage_list)
        [y1,x1] = find(temp_slice == stage_list(stage_i));
        if  y1<length(temp_slice(:,1)) && temp_slice(y1+1,x1)>0
            temp_boundary{y1,x1}(1) = temp_coordinate_stage{y1+1,x1}(1)-1;
        else
            temp_boundary{y1,x1}(1) = temp_coordinate_stage{y1,x1}(1)+image_size(1)-1;
        end
        if x1<length(temp_slice(1,:)) && temp_slice(y1,x1+1)>0
            temp_boundary{y1,x1}(2) = temp_coordinate_stage{y1,x1+1}(2)-1;
        else
            temp_boundary{y1,x1}(2) = temp_coordinate_stage{y1,x1}(2)+image_size(2)-1;
        end
    end
    
    for stage_i = 1:length(stage_list)   
        [y1,x1] = find(temp_slice == stage_list(stage_i));
        dapi_img = double(imread(['Z:\CL\20221003_jisui_no18\z_projection\cell_img\4_ex_z_cell\max_20221003_jisui18_no1_round11_DAPI_s' num2str(stage_list(stage_i,1))  '.ome.tif_stardist.tiff'   ]));
        if max(max(dapi_img))>0
        temp_cell_list = sort(unique(dapi_img));
        temp_cell_list(1,:) = [];
        temp_cell_list(:,2) = stage_list(stage_i,1);
        temp_cell_list(:,5) = temp_cell_list(:,1)+num_cell;
        temp_cell_list(:,6) = 0; 
        temp_cell_list(:,7) = slice_i;
        num_cell = num_cell+max(temp_cell_list(:,1));
        for cell_i = 1:length(temp_cell_list(:,1))
        [y,x] = find(dapi_img == temp_cell_list(cell_i,1));
        temp_cell_list(cell_i,8) = length(y);
        temp_cell_list(cell_i,3) = round(mean(x))+temp_coordinate_stage{temp_slice==stage_list(stage_i,1)}(2); %x
        temp_cell_list(cell_i,4) = round(mean(y))+temp_coordinate_stage{temp_slice==stage_list(stage_i,1)}(1); %y
        if min(x)+temp_coordinate_stage{y1,x1}(2)>=temp_boundary{y1,x1}(2) || min(y)+temp_coordinate_stage{y1,x1}(1)>= temp_boundary{y1,x1}(1)
            temp_cell_list(cell_i,6) = 1;       
        end
        end
        
        cell_list = [cell_list ; temp_cell_list];
        end
    end
    
   
    toc
end

save('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\3_cell_data\cell_list1008_dilate.mat','cell_list' );




























