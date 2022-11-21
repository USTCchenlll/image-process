




clc
clear

drift_location  = 'D:\CL\jisui\jisui_image_analysis\20220912_img_flow\intermediate_data\1_spot_data\offset_drift\';

% % Input parameters
ztack_image_file_location = 'Z:\CL\20220910_jisui_no7\z_projection\';
position_file = 'D:\CL\jisui\jisui_image_analysis\20220912_img_flow\intermediate_data\2_stage_data\position_0910_final.txt'; %only list of location
image_size = [2048,2048]; %[x,y]
overlap_percent = 0.1;
% get each slice fov
file_output_location = 'D:\CL\jisui\jisui_image_analysis\20220912_img_flow\intermediate_data\2_stage_data\';
file_name_qianzhui = '20220910jisuino7';
slice_cut_parameter = importdata([file_output_location 'all_slice_position0910.xlsx']);
slice_cut_parameter = slice_cut_parameter.data.Sheet1;
slice_num = length(slice_cut_parameter(:,1));
 
% get the overlap size
overlap_size = round(image_size*overlap_percent);

% Determine the coordinates(确定所有图像的坐标)
raw_position = importdata(position_file);
for i1 = 1:length(raw_position)    
   temp_rawlocation(1,:) = strsplit(raw_position{i1},',');   
   x_position = strfind(temp_rawlocation{1,1},'x');
   y_position = strfind(temp_rawlocation{1,1},'y');
   d_position = strfind(temp_rawlocation{1,1},'"');
   position_list(i1,1:length(temp_rawlocation)) = temp_rawlocation(1,:);
   position_list{i1,14} = str2double( temp_rawlocation{1,1}(x_position+1:y_position-1) );
   position_list{i1,15} = str2double( temp_rawlocation{1,1}(y_position+1:d_position(2)-1) );
end

clear raw_position x_position y_position d_position

% mapping the image location 
xmax = max(cell2mat(position_list(:,14)));
ymax = max(cell2mat(position_list(:,15)));
weizhi_list = [cell2mat(position_list(:,14)) cell2mat(position_list(:,15))]; % [x,y]
for i1 = 1:xmax
    for i2 = 1:ymax
        list_sum(i2,i1) = ismember([i1,i2],weizhi_list,'rows') ; 
        if  ismember([i1,i2],weizhi_list,'rows')   
            map_location(i2,i1) = intersect(find(weizhi_list(:,1)==i1 ),find(weizhi_list(:,2)==i2 ));
        end  
    end
end
map_location(9,2) = 0;

for i1 = 1:slice_num
    slice_part = map_location(slice_cut_parameter(i1,3):slice_cut_parameter(i1,3)+slice_cut_parameter(i1,5)-1,slice_cut_parameter(i1,2):slice_cut_parameter(i1,2)+slice_cut_parameter(i1,4)-1);
    [y1,x1] = size(slice_part);
    all_slice_file_list(1:y1,1:x1,i1) = slice_part;   
    
end


% save([file_output_location '\all_slice_file_list.mat'], 'all_slice_file_list')
%%  stitch each slice

    coordinate_stage = cell(size(all_slice_file_list));
for slice_i = 1 : slice_num 
    disp([ 'slice_' num2str(slice_i)])
    %start from y=4 up and down
    temp_slice_stage = all_slice_file_list(:,:,slice_i);
    big_slice_image = zeros((image_size-overlap_size).*size(temp_slice_stage)+400);

    
    coordinate_stage{4,1,slice_i} = [  (image_size(2)-overlap_size(2)).*4+201 ,201 ]; 
    img_start = double(imread([ztack_image_file_location 'z_round1\max_' file_name_qianzhui '_round11_DAPI_s'  num2str(temp_slice_stage(4,1))   '.ome.tif'  ]));
    img_start_gfp = double(imread([ztack_image_file_location 'z_round5\max_' file_name_qianzhui '_round51_DAPI_s'  num2str(temp_slice_stage(4,1))   '.ome.tif'  ]));
    drift_c = load([ drift_location 'noncoding_offset_drift_s' num2str(temp_slice_stage(4,1))  '.txt'  ]);
    img_start_gfp = imtranslate(img_start_gfp,[-drift_c(5,2),-drift_c(5,1)]);
    
	big_slice_image(coordinate_stage{4,1,slice_i}(1)+1:coordinate_stage{4,1,slice_i}(1)+image_size(2),coordinate_stage{4,1,slice_i}(2)+1:coordinate_stage{4,1,slice_i}(2)+image_size(2)) =	img_start_gfp;
	for i1 = 2:length(temp_slice_stage(1,:))
        if temp_slice_stage(4,i1)>0
        img_next = double(imread([ztack_image_file_location 'z_round1\max_' file_name_qianzhui '_round11_DAPI_s'  num2str(temp_slice_stage(4,i1))   '.ome.tif'  ]));  
        img_next_gfp = double(imread([ztack_image_file_location 'z_round5\max_' file_name_qianzhui '_round51_DAPI_s'  num2str(temp_slice_stage(4,i1))   '.ome.tif'  ]));
        drift_c = load([ drift_location 'noncoding_offset_drift_s' num2str(temp_slice_stage(4,i1))  '.txt'  ]);
        img_next_gfp = imtranslate(img_next_gfp,[-drift_c(5,2),-drift_c(5,1)]);
        
        C = normxcorr2(img_start(:,image_size(1)-overlap_size(1)+1:image_size(1)), img_next(:,1:overlap_size(1)));   % 相关性计算 
        cy = round(size(C,1)/2); cx = round(size(C,2)/2);
        maxd = 100;
        CC = C((cy-maxd):(cy+maxd), (cx-maxd):(cx+maxd)); 
        [max_cc, imax] = max(abs(CC(:)));       
        [ypeak, xpeak] = ind2sub(size(CC),imax(1));
        ypeak = ypeak + cy-maxd; xpeak = xpeak + cx - maxd;
        temp = [(ypeak-image_size(1))-1 (xpeak-overlap_size(2))-1]; %漂移校准坐标差值
               
        
        coordinate_stage{4,i1,slice_i}(1) = coordinate_stage{4,i1-1,slice_i}(1)-temp(2); 
        coordinate_stage{4,i1,slice_i}(2) = coordinate_stage{4,i1-1,slice_i}(2)+image_size(2)-overlap_size(2)-temp(1);
        big_slice_image(coordinate_stage{4,i1,slice_i}(1)+1:coordinate_stage{4,i1,slice_i}(1)+image_size(2),coordinate_stage{4,i1,slice_i}(2)+1:coordinate_stage{4,i1,slice_i}(2)+image_size(2)) = img_next_gfp;
        img_start = img_next;
        end
    
	end
    
    

	for i1 = 1:length(temp_slice_stage(1,:))
        for i2 = 1:3
        i2 = 3-i2+1;
        if temp_slice_stage(i2,i1)>0
            disp([num2str(i2) '-' num2str(i1) ])
        img_start = double(imread([ztack_image_file_location 'z_round1\max_' file_name_qianzhui '_round11_DAPI_s'  num2str(temp_slice_stage(i2+1,i1))   '.ome.tif'  ]));
        img_next = double(imread([ztack_image_file_location 'z_round1\max_' file_name_qianzhui '_round11_DAPI_s'  num2str(temp_slice_stage(i2,i1))   '.ome.tif'  ]));
        img_next_gfp = double(imread([ztack_image_file_location 'z_round5\max_' file_name_qianzhui '_round51_DAPI_s'  num2str(temp_slice_stage(i2,i1))   '.ome.tif'  ]));
        drift_c = load([ drift_location 'noncoding_offset_drift_s' num2str(temp_slice_stage(i2,i1))  '.txt'  ]);
        img_next_gfp = imtranslate(img_next_gfp,[-drift_c(5,2),-drift_c(5,1)]);
        
        C = normxcorr2(img_start(1:overlap_size(1),:), img_next(image_size(1)-overlap_size(1)+1:image_size(1),:));   % 相关性计算 
        cy = round(size(C,1)/2); cx = round(size(C,2)/2);
        maxd = 100;
        CC = C((cy-maxd):(cy+maxd), (cx-maxd):(cx+maxd)); 
        [max_cc, imax] = max(abs(CC(:)));       
        [ypeak, xpeak] = ind2sub(size(CC),imax(1));
        ypeak = ypeak + cy-maxd; xpeak = xpeak + cx - maxd;
        temp = [(ypeak-overlap_size(1))-1 (xpeak-image_size(2))-1]; %漂移校准坐标差值
               
        coordinate_stage{i2,i1,slice_i}(1) = coordinate_stage{i2+1,i1,slice_i}(1)-(image_size(2)-overlap_size(2))-temp(2); 
        coordinate_stage{i2,i1,slice_i}(2) = coordinate_stage{i2+1,i1,slice_i}(2)-temp(1);
        
        big_slice_image(coordinate_stage{i2,i1,slice_i}(1)+1:coordinate_stage{i2,i1,slice_i}(1)+image_size(2),coordinate_stage{i2,i1,slice_i}(2)+1:coordinate_stage{i2,i1,slice_i}(2)+image_size(2)) = img_next_gfp;
        end
        end
        
	end
 
	for i1 = 1:length(temp_slice_stage(1,:))
        for i2 = 5:length(temp_slice_stage(:,1))
            
        if temp_slice_stage(i2,i1)>0
        img_start = double(imread([ztack_image_file_location 'z_round1\max_' file_name_qianzhui '_round11_DAPI_s'  num2str(temp_slice_stage(i2-1,i1))   '.ome.tif'  ]));
        img_next = double(imread([ztack_image_file_location 'z_round1\max_' file_name_qianzhui '_round11_DAPI_s'  num2str(temp_slice_stage(i2,i1))   '.ome.tif'  ]));
        img_next_gfp = double(imread([ztack_image_file_location 'z_round5\max_' file_name_qianzhui '_round51_DAPI_s'  num2str(temp_slice_stage(i2,i1))   '.ome.tif'  ]));
        drift_c = load([ drift_location 'noncoding_offset_drift_s' num2str(temp_slice_stage(i2,i1))  '.txt'  ]);
        img_next_gfp = imtranslate(img_next_gfp,[-drift_c(5,2),-drift_c(5,1)]);
        C = normxcorr2(img_start(image_size(1)-overlap_size(1)+1:image_size(1),:), img_next(1:overlap_size(1),:));   % 相关性计算 
        cy = round(size(C,1)/2); cx = round(size(C,2)/2);
        maxd = 100;
        CC = C((cy-maxd):(cy+maxd), (cx-maxd):(cx+maxd)); 
        [max_cc, imax] = max(abs(CC(:)));       
        [ypeak, xpeak] = ind2sub(size(CC),imax(1));
        ypeak = ypeak + cy-maxd; xpeak = xpeak + cx - maxd;
        temp = [(ypeak-overlap_size(1))-1 (xpeak-image_size(2))-1]; %漂移校准坐标差值
               
        coordinate_stage{i2,i1,slice_i}(1) = coordinate_stage{i2-1,i1,slice_i}(1)+(image_size(2)-overlap_size(2))-temp(2); 
        coordinate_stage{i2,i1,slice_i}(2) = coordinate_stage{i2-1,i1,slice_i}(2)-temp(1);
        big_slice_image(coordinate_stage{i2,i1,slice_i}(1)+1:coordinate_stage{i2,i1,slice_i}(1)+image_size(2),coordinate_stage{i2,i1,slice_i}(2)+1:coordinate_stage{i2,i1,slice_i}(2)+image_size(2)) = img_next_gfp;
        end
        end
        
	end    
    
    big_slice_image = imtranslate(big_slice_image,[-16,1]);
    imwrite(uint16(big_slice_image),[ file_output_location 'slice_'  num2str(slice_i)  'stitch_dapi_round5.tif' ])
    
    
end









