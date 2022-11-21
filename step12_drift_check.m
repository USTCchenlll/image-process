clc
clear


stage_list = [193,310,423,561,606,909,912,918,921,932,936,944];
dim = [2048,2048];
file_base_location = 'Z:\CL\20220806_coding_jisui\z_projection_0806\gfp_drift\';
maxd = 100;
for stage_i = 1:12
% input the dapi image  
dapi = zeros(dim(1),dim(2),3);
gfp = zeros(dim(1),dim(2),3);

for i= 1:3
    
    temp_dapi=double(imread([ file_base_location 'dapi\max_jisui0806_round' num2str(i-1) '1_DAPI_s' num2str(stage_list(stage_i)) '.ome.tif'  ]));    
    temp_gfp=double(imread([ file_base_location 'gfp\max_jisui0806_round' num2str(i-1) '1_GFP-only-3_s' num2str(stage_list(stage_i)) '.ome.tif'  ]));    
    
     temp_gfp = imtranslate(temp_gfp,[15,-3]);
    dapi(:,:,i) = temp_dapi;
    gfp(:,:,i) = temp_gfp;
    
end

% drift correct and save it
offset_raw = [0,0,1];
offset_drift_dapi(1,:) = offset_raw;
for i = 2:3
    C = normxcorr2(dapi(:,:,1), dapi(:,:,i));   % 相关性计算
    % find peak in the center of the correlation map
    cy = round(size(C,1)/2); cx = round(size(C,2)/2);
    CC = C((cy-maxd):(cy+maxd), (cx-maxd):(cx+maxd)); 
    [max_cc, imax] = max(abs(CC(:)));       
    [ypeak, xpeak] = ind2sub(size(CC),imax(1));
    ypeak = ypeak + cy-maxd; xpeak = xpeak + cx - maxd;
    temp = [(ypeak-dim(1))-1 (xpeak-dim(2))-1 i]; %漂移校准坐标差值
    offset_raw = [offset_raw;temp];
    offset_drift_dapi(i,:) = temp ;
end
for i1 = 1:3
    drift_img(:,:) = dapi(:,:,i1);
  %  imwrite(uint16(drift_img),['Z:\CL\20220806_coding_jisui\z_projection_0806\gfp_drift\nodrift_s' num2str(stage_i) '_dapi.tif'],'WriteMode','append');    
    drift_img(:,:) = imtranslate(drift_img(:,:),[-offset_drift_dapi(i1,2),-offset_drift_dapi(i1,1)]);
    imwrite(uint16(drift_img),['Z:\CL\20220806_coding_jisui\z_projection_0806\gfp_drift\drift_s' num2str(stage_i) '_dapi.tif'],'WriteMode','append');
end


offset_raw = [0,0,1];
offset_drift_gfp(1,:) = offset_raw;
for i = 2:3
    C = normxcorr2(dapi(:,:,1), dapi(:,:,i));   % 相关性计算
    % find peak in the center of the correlation map
    cy = round(size(C,1)/2); cx = round(size(C,2)/2);
    CC = C((cy-maxd):(cy+maxd), (cx-maxd):(cx+maxd)); 
    [max_cc, imax] = max(abs(CC(:)));       
    [ypeak, xpeak] = ind2sub(size(CC),imax(1));
    ypeak = ypeak + cy-maxd; xpeak = xpeak + cx - maxd;
    temp = [(ypeak-dim(1))-1 (xpeak-dim(2))-1 i]; %漂移校准坐标差值
    offset_raw = [offset_raw;temp];
    offset_drift_gfp(i,:) = temp ;
end
for i1 = 1:3
    drift_img(:,:) = gfp(:,:,i1);
  %  imwrite(uint16(drift_img),['Z:\CL\20220806_coding_jisui\z_projection_0806\gfp_drift\nodrift_s' num2str(stage_i) '_gfp.tif'],'WriteMode','append');    
    drift_img(:,:) = imtranslate(drift_img(:,:),[-offset_drift_gfp(i1,2),-offset_drift_gfp(i1,1)]);
    imwrite(uint16(drift_img),['Z:\CL\20220806_coding_jisui\z_projection_0806\gfp_drift\drift_s' num2str(stage_i) '_gfp.tif'],'WriteMode','append');
end




end