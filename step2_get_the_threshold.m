%test for the threshold

clc
clear
nstack = 14;
stage_num = 949;
file_base_location = 'Z:\CL\20221003_jisui_no18\raw_data\';
file_name_qianzhui = '20221003_jisui18';
stage_i = 864;





Round = 5;
coding_round = 10;
maxd = 100; %相关性漂移区域
dim = [2048,2048];
bpass_lnoise = 2.5;
bpass_lobject = 4;
pkfnd_size = 4;    %play with this - try 3.
cntrd_region = 9;


%pkfnd_thresh = [200,200,160,160,120,120,80,80,50,50].';   %单分子定义的阈值，这个不同通道需要调整
% generate file name list of this stage
for i = 1:Round
img_file{i*2-1,1} = [file_base_location 'round' num2str(i) '\' file_name_qianzhui '_round' num2str(i) '1_Cy5_s' num2str(stage_i) '.ome.tif'  ];
img_file{i*2,1}   = [file_base_location 'round' num2str(i) '\' file_name_qianzhui '_round' num2str(i) '1_Cy7_s' num2str(stage_i) '.ome.tif'  ];
end

% input the dapi image  
dapi = zeros(dim(1),dim(2),Round);
for i=1:Round
    file = fullfile([file_base_location 'round' num2str(i) '\' file_name_qianzhui '_round' num2str(i) '1_DAPI_s' num2str(stage_i) '.ome.tif'  ]);
    fileinfo = imfinfo(file);
    temp = zeros(dim(1),dim(2),nstack);
    for frame=1:nstack
        temp(:,:,frame)=imread(file,frame,'Info', fileinfo);
    end
    temp = max(temp,[],3); %maximal projection
    dapi(:,:,i) = temp(:,:,1);
end

% drift correct and save it
offset_raw = [0,0,1];
offset_drift(1,:) = offset_raw;
offset_drift(2,:) = offset_raw;
for i = 2:Round
    C = normxcorr2(dapi(:,:,1), dapi(:,:,i));   % 相关性计算
    % find peak in the center of the correlation map
    cy = round(size(C,1)/2); cx = round(size(C,2)/2);
    CC = C((cy-maxd):(cy+maxd), (cx-maxd):(cx+maxd)); 
    [max_cc, imax] = max(abs(CC(:)));       
    [ypeak, xpeak] = ind2sub(size(CC),imax(1));
    ypeak = ypeak + cy-maxd; xpeak = xpeak + cx - maxd;
    temp = [(ypeak-dim(1))-1 (xpeak-dim(2))-1 i]; %漂移校准坐标差值
    offset_raw = [offset_raw;temp];
    offset_drift(i*2-1,:) = temp ;
    offset_drift(i*2,:) = temp;
end

%%
pkfnd_thresh = [200,220,200,220,200,220,180,200,140,160].';   %单分子定义的阈值，这个不同通道需要调整
% find spot by bpass and pkfnd
res = [];
r1 = 10;

for r= r1:r1
    clear loc
    fileinfo = imfinfo(img_file{r,1});
    im0 = zeros(dim(1),dim(2),nstack);
    for frame = 1:nstack
        im0(:,:,frame)=imread(img_file{r,1},frame,'Info', fileinfo);
    end
    im = zeros(dim(1),dim(2));
    loc = [];
    for i =8:8
        im(:,:) = im0(:,:,i);
        B0 = bpass(im,bpass_lnoise,bpass_lobject);
        Peak = pkfnd(B0,pkfnd_thresh(r,1),pkfnd_size);
        if ~isempty(Peak)
            margin_peak = [find(Peak(:,1)<cntrd_region/2+2) ; find(Peak(:,1)>dim(1)- cntrd_region/2-2) ; find(Peak(:,2)>dim(1)- cntrd_region/2-2  ) ; find(Peak(:,2)<cntrd_region/2+2) ];
            Peak(margin_peak,:) = [];
            P = get_averagedcenter(B0, Peak, cntrd_region, 2);
            P(:,1) = P(:,1)-offset_drift(r,2);
            P(:,2) = P(:,2)-offset_drift(r,1);
            temp = [P(:,1:2) P(:,1).*0+i P(:,1).*0 P(:,1).*0+r];
            loc = [loc;temp];
        end
    end  
    
    
  figure   
  imshow(imtranslate(im,[-offset_drift(r,2),-offset_drift(r,1)]),[100,max(max(im))*0.3]) 
  hold on 
  scatter(loc(:,1),loc(:,2),5,'green')      

end



 



    % 目前的方法有潜在缺陷：结果很可能与遍历的顺序相关
    % 理想情况下，可对同一批数据进行两种或多种不同顺序的遍历，丢弃任何一种情况下被判断为无法解码的分子，以提高正确率
 


%%
% --------------------------------------------------------------------------------------------------
%%

%test for all image spot noncoding round12345 561channel

clc
clear
nstack = 14;
stage_num = 949;
file_base_location = 'Z:\CL\20221003_jisui_no18\raw_data\';
file_name_qianzhui = '20221003_jisui18';
stage_i = 360;



Round = 5;
maxd = 100; %相关性漂移区域
dim = [2048,2048];

% coding_round = 10;


bpass_lnoise = 2.5;
bpass_lobject = 4;
pkfnd_size = 4;    %play with this - try 3.
cntrd_region = 9;
%pkfnd_thresh = [100,90,80,75,70].';   %单分子定义的阈值，这个不同通道需要调整




% generate file name list of this stage jisui_0830_no1_round11_RFP-filter_s1.ome
for i = 1:5
img_file{i,1}   = [file_base_location 'round' num2str(i) '\' file_name_qianzhui '_round' num2str(i) '1_RFP-filter_s' num2str(stage_i) '.ome.tif'  ];
end
% input the dapi image  
dapi = zeros(dim(1),dim(2),5);
for i=1:5
    file = fullfile([file_base_location 'round' num2str(i) '\' file_name_qianzhui '_round' num2str(i) '1_DAPI_s' num2str(stage_i) '.ome.tif'  ]);
    fileinfo = imfinfo(file);
    temp = zeros(dim(1),dim(2),nstack);
    for frame=1:nstack
        temp(:,:,frame)=imread(file,frame,'Info', fileinfo);
    end
    temp = max(temp,[],3); %maximal projection  
    dapi(:,:,i) = temp(:,:,1);
end
% drift correct and save it
offset_raw = [0,0,1];
offset_drift(1,:) = offset_raw;
for i = 2:5
    C = normxcorr2(dapi(:,:,1), dapi(:,:,i));   % 相关性计算
    % find peak in the center of the correlation map
    cy = round(size(C,1)/2); cx = round(size(C,2)/2);
    CC = C((cy-maxd):(cy+maxd), (cx-maxd):(cx+maxd)); 
    [max_cc, imax] = max(abs(CC(:)));       
    [ypeak, xpeak] = ind2sub(size(CC),imax(1));
    ypeak = ypeak + cy-maxd; xpeak = xpeak + cx - maxd;
    temp = [(ypeak-dim(1))-1 (xpeak-dim(2))-1 i]; %漂移校准坐标差值
    offset_raw = [offset_raw;temp];
    offset_drift(i,:) = temp ;
end

%%



pkfnd_thresh = [70,75,70,70,70].';   %单分子定义的阈值，这个不同通道需要调整
% find spot by bpass and pkfnd
res = [];
rr= 3 ;
for r=rr:rr
    fileinfo = imfinfo(img_file{r,1});
    im0 = zeros(dim(1),dim(2),nstack);
    for frame = 1:nstack
        im0(:,:,frame)=imread(img_file{r,1},frame,'Info', fileinfo);
    end
    im = zeros(dim(1),dim(2));
    loc = [];
    for i =7:7
        im(:,:) = im0(:,:,i);
        B0 = bpass(im,bpass_lnoise,bpass_lobject);
        Peak = pkfnd(B0,pkfnd_thresh(r,1),pkfnd_size);
        if ~isempty(Peak)
            margin_peak = [find(Peak(:,1)<cntrd_region/2+2) ; find(Peak(:,1)>dim(1)- cntrd_region/2-2) ; find(Peak(:,2)>dim(1)- cntrd_region/2-2  ) ; find(Peak(:,2)<cntrd_region/2+2) ];
            Peak(margin_peak,:) = [];
            P = get_averagedcenter(B0, Peak, cntrd_region, 2);
            P(:,1) = P(:,1)-offset_drift(r,2);
            P(:,2) = P(:,2)-offset_drift(r,1);
            temp = [P(:,1:2) P(:,1).*0+i P(:,1).*0 P(:,1).*0+r];
            loc = [loc;temp];
        end
    end
% loc有5列：x, y, #z-stack, 离焦标记，轮次
end
 figure   
  imshow(imtranslate(im,[-offset_drift(r,2),-offset_drift(r,1)]),[100,max(max(im))*0.3]) 
  hold on 
  scatter(loc(:,1),loc(:,2),5,'green') 
% for r=1:1
%     fileinfo = imfinfo(img_file{r,1});
%     im0 = zeros(dim(1),dim(2),nstack);
%     for frame = 1:nstack
%         im0(:,:,frame)=imread(img_file{r,1},frame,'Info', fileinfo);
%     end
%     im = zeros(dim(1),dim(2));
%     loc = [];
%     for i =7:7
%         im(:,:) = im0(:,:,i);
%         B0 = bpass(im,bpass_lnoise,bpass_lobject);
%         Peak = pkfnd(B0,pkfnd_thresh(r,1),pkfnd_size);
%         if ~isempty(Peak)
%             margin_peak = [find(Peak(:,1)<cntrd_region/2+2) ; find(Peak(:,1)>dim(1)- cntrd_region/2-2) ; find(Peak(:,2)>dim(1)- cntrd_region/2-2  ) ; find(Peak(:,2)<cntrd_region/2+2) ];
%             Peak(margin_peak,:) = [];
%             P = get_averagedcenter(B0, Peak, cntrd_region, 2);
%             P(:,1) = P(:,1)-offset_drift(r,2);
%             P(:,2) = P(:,2)-offset_drift(r,1);
%             temp = [P(:,1:2) P(:,1).*0+i P(:,1).*0 P(:,1).*0+r];
%             loc = [loc;temp];
%         end
%     end
% % loc有5列：x, y, #z-stack, 离焦标记，轮次
% end
% 
%   scatter(loc(:,1),loc(:,2),5,'red') 
 

% save(fullfile(['s' num2str(stage_i) '_filtered_decoding.txt']), 'loc', '-ascii', '-tabs')
%  scatter(loc_final(:,1),loc_final(:,2),5,'green')          




















