clc
clear 

load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\pair_coding_spot_assignment.txt');
load('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\noncoding_spot_assignment.txt');

all_spot_unfilter1008 = [pair_coding_spot_assignment;noncoding_spot_assignment];
save('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\all_spot_unfilter1008.mat','all_spot_unfilter1008' );


pair_coding_spot_assignment(pair_coding_spot_assignment(:,9)==0,:)=[];
noncoding_spot_assignment(noncoding_spot_assignment(:,9)==0,:)=[];

all_spot_1008 =  [pair_coding_spot_assignment;noncoding_spot_assignment];
save('D:\CL\jisui\jisui_image_analysis\20221003_img_flow\intermediate_data\1_spot_data\all_spot_1008.mat','all_spot_1008' );

all_spot_unfilter1008(all_spot_unfilter1008(:,6)>all_spot_unfilter1008(:,7) ,:)=[];