%this program is for z_projection


clc
clear

input_file_location = 'Z:\CL\20221003_jisui_no18\raw_data\round4\';
output_file_location = 'Z:\CL\20221003_jisui_no18\z_projection\round4\';

z = 14;

raw_data = dir([input_file_location '*.tif']);
%raw_data = raw_data(1899:2847);
%p=parpool(4);
raw_fig = zeros(2048,2048,z);
%parfor i1=1:length(raw_data)
for i1=1:length(raw_data)
    disp(i1)  
    for j1 = 1:z
       raw_fig(:,:,j1) = imread([input_file_location raw_data(i1).name],j1);         
    end
    max_f = uint16(max(raw_fig,[],3));  
    imwrite(max_f,[output_file_location 'max_' raw_data(i1).name])
end
%delete(p);
