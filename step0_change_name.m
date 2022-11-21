clc
clear
cd Z:\CL\20221003_jisui_no18\raw_data\round1
raw_name = dir('Z:\CL\20221003_jisui_no18\raw_data\round1\*.tif');

for i1 = 1:length(raw_name)
disp(i1)
a  = raw_name(i1).name;
b = length(a);
copyfile(raw_name(i1).name, [a(1:16) '_round11' a(18:b)  ]);  % 把2.txt复制成22.txt(2.txt依然存在)  



end