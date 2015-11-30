clc; clear all;
%for loading 3d into averages and dynamic data files
tic();

for i=[35];
    close all;
    %Load3d(i);
    %tts(strcat('Run',' ',num2str(i),' has been loaded and averaged'));
    [Rcore(i),Haz(i)]=Vanalysis(i,'y','n','n');
    %tts(strcat('Run',' ',num2str(i),' has been analyzed'));
end
% output(:,2)=Haz';output(:,1)=Rcore';
% xlswrite('output.xls',output);

elapsed_time=toc()