% Script functions as the basic interface to process 3d vector files into
% fully analyzed data sets and figures. This is effectively a control
% panel.

clc; clear all;
tic();

Runs=[1:38 41:70];   % the runs to process.
drive ='E';         % the drive leter where the data is being stored.

VaxProf=zeros(10000,2,70);


for i=Runs;
    close all; 
 
    %Function contains aditional user defined parameters, if errors occur,
    %examine the filepaths compiled towards the top of the code.
    [Rcore(i),Haz(i),Hax(i),Lax(i),Vaz{i},VazProfile{i},...
        Vax{i},VaxProfile{i},Vra{i}]=Vanalysis(i,'y','n','n',drive);
                                 
 fprintf('Run number %2.0f Finished! \n\n',i);   
end

output(:,2)=Haz';
output(:,1)=Rcore'; 
output(:,3)=Hax'; 
output(:,4)=Lax';
output(:,5)=Wand';
xlswrite('output.xls',output);

elapsed_time=toc()/60;
fprintf('Completed in %3.2f minutes!',elapsed_time);