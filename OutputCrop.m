%multithreading magic
    % matlabpool open 4;

%script to run custom crop function on all output data.
    clc;clear all;
    
%specify path
    path='E:\Data2\Ely_May28th\Output charts\E:\Data2\Ely_May28th\Output charts\Profile Axial Velocity\1\';
    cd(path);
    d=dir('*.jpg');
    
parfor z=1:length(d)
    BatchCropDirectory(strcat(path,d(z).name),20,'jpg');
    cd(path);
end

disp('Done!');
    % matlabpool close;
