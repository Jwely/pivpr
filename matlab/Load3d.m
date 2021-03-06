function [timer]=Load3d(Run,drive)

%Code to process 3d vector files, take averages, remove bad vectors.etc
clc;tic(); fclose('all');

%inputs
    %Run=70;                  %Set the run to process valid numbers are 
                             %integers 1:17,1.5 and 2.5

%make directory (d) of only .v3d files
    path=strcat(drive,':\Data2\Ely_May28th\Vector\',num2str(Run));
    cd(path);
    d=dir('*.v3d');
    
%load Velocity index file, filter velocities more than 120% free stream.
    fid=fopen(strcat(drive,':\Data2\Vindex.txt'));
    Vin=importdata(strcat(drive,':\Data2\Vindex.txt'));
    Vfree=Vin(Run,2);

%If the AvgData file already exists, attempt to delete and replace it.
    if exist('AvgData.xls','file')==2
        delete('AvgData.xls');
    else end

%Information update
fprintf('Run # %2.1f \n Free Stream Velocity = %2.2f m/s\n\n',Run,Vfree);

%load and import all the .v3d text files into one matrix.    
for k=1:length(d)
    A=importdata(d(k).name);
    B=A.data;

    %filter inifnity values(errors)
    for i=1:size(B,1)
        for j=4:6  
            if B(i,j)>=Vfree*1.60              
                B(i,j)=0;
            else end
        end
    end
    
    %format (x,y,z,u,v,w) in units (mm ,mm ,mm ,m/s ,m/s ,m/s )
        Data(:,:,k)=B(:,1:6);
        if rem(k/20,1)==0
        fprintf('processed %4.0f /%4.0f \n',k,length(d));
        else end
end

%saves dynamic data matrix as a reshape file. this includes all values for
%every data point
    DynData=reshape(Data,size(B,1),6*length(d));
    csvwrite('DynamicData.txt',DynData);
    disp('Finished writing dynamic data');

%status update
    disp('Beginning averaging operation'); 

%finding the averages of the vector data

%fetch coordinates from first k of Data
    AvgData(:,1:3)=Data(:,1:3,1);
      
%AvgData(a,:) (columnas) format ==> (x,y,z,u,v,w,Uu,Uv,Uw,Npoints)
%where Us are standard deviations (U for Uncertainty)
   count(1:size(B,1),4:6)=zeros();count(1:size(B,1),1:3)=200;
   SDlist=zeros(1);
    
%counter used to ignore 0 values, count the number of good values for
%taking the standard deviation and mean
    for i=1:size(B,1)
        for j=4:6                             %for u,v,w
            for k=1:length(d)                 %for all runs (usually 200)
                if Data(i,j,k) ~=0
                    count(i,j)=count(i,j)+1;
                    SDlist(count(i,j))=Data(i,j,k);
                else end
            end
            AvgData(i,j+3)=std(SDlist);
            SDlist=zeros(1);
        end  
    end
 %save the count to AvgData matrix  
 AvgData(:,10)=count(:,6);
 
 %taking the mean of all good data points for each  u,v,w (j=4:6)   
    for i=1:size(Data,1)
        for j=4:6
            if count(i,j)>=3    %at least 3 good data points are needed
                AvgData(i,j)=sum(Data(i,j,:))/count(i,j);
                if AvgData(i,j)>=1.5*Vfree
                    AvgData(i,j)=-0.01;
                end
            else 
                AvgData(i,j)=-0.01;     %to avoid  those ugly NaN values.
                AvgData(i,j+3)=-00.01;
            end
        end
    end
    
%add another 3 columnsfor uncertainty/value. percentage of measurement.
    disp('Refining uncertainty values');

for j=4:6
    for i=1:size(AvgData,1)
        if AvgData(i,j+5)~=0
            AvgData(i,j+7)=abs(AvgData(i,j+3)/AvgData(i,j));
        else
            AvgData(i,j+7)=0;
        end
        
        if abs(AvgData(i,j+7))>=100.0
            AvgData(i,11:13)=0;
        end
        
    end
    
end

%Use some simple logic for filtering out bad vectors.

    
%save AvgData information
    xlswrite('AvgData.xls',AvgData);
    
%plot AvgData    
%     figure(2); axis equal;grid on;
%     quiver(AvgData(:,1),AvgData(:,2),AvgData(:,4),AvgData(:,5));
%     xlabel('mm');ylabel('mm');title('Flowfield')

%status update
    disp(path);disp('Finished'); fclose('all')
%end timer
    timer=toc();
end


