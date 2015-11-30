%image pairs

drive='E';
d=':\Data2\Ely_May28th\Image\Dual\';
Run=70;
R=num2str(Run);


for i=1:1
    ii=num2str(i);
    if j <10
        path=strcat(drive,d,R,'\Ely_May28th',R,'00',ii; 
    elseif j<100
        path=strcat(drive,d,R,'\Ely_May28th',R,'0',ii);   
    elseif j >=100
        path=strcat(drive,d,R,'\Ely_May28th',R,ii); 
    end
    
    
    
    
    
    output=ImageTo2dVectors(path,calibration,16,3);
    
    
    
    
    
end