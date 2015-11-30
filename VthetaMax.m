%script to
clc;clear all; close all;

%inputs
    Run=8;                  % Set the run to process
    r=15;                    % Range from max point to look for min point 
                            % (in GRID steps, not mm)
                            
%open up all the processed data files   
    path=strcat('E:\Jeff\Ely_Apr5th\Vector\',num2str(Run));cd(path);
    fid=fopen('AvgData.xls','r');
    AvgData=xlsread('AvgData.xls');

%In a perfect world, all in plane velocities (Vp) should be perfectly 
%tangent to the center of the vortex, this will be assumed untill the core 
%is found. Vp=(:,:,[x,y,u,v,VpMag,w])
    Vp(:,1:2)= AvgData(:,1:2);
    Vp(:,3)  = AvgData(:,4);
    Vp(:,4)  = AvgData(:,5);
    Vp(:,5)  = (AvgData(:,4).^2+AvgData(:,5).^2).^0.5;
    Vp(:,6)  = AvgData(:,6);
    Vp

    %transform Vp into a matrix format that matlab likes more (Vm)
    k=1;j=1;
 for i=2:size(Vp,1)
     if Vp(i,2)==Vp(i-1,2)
         j=j+1;
     else
         k=k+1; j=1;
         Y(k-1)=Vp(i-1,2);
         Y(k)=Vp(i,2);
     end
         Vm(k,j)=Vp(i,5);
         %these 3 were added after I realized it was usefull to preserve
         %independent X Y Z components in the new format.
         Vmx(k,j)=Vp(i,3);
         Vmy(k,j)=Vp(i,4);
         Vmz(k,j)=Vp(i,6);
 end

 %get rid of NAN values. set them to 0
 for a=1:size(Vm,1)
     for b=1:size(Vm,2)
         if isnan(Vm(a,b))
             Vm(a,b)=0;
             Vmx(a,b)=0;
             Vmy(a,b)=0;
             Vmz(a,b)=0;
         else end
     end
 end

%extract X and Y coordinates as arrays from the data matrix.
    X=Vp(1:j,1); Y=Y';
 
%finding the maximum of Vm
    Max=max(max(Vm'));
    [Ymax Xmax]=find(Vm==Max);
    Yh=Y(Ymax); Xh=X(Xmax);

%looking for a minimum wich is near the maximum, attempting to hone in on
%the vortex core on the first try. range (r) is an input value at top.
    [Min]=min(min(Vm(Ymax-r:Ymax+r,Xmax-r:Xmax+r)'));
    [Ymin Xmin]=find(Vm==Min);
    Yl=Y(Ymin);Xl=X(Xmin);
    
%getting a more precise core location, since the core need not lay on top 
%of a preexisting grid point. Examines the adjacent points to Min point
%Core Zone (cz)             CZ(:,:,#)===> 1=V, 2=X coord 3=Y coord
    CZ(:,:,1)=Vm(Ymin-1:Ymin+1,Xmin-1:Xmin+1);
    for i=1:3
        for j=1:3
            CZ(i,j,2)=X(Xmin-2+j);
            CZ(i,j,3)=Y(Ymin-2+i);
        end
    end
  
%subtract the max value from velocity and take a weighted local average.    
    A=max(max(CZ(:,:,1)));
    CZ(:,:,1)=CZ(:,:,1)-A;
    CZ(:,:,4)=(CZ(:,:,1).*-CZ(:,:,2));
    CZ(:,:,5)=(CZ(:,:,1).*-CZ(:,:,3));
    
    Xcore=-sum(sum(CZ(:,:,4)))/sum(sum(CZ(:,:,1)));
    Ycore=-sum(sum(CZ(:,:,5)))/sum(sum(CZ(:,:,1)));
    
%Back to square one, now that the core has been found, we can calculate
%true and precise azimuthal velocities(Vaz) (ASUMING THE CORE IS CIRCULAR),
%this will need to be examined more closely if eliptical approximations
%prove necessary.
for xi=1:length(X)
    for yi=1:length(Y)
        %calculate azimuthal velocities (Vaz)
        A=[-(Y(yi)-Ycore);(X(xi)-Xcore);0]./...
                                (((Y(yi)-Ycore)^2+(X(xi)-Xcore)^2)^0.5);
        B=[Vmx(yi,xi);Vmy(yi,xi);0];
        Vaz(yi,xi)=-(dot(A,B));
        %calculate radial velocities (Vra)
        C=[(X(xi)-Xcore);(Y(yi)-Ycore);0]./...
                                (((Y(yi)-Ycore)^2+(X(xi)-Xcore)^2)^0.5);
        D=[Vmx(yi,xi);Vmy(yi,xi);0];
        Vra(yi,xi)=-(dot(C,D));
    end
end

%finding the maximum of Vaz, since it may not be the same as max Vm
    Maxaz=max(max(Vaz'));
    [Ymax Xmax]=find(Vaz==Maxaz);
    Yazh=Y(Ymax); Xazh=X(Xmax);
%finding the maximum of Vra
    Maxra=max(max(Vra'));
    [Ymax Xmax]=find(Vra==Maxra);
    Yrah=Y(Ymax); Xrah=X(Xmax);
%finding the maximum of Vaxial (Vmz)
    Maxax=max(max(Vmz'));
    [Ymax Xmax]=find(Vmz==Maxax);
    Yaxh=Y(Ymax); Xaxh=X(Xmax);
    
%calculate the core radius
    Rcore=((Xh-Xcore)^2+(Yh-Ycore)^2)^0.5;
    
%print out some of this data.   
    fprintf('Run number %1.1f \n',Run);
    disp('------------------------------------------------------');
    fprintf('Max point at \t X= %3.4f, \tY= %3.4f \n',Xh,Yh);
    fprintf('Vmin at \t\t X= %3.4f, \tY= %3.4f \n',Xl,Yl);
    fprintf('Core found at \t X= %3.4f, \tY= %3.4f \n\n',Xcore,Ycore);
    fprintf('Core raduis = %3.4f (mm) \n',Rcore);
    fprintf('Vtheta_max = %3.4f (m/s) \n',Maxaz);

%calculate core circle for plotting purposes.
    ang=0:0.01:2*pi; 
    CircX=Rcore*cos(ang);
    CircY=Rcore*sin(ang);
    CircZ(1:length(CircX))=Max;

%Make pretty pictures

%Azimuthal velocity contour map
figure (2); hold on;axis equal; 
    title('Azimuthal velocity component');
    scale=[0:Maxaz/30:Maxaz];
    [cs,h]=contourm(Y,X,Vaz,scale); 
    leg=clegendm(cs,h,1,' m/s');
    contourf(X,Y,Vaz,scale); contour(X,Y,Vaz,scale);
    plot3(Xazh,Yazh,Maxaz,'xw');
    plot3(Xcore,Ycore,Min,'ow');
    plot(Xcore+CircX,Ycore+CircY,':w');
    grid on; %set(gca,'Color',[.01 .01 .01]);
    
%Radial velocity contour map    
figure (3); hold on;axis equal;title('Radial velocity component');
    scale=[-Maxra:Maxra/15:0 Maxra/15:Maxra/15:Maxra];
    [cs,h]=contourm(Y,X,Vra,scale); 
    leg=clegendm(cs,h,1,' m/s');
    contourf(X,Y,Vra,scale); contour(X,Y,Vra,scale);
    plot3(Xazh,Yazh,Maxaz,'xw');
    plot3(Xcore,Ycore,Min,'ow');
    plot(Xcore+CircX,Ycore+CircY,':w');
    grid on; %set(gca,'Color',[0.01 0.01 0.01]);
    
%Axialvelocity contour map    
figure (4); hold on;axis equal;title('Axial velocity component');
    scale=30;%[0 (4/5)*Maxax:Maxax
    [cs,h]=contourm(Y,X,Vmz,scale); 
    leg=clegendm(cs,h,1,' m/s');
    contourf(X,Y,Vmz,scale); contour(X,Y,Vmz,scale);
    plot3(Xazh,Yazh,Maxaz,'xw');
    plot3(Xcore,Ycore,Min,'ow');
    plot(Xcore+CircX,Ycore+CircY,':w');
    grid on; %set(gca,'Color',[0.01 0.01 0.01]);  
    
fclose('all');
    