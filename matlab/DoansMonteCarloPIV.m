clear;
clc;

%Calibration Equations
Left_Camera_x_mm=@(xpix,ypix,zmm)-119.955+0.181739.*xpix+0.00114229.*ypix-0.117082.*zmm+5.02817e-006.*xpix.^2-1.9942e-007.*xpix.*ypix-1.29221e-007.*ypix.^2-0.000140106.*xpix.*zmm+2.87984e-006.*ypix.*zmm;
Left_Camera_y_mm=@(xpix,ypix,zmm)-92.2816-0.00345827.*xpix+0.181219.*ypix+0.0776158.*zmm+1.73291e-009.*xpix.^2+4.89664e-006.*xpix.*ypix+2.65745e-008.*ypix.^2-5.31084e-006.*xpix.*zmm-0.000118238.*ypix.*zmm;
Right_Camera_x_mm=@(xpix,ypix,zmm)-124.168+0.191966.*xpix+0.00303101.*ypix+0.289226.*zmm-4.83314e-006.*xpix.^2+7.78692e-008.*xpix.*ypix+1.00813e-007.*ypix.^2-0.000142603.*xpix.*zmm-3.04742e-006.*ypix.*zmm;
Right_Camera_y_mm=@(xpix,ypix,zmm)-90.9543-0.000829306.*xpix+0.183936.*ypix+0.0857377.*zmm+3.32252e-008.*xpix.^2-4.65712e-006.*xpix.*ypix+2.06932e-007.*ypix.^2-1.27237e-005.*xpix.*zmm-0.000124759.*ypix.*zmm;

%Simulation Constants
dt=0.2; % this is actually in pixels?
q=100;
delta_zo=3;
xpixels=1280/4;
ypixels=1024/4;
particle_number=20;
xdisp=1.0;%mm
ydisp=1.0;%mm
zdisp=5.0;%mm
%forimage_num=0:1:39;
[x,y] = meshgrid(1:1:xpixels, 1:1:ypixels);
xmin=min(min(x));
xmax=max(max(x));
ymin=min(min(y));
ymax=max(max(y));

[xlength,ylength]=size(x);

xmm_max1=max(max(Left_Camera_x_mm(x,y,0)));
xmm_min1=min(min(Left_Camera_x_mm(x,y,0)));
ymm_max1=max(max(Left_Camera_y_mm(x,y,0)));
ymm_min1=min(min(Left_Camera_y_mm(x,y,0)));

xo1=(xmm_max1+xmm_min1)/2;
yo1=(ymm_max1+ymm_min1)/2;

xmm_max2=max(max(Right_Camera_x_mm(x,y,0)));
xmm_min2=min(min(Right_Camera_x_mm(x,y,0)));
ymm_max2=max(max(Right_Camera_y_mm(x,y,0)));
ymm_min2=min(min(Right_Camera_y_mm(x,y,0)));

xo2=(xmm_max2+xmm_min2)/2;
yo2=(ymm_max2+ymm_min2)/2;

Io1=@(z)q.*exp(-(z.^2)/((delta_zo.^2)./8));
I1=@(x,y,z)Io1(z).*exp((-(x-xo1).^2-(y-yo1).^2)./((dt.^2)./8));
Io2=@(z)q.*exp(-(z.^2)/((delta_zo.^2)./8));
I2=@(x,y,z)Io2(z).*exp((-(x-xo2).^2-(y-yo2).^2)./((dt.^2)./8));

final_pic1(1:xlength,1:ylength)=0;
final_pic2(1:xlength,1:ylength)=0;
final_pic3(1:xlength,1:ylength)=0;
final_pic4(1:xlength,1:ylength)=0;
%GENERATE PARTICLES
%inner_count=1;

for k=1:1:particle_number;
    %fprintf('Generating Spot %i\n',k);
    pert_x=rand*xmax-xmax/2;
    pert_y=rand*ymax-ymax/2;
    pert_z=rand*delta_zo-delta_zo/2;
    pic1=I1(Left_Camera_x_mm(x+pert_x,y+pert_y,pert_z), Left_Camera_y_mm(x+pert_x,y+pert_y,pert_z),pert_z);
    pic2=I1(Left_Camera_x_mm(x+pert_x,y+pert_y,pert_z+zdisp)-xdisp, Left_Camera_y_mm(x+pert_x,y+pert_y,pert_z+zdisp)-ydisp, pert_z+zdisp);
    pic3=I2(Right_Camera_x_mm(x+pert_x,y+pert_y,pert_z), Right_Camera_y_mm(x+pert_x,y+pert_y,pert_z),pert_z);
    pic4=I2(Right_Camera_x_mm(x+pert_x,y+pert_y,pert_z+zdisp)-xdisp, Right_Camera_y_mm(x+pert_x,y+pert_y,pert_z+zdisp)-ydisp,pert_z+zdisp);
    final_pic1=final_pic1+pic1;
    final_pic2=final_pic2+pic2;
    final_pic3=final_pic3+pic3;
    final_pic4=final_pic4+pic4;
end
%Limit Particle Intensity
for i=1:xlength
    for j=1:ylength
        if final_pic1(i,j)>q
            final_pic1(i,j)=q;
        end
        if final_pic2(i,j)>q
            final_pic2(i,j)=q;
        end
        if final_pic3(i,j)>q
            final_pic3(i,j)=q;
        end
        if final_pic4(i,j)>q
            final_pic4(i,j)=q;
        end
    end
end
%Output Image File
image_num='1';
path='C:\Users\Jeff\Desktop\Github\pivpr\'; cd(path);
filename1=strcat('C:\Users\Jeff\Desktop\Github\pivpr\La',image_num, '.tif');
filename2=strcat('C:\Users\Jeff\Desktop\Github\pivpr\Lb',image_num, '.tif');
filename3=strcat('C:\Users\Jeff\Desktop\Github\pivpr\Ra',image_num, '.tif');
filename4=strcat('C:\Users\Jeff\Desktop\Github\pivpr\Rb',image_num, '.tif');

imwrite(final_pic1,filename1,'tif');
imwrite(final_pic2,filename2,'tif');
imwrite(final_pic3,filename3,'tif');
imwrite(final_pic4,filename4,'tif');
