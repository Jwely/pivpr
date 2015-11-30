%script to provide an example of how the PIV system extracts velocities
close all force; clear all;clc;

order=3;


%use 56010 for examples also use +156 for examples.

L1in=imread('E:\Data2\Ely_May28th\Image\Dual\63\Ely_May28th63010La.tif');
L2in=imread('E:\Data2\Ely_May28th\Image\Dual\63\Ely_May28th63010Lb.tif');


figs='on';
W=size(L1in,1)/2+300;
L=size(L1in,2)/2+300;
u=16;

A=L1in((W-u):(W+u-1),(L-u):(L+u-1));
A= A-mean2(A);

    
B=L2in((W-u):(W+u-1),(L-u):(L+u-1));
B= B-mean2(B);

%upsample the image data to allow higher precision displacement resolution
tic();

if order~=0
for i=1:order
    A=expandimage(A);
    B=expandimage(B);
end
end

a=fft2(A); b=fft2(B);
C=real(ifft2((b.*conj(a))));
Ctrim=C(1:end/2,1:end/2);

m=max(max(Ctrim));
[y,x]=find(C==m);
xfft=(x(1)-1)/((2^order));
yfft=(y(1)-1)/((2^order));

%plotting a heatmap of the two refined variables.      

if figs(1:2)=='on'
figure (1); hold on; 
    colormap(cbrewer('seq','Reds',32,'cubic')); grid off;
    contourf(A,16);contour(A,16); axis off;
    title('Colorized Monochromatic PIV image sector');
figure (2); hold on;
    colormap(cbrewer('seq','Blues',32,'cubic')); grid off;
    contourf(B,16);contour(B,16); axis off;
    title('Colorized Monochromatic PIV image sector');
    
figure(3); hold on;
    colormap(cbrewer('div','RdYlGn',32,'cubic')); grid off;
    contourf(Ctrim,16); contour(Ctrim,16);colorbar
    title('2d InverseFFT(F1*conj(F2)) Map');grid minor; box on;
    
    print(figure (1),strcat('PIVe figA order',num2str(order)),'-djpeg85','-r800');
    print(figure (2),strcat('PIVe figB order',num2str(order)),'-djpeg85','-r800');
    print(figure (3),strcat('PIVe FFT order',num2str(order)),'-djpeg85','-r800');
end    

toc()
fprintf('%1.4f \n',xfft);
fprintf('%1.4f \n\n',yfft);

