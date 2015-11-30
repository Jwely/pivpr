%==========================================================================
%           [VecFile2d]= ImageTo2dVectors(path,Sector_size,order)
%========================================================================== 
% This Function analyzes PIV images and creates 2d vector fields.

%========================================================================== 
%                           Function Syntax
%========================================================================== 
function [VecFile2d]= ImageTo2dVectors(path,,calibration,Sector_size,order)

figs = 'off';       %no figures, except to generate example plots.

La=imread(strcat(path,'La.tif'));
Lb=imread(strcat(path,'Lb.tif'));
Ra=imread(strcat(path,'La.tif'));
Rb=imread(strcat(path,'La.tif'));




%upsample the image data to allow higher precision displacement resolution
if order~=0
    for i=1:order
        A=expandimage(A);
        B=expandimage(B);
    end
end
       
[xmove,ymove]=findoff(B,A);
xmove/(2^order)
ymove/(2^order)


%plotting a heatmap of the two refined variables.      
if figs(1:2)=='on'
figure (1);
    ha=imagesc(A); axis off; grid on; grid minor;
    colormap(cbrewer('seq','Reds',32,'cubic'));
    title('Colorized Monochromatic PIV image sector');
figure (2);
    hb=imagesc(B); axis off; grid on; grid minor;
    colormap(cbrewer('seq','Blues',32,'cubic'));
    title('Colorized Monochromatic PIV image sector');
    
    print(figure (1),'PIV example figureA','-djpeg85','-r800');
    print(figure (2),'PIV example figureB','-djpeg85','-r800');
end    
end