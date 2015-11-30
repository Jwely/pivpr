% Find displacement between 2 grayscale images by using cross-correlation
% output is in units of pixels, so use concurently with 'ExpandImage'
% function if desired displacement accuracy is greater than a pixel.


function [yoffset,xoffset]=findoff(unreg , base)

c = normxcorr2(unreg,base);

[~, imax] = max(abs(c(:)));
[ypeak, xpeak] = ind2sub(size(c),imax(1));

corr_offset = [(xpeak-(size(c,2)+1)/2) (ypeak-(size(c,1)+1)/2)];

offset = corr_offset;
xoffset = offset(1);
yoffset = offset(2);
