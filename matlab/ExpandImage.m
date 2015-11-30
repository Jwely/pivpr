%created to interpolate points within an image to double its dimensions
function [Iout]=ExpandImage(Iin)


for i=1:size(Iin,1)
    for j=1:size(Iin,2)
        Image(2*i,2*j)=Iin(i,j);
    end
end

Image=Image(2:end,2:end);

%interpolate surrounding points to create values.
    for i=1:size(Image,1)
        for j=1:size(Image,2)
            if and(rem(i,2)==0,rem(j,2)==0)
                Image(i,j)=(Image(i-1,j-1)+Image(i-1,j+1)+...
                        Image(i+1,j-1)+Image(i+1,j+1))/4;
            elseif and(rem(i,2)==1,rem(j,2)==1)
                %does nothing
            elseif and(rem(i,2)==1,rem(j,2)==0) 
                Image(i,j)= (Image(i,j-1)+Image(i,j+1))/2;                  
            elseif and(rem(i,2)==0,rem(j,2)==1)
                Image(i,j)= (Image(i-1,j)+Image(i+1,j))/2;
            end
        end
    end

Iout=Image;

end