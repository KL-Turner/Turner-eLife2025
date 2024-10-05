
function [data]=gaussian2d_2P_nNOS(xstd,ystd,xsize,ysize)
%2-d gaussian kernel for filtering
data=zeros(xsize-1,ysize-1);
x0=xsize/2;
y0=ysize/2;
for x=1:xsize-1
    for y=1:ysize-1
        data(x,y)=exp(-((x-x0)^2)/(2*xstd)-((y-y0)^2)/(2*ystd));
    end
end
data=data/sum(data(:));
end