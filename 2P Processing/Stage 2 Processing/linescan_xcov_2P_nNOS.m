function [maxspot,xc_image,xshift]=linescan_xcov_2P(MScanData)
%maxspot-amplitude of the cross correlation peak.
%xcorr_image- power-spectra normalized cross correlation image
%xshift- shift of the peak away from the center, in pixels
x_spread=round(max(1,.5/MScanData.notes.Xfactor));%spatial gaussian with 0.5um std
t_spread=round(MScanData.notes.Tfactor/10);%temporal gaussian with 10ms std
the_kernel=gaussian2d_2P(x_spread,t_spread,3*x_spread,3*t_spread);
theimage=double(MScanData.Blood_flow.Image');%   DOUBLE CHECK IS THIS REALLY THE RIGHT IMAGE????


nlines=size(theimage,2);
npoints=size(theimage,1);
average_line=mean(theimage,2);
xc_image=zeros(size(theimage));
for t=1:nlines
    theimage(:,t)=theimage(:,t)-average_line;%subtract out the background
end
for t=3:nlines
    %take the convolution of the two line-scans normalized by the joint
    %power spectrums
    xc_image(:,t)=ifft(fft(theimage(:,t)).*conj(fft(theimage(:,t-1)))./sqrt(abs(fft(theimage(:,t))).*abs(fft(theimage(:,t-1)))));
end
%convolve with a gaussian space-time kernel to average velocity
xc_image=conv2(fftshift(xc_image,1),the_kernel,'same');%fftshift puts the lower absolute values for the velocities togetehr in the middel of the matrix,then we filter with a matched kernel
[~,maxspot]=max(xc_image(:,:));%find the peak in the cross correlation
xshift=(round(npoints/2)-maxspot);
end