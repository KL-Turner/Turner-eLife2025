function [m2error] = DoubleGammaConvolve_IOS_nNOS(x,inp,outp,Fs,HRFDur)
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Dec 2011
%   Version 2 - Updated Jan 2012
%   SUMMARY: Creates a gamma function from a set of provided parameters,
%   convolves the gamma function with the variable called "input" and
%   calculates the mean squared error between the result of the convolution
%   and the variable called "output".
%_______________________________________________________________
%   INPUTS:
%                   x - inital gamma function parameter values as an array: 
%                       [x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8)]. 
%                           x(1) = amplitude of the gamma,
%                           x(2) = time to peak, 
%                           x(3) = peak width @ 75% max 
%                           x(4) = dc offset of convolution result
%                           (optional)
%
%                           x(5) = amplitude of the negative gamma,
%                           x(6) = time to peak, 
%                           x(7) = peak width @ 75% max 
%                           x(8) = dc offset of convolution result
%                           (optional)
%                   inp - array of data to be convolved with the gamma
%                   function.
%                   outp - array of data to be compared to the result of
%                   the convolution by mean squared error.
%                   Fs - the sampling frequency for the data
%_______________________________________________________________
%   OUTPUTS:
%                   m2error - The mean squared error between the vector
%                   called "output" and the result of convolution between
%                   vector called "input" and gamma function.
%_______________________________________________________________
%   REQUIRED SCRIPTS: None
%_______________________________________________________________
%   CALLED BY: 
%               IRFsearch_gamma.m
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________
                     
% Create Gamma Function 
t = 0:1/Fs:HRFDur;
a = ((x(2)/x(3))^2*8*log10(2));
a_undr = ((x(6)/x(7))^2*8*log10(2)); % second alpha value used to estimate undershoot 04/30/20 KWG
beta = ((x(3)^2)/x(2)/8/log10(2));
beta_undr = ((x(7)^2)/x(6)/8/log10(2)); % second beta value used to estimate undershoot 04/30/20 KWG
gamma = x(1)*(t/x(2)).^a.*exp((t-x(2))/(-1*beta)); % this is your original positive gamma function from Winder et al 2016.
gamma_undr = x(5)*(t/x(6)).^a_undr.*exp((t-x(6))/(-1*beta_undr)); % second gamma function using values x(5:8)
gamma = gamma - gamma_undr; %minimize the function that is the first positive gamma function and a second negative gamma function to model positive and negative hemodynamic response.
% Perform Convolution
convl = conv(inp, gamma);
% Add DC offset to result
% convout = convl(1:length(outp))+x(4);
convout = convl(1:length(outp));
% Calculate mean squared error
m2_ind = round(Fs:length(convout)-Fs);
r2 = 1- sum((outp(m2_ind) - convout(m2_ind)).^2)/sum((outp(m2_ind) - mean(outp(m2_ind))).^2); %#ok<NASGU>
m2error = mean((convout(m2_ind) - outp(m2_ind)).^2);
end
