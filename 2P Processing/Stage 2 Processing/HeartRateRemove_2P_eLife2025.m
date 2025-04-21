function [mpP, fixed_v_noHR]=HeartRateRemove_2P_eLife2025(mpP)
%removes high freqency (>8hz) noise with 2nd order butterworth filter
f_high=5;
fs=mpP.Blood_flow.the_Fs;
[hr_B,hr_A]=butter(2,2*[f_high]/fs,'low');
for k=1:min(size(mpP.Blood_flow.fixed_v))
    fixed_v_noHR(k,:)=filter(hr_B,hr_A,mpP.Blood_flow.fixed_v(k,:));%filter the velocity signal
end
mpP.Blood_flow.fixed_v_noHR=fixed_v_noHR;
end