 function fitObj = lfp_findPowerFit(x)

gammafreq = [0.1 600];
in = bz_Filter(x,'passband',gammafreq,'filter','fir1');

Fs = 1250;            % Sampling frequency                    
L = length(in);             % Length of signal

Y = fft(in);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

fAll = Fs*(0:(L/2))/L;
yupper = movmean(P1,1000);
fitObj=fit(fAll(fAll>=6)', yupper(fAll>=6),'power1');

% subplot(2,1,2);cla;
% plot(fitObj,fAll, yupper); xlim([3,200]);
% hold on; 
