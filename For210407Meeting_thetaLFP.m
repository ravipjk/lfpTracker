clear; clc; 
load thetaStruct.mat
%%
 th_freq = movmean(thetaStruct.freq, 1);
th_amp = movmean(thetaStruct.amp, 1);
th_gain = thetaStruct.gain;
th_vel = thetaStruct.vel;
th_accel = thetaStruct.accel;
th_sym = thetaStruct.sym;
th_ratNum = thetaStruct.ratNum;
th_autoLM = thetaStruct.autoControl;
th_progress = thetaStruct.progress; 

%%
gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];
gainWin = [0 1.8];
% gainWin = [1 1];

figure(8); clf;
for i = 1:size(gainWin,1) 
        clf; 
            idx =  (th_gain) >= gainWin(i,1) & (th_gain) <= gainWin(i,2) & ...
                th_vel > 25  & th_vel < 50 & abs(th_accel) < 200 & th_freq > 6 & th_freq < 8.5 & th_ratNum == 515;% & th_amp > 9000 & th_amp < 18000;
%             var1 = th_accel(idx);
            var1 = th_vel(idx);
%             var1 = th_gain(idx).*th_vel(idx);
%             var1 = th_progress(idx);
            var1 = th_gain(idx);
            var2 = th_freq(idx);


           [N,c{1}, c{2}] =histcounts2(var1, var2, [10 20]);
           subplot(1,size(gainWin,1),i); 
           
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
            normMatDiff(isnan(normMatDiff)) = 0; 
%             normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
             c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
            imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
            set(gca,'YDir','normal');
            xlabel('Stability');ylabel('Theta frequency');
             colormap jet;
end

%%
%%
gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];
% gainWin = [1.1 2];
% gainWin = [0 0.9];
gainWin = [1 1];

figure(8); clf;
gainWin = [0.2 0.6]; 
while gainWin(2) < 1.8
        clf; 
            idx =  (th_gain) >= gainWin(i,1) & (th_gain) <= gainWin(i,2) & ...
                th_vel > 5  & th_vel < 50 & abs(th_accel) < 200 & th_freq > 6 & th_freq < 11 & th_ratNum == 692 & th_amp > 9000 & th_amp < 18000 & th_accel < 0;
%             var1 = th_accel(idx);
            var1 = th_vel(idx);
%             var1 = th_gain(idx).*th_vel(idx);
%             var1 = th_progress(idx);

            var2 = th_freq(idx);
%             var2 = th_amp(idx);


           [N,c{1}, c{2}] =histcounts2(var1, var2, [20 20]);
           subplot(1,size(gainWin,1),i); 
           
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
            normMatDiff(isnan(normMatDiff)) = 0; 
%             normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
             c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
            imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
            set(gca,'YDir','normal');
            xlabel('Stability');ylabel('Theta frequency');
            title(sprintf('Gain bin: %8.2f - %8.2f',gainWin(1), gainWin(2))); 
             colormap jet; 
        ylim([6.5 9]); 
                xlim([5 50]);
        gainWin = gainWin+0.01;     
%         progressLim = progressLim + 0.01; 
        pause(0.5)
end