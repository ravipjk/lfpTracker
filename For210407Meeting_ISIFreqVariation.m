clear; clc; 
load('spikesStruct.mat')
%%

spikeGain =  spikesStruct.gain; 
scaledRelAngle= spikesStruct.scaledRelAngle;
spikeTimes = spikesStruct.spikeTimes_relToFirstSpikeInSweep;
thetaPhase = spikesStruct.thetaPhase;
encAngle = spikesStruct.encAngle_relToFirstSpikeInSweep;
accel = spikesStruct.accel;
vel = spikesStruct.vel;
freqISI = spikesStruct.ISI;
% diffGain = spikesStruct.diffGain;
stability  = spikesStruct.stability;
progress  = spikesStruct.progress;

%%
gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];
figure()
clf; 

n = 20; 
freqEdges = linspace(20,350, n); 
thetaEdges = linspace(0,720, n); 
scaleRelEdges = linspace(0,2, n); 

for i = 1:size(gainWin,1) 
% for i = 3:3    
    scalePosLim = [0 1];
%     thetaLim_init = [240 360];
    thetaLim_init = [0 340];
    thetaLim = thetaLim_init;
    thetaLim = [0 360];
    theta_Window = 90;
    theta_Incr = -.5; 

%     v = VideoWriter(['Gain_', num2str(gainWin(i,1)) '-', num2str(gainWin(i,2)),'_field_windowTheta.avi']);
%     open(v);


    % for relPosLoop = 
    
    
%     while thetaLim(2) >= 0
%         clf; 
%         for i = 5:5     
            idx =  (spikeGain) >= gainWin(i,1) & (spikeGain) <= gainWin(i,2) & ...
                scaledRelAngle < max(scalePosLim) & scaledRelAngle > min(scalePosLim) & ...
                freqISI > 20 & freqISI < 350 & ...
                thetaPhase > max([0 min(thetaLim)]) & thetaPhase < max(thetaLim) & ...
                vel > 5; 

            var1 = scaledRelAngle(idx);
%             var1 = spikeGain(idx);
            var1 = [var1; var1+1];

            var2 = freqISI(idx);
            var2 = [var2; var2];

%            [N,c] =hist3([var1, var2],'CdataMode','auto','Nbins', [20 20]);
           [N,c{1}, c{2}] =histcounts2(var1, var2, scaleRelEdges, freqEdges);
           subplot(2,5,i+5); 
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
        %     normMatDiff = normMatDiff./repmat(sum(normMatDiff,2), 1,size(N,1));
            normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
             c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
            imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
%             imagesc(c{1}, (c{2}),filteredIM ); 
            set(gca,'YDir','normal');
            xlabel('Progress through field');ylabel('Spiking frequency');
%              title({['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2))], ['  Theta Phase: ', num2str(thetaLim(1)) ' - ', num2str(thetaLim(2)) ]}); 
             title({['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2))], sprintf('Theta Phase: %8.2f - %8.2f',thetaLim(1), thetaLim(2))}); 
             title(['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2))]); 
             colormap jet;
             ylim([20 350]);
%         end
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%         thetaLim = thetaLim+theta_Incr; 
%     end
%     close(v);
    %
% 
    thetaLim = [0 360];
    scalePosLim_init = [0 0.30];
%     scalePosLim_init = [0 1];
    scalePosLim = scalePosLim_init; 
    scalePosLim = [0 1]; 
    scalePosLim_Incr = 0.001; 

%     v = VideoWriter(['Gain_', num2str(gainWin(i,1)) '-', num2str(gainWin(i,2)),'_theta_windowField.avi']);
%     open(v);

    % for i = 1:size(gainWin,1) 
    % for relPosLoop = 
%    pause(1)
    
%     while scalePosLim(1) <= 1
%         clf; 
%         for i = 1:1     
            idx =  (spikeGain) >= gainWin(i,1) & (spikeGain) <= gainWin(i,2) & ...
                scaledRelAngle < min([1 max(scalePosLim)]) & scaledRelAngle > min(scalePosLim) &  ...
                freqISI > 25 & freqISI < 400 & ...
                thetaPhase > min(thetaLim) & thetaPhase < max(thetaLim) & ...
                vel > 5; 

            var1 = thetaPhase(idx); 
            var1 = [var1; var1+360];

            var2 = freqISI(idx);
            var2 = [var2; var2];

           [N,c{1}, c{2}] =histcounts2(var1, var2, thetaEdges, freqEdges);
           subplot(2,5,i); 
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
        %     normMatDiff = normMatDiff./repmat(sum(normMatDiff,2), 1,size(N,1));
            normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
            c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
            imTheta = imresize(filteredIM,20); 
            imagesc(c1Interp, c2Interp,imTheta); colormap(jet)
%             imagesc(c{1}, (c{2}),filteredIM ); 
            set(gca,'YDir','normal');
            xlabel('Theta Phase');ylabel('Spiking frequency');
%              title({['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2))],[ '  Field Progress: ', num2str(scalePosLim(1)) ' - ', num2str(scalePosLim(2)) ]}); 
             title({['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2))], sprintf('Field Progress: %8.2f - %8.2f',scalePosLim(1), scalePosLim(2))}); 
             title(['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2))]); 
             colormap jet;
             ylim([20 350]);
%         end
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%         scalePosLim = scalePosLim+scalePosLim_Incr; 
%     end
%     close(v);
% pause(1)
end
%%
% gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];
gainWin = [0 0.99;1 1; 1.01 2];
gainWin = [0 0.6; 1 1; 1.4 2];
% gainWin = [1 1];
% gainWin = [0 0.99];
% gainWin = [1.01 2];

diffPhase = [0; diff(thetaPhase)]; 
figure(7); 
clf;
% subplot(1,2,1);
for i = 1:size(gainWin,1) 
  
% for i = 3:3    
    scalePosLim = [0 1];
%     thetaLim_init = [0 360];
%     thetaLim_init = [120 270];
%     thetaLim_init = [270 360];
%     thetaLim = thetaLim_init;
%     thetaLim = [0  220]; 
%         thetaLim = [60 360];
%     thetaLim = [220 360];
    thetaLim = [0 360];
%     theta_Window = 90;
%     theta_Incr = 1; 
% progressLim = [0.0 0.1];

%     v = VideoWriter(['spk_accel_ThetaFreq.avi']);
%     open(v);


    % for relPosLoop = 
    
    progressLim = [0.0 1]; 
%     gainWin = [0.2 1.8]; 
%     while gainWin(2) < 1.8
    while progressLim(2) <= 1
%         clf; 
%         for i = 5:5     
            idx =  (spikeGain) >= gainWin(i,1) & (spikeGain) <= gainWin(i,2) & ...
                scaledRelAngle < scalePosLim(2) & scaledRelAngle > scalePosLim(1) & ...
                freqISI >6 & freqISI <12 & ...
                thetaPhase > thetaLim(1) & thetaPhase < thetaLim(2) & ...
                vel > 5 & progress >= progressLim(1) & progress <= progressLim(2) & vel < 50 & abs(accel) < 200 & stability<= 1;% & spikeGain==1; %progressLim
%             & spikeTimes/1e6<8 ; 
%             var1 = accel(idx);
%             var1 = spikeGain(idx).*vel(idx);
            var1 = vel(idx);
%                 var1 = scaledRelAngle(idx); 
%                var1 = spikeGain(idx);
%                var1 = progress(idx)*100;
% %             var1 = stability(idx);
%             var1 = thetaPhase(idx);
%                 var1 = diffPhase(idx);
%             var1 = thetaPhase(idx);
%             var1 =  encAngle(idx);
%             var1 = [var1; var1+360];

%             var2 = scaledRelAngle(idx); 
            var2 = (freqISI(idx));
%             var2 = thetaPhase(idx); 
%             var2 = [var2; var2];

%            [N,c] =hist3([var1, var2],'CdataMode','auto','Nbins', [20 20]);
           [N,c{1}, c{2}] =histcounts2(var1, var2, [10 10]);
           subplot(1,size(gainWin,1),i); 
           
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
            normMatDiff(isnan(normMatDiff)) = 0; 
%             normMatDiff = normMatDiff./repmat(sum(normMatDiff,2), 1,size(N,1));
%             normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
             c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
%             subplot(1,2,1);
            imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
%             imagesc(c{1}, (c{2}),filteredIM ); 
            set(gca,'YDir','normal');
            xlabel('Speed');ylabel('Spiking frequency');
%              title(['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2)), '  Theta Phase: ', num2str(thetaLim(1)) ' - ', num2str(thetaLim(2)) ]); 
             title(['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2))]); 
             colormap jet;
             hold on; plot([0.0 0.0], ylim, 'w--', 'Linewidth',3); hold off; 
%              title(sprintf('Progress: %8.2f - %8.2f',progressLim(1,1), progressLim(1,2))); 
%             title(sprintf('Velocity bin: %8.2f - %8.2f',5, 25)); 
             ylim([6.5 12]); 
%             xlim([-100 100]);
%         end
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%         thetaLim = thetaLim+theta_Incr; 
%         gainWin = gainWin+0.01;     
        progressLim = progressLim + 0.01; 
%         pause(0.5)
    end
%     close(v);
end
%% 
% progress through epochs
gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];
figure()
clf; 

% n = 30; 
freqLim = [20 200]
freqEdges = linspace(freqLim(1),freqLim(2), 10); 
thetaEdges = linspace(0,720, 20); 
scaleRelEdges = linspace(0,2, n); 

% for i = 1:size(gainWin,1) 
% for i = 3:3    
    scalePosLim = [0 1];
%     thetaLim_init = [240 360];
    thetaLim_init = [0 360];
    thetaLim = thetaLim_init;
    thetaLim = [0 360];
    theta_Window = 90;
    theta_Incr = -.5; 

%     v = VideoWriter(['Gain_', num2str(gainWin(i,1)) '-', num2str(gainWin(i,2)),'_field_windowTheta.avi']);
%     open(v);


    % for relPosLoop = 
    
    
%     while thetaLim(2) >= 0
%         clf; 
%         for i = 5:5   
    progress_init = [0 0.3]; 
    
    while progress_init(1) <= 0.9
        mean(progress_init)*100
            idx = scaledRelAngle < max(scalePosLim) & scaledRelAngle > min(scalePosLim) & ...
                freqISI > freqLim(1) & freqISI < freqLim(2) & ...
                thetaPhase > max([0 min(thetaLim)]) & thetaPhase < max(thetaLim) & ...
                vel > 5 & progress >= progress_init(1) &  progress <= progress_init(2) & spikeGain == 1; 

            var1 = thetaPhase(idx);
%             var1 = spikeGain(idx);
            var1 = [var1; var1+360];

            var2 = freqISI(idx);
            var2 = [var2; var2];

%            [N,c] =hist3([var1, var2],'CdataMode','auto','Nbins', [20 20]);
           [N,c{1}, c{2}] =histcounts2(var1, var2, thetaEdges, freqEdges);
%            subplot(1,5,i); 
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
            normMatDiff(isnan(normMatDiff)) = 0; 
        %     normMatDiff = normMatDiff./repmat(sum(normMatDiff,2), 1,size(N,1));
            normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
             c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
            imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
%             imagesc(c{1}, (c{2}),filteredIM ); 
            set(gca,'YDir','normal');
            xlabel('Progress through field');ylabel('Spiking frequency');
%              title({['Gain: ', num2str(gainWin(i,1)) ' - ', num2str(gainWin(i,2))], ['  Theta Phase: ', num2str(thetaLim(1)) ' - ', num2str(thetaLim(2)) ]}); 
             colormap jet;
             ylim(freqLim);
%         end
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%         thetaLim = thetaLim+theta_Incr; 
        pause(0.1)
            progress_init = progress_init + 0.01; 
    end
%     close(v);
    %
% 
 
% end

%%
clc;
% gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];
gainWin = [0 2];
% X = [thetaPhase scaledRelAngle spikeGain spikeTimes encAngle vel accel diffGain];
X = [spikeGain];
scalePosLim = [0 1];
thetaLim= [0 360];
for i = 1:size(gainWin,1)   
    idx =  (spikeGain) >= gainWin(i,1) & (spikeGain) <= gainWin(i,2) & ...
        scaledRelAngle < scalePosLim(2) & scaledRelAngle > scalePosLim(1) & ...
        freqISI > 25 & freqISI < 350 & ...
        thetaPhase > thetaLim(1) & thetaPhase < thetaLim(2) & ...
        vel > 5; 
    
    mdl = fitglm(X(idx,:),freqISI(idx))
    mdl1 = step(mdl,'NSteps',5)
end
%%
%%

gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];

% gainWin = [0 0.8;1 1; 1.2 2];
freqLim = [20 120]; 
% for i = 1:size(gainWin,1) 
n = 10; 
freqEdges = linspace(freqLim(1),freqLim(2), n); 
thetaEdges1 = linspace(0-90,270+90, 180); 
thetaEdges2 = linspace(90-90,360+90, 180); 
scaleRelEdges = linspace(0,1, n); 
D_scaledRel = [];
binNum = 3; 

for i = binNum:binNum  
    scalePosLim = [0 1];
    for j = 1:length(thetaEdges1)
            idx =  (spikeGain) >= gainWin(i,1) & (spikeGain) <= gainWin(i,2) & ...
                scaledRelAngle < max(scalePosLim) & scaledRelAngle > min(scalePosLim) & ...
                freqISI > freqLim(1) & freqISI < freqLim(2) & ...
                thetaPhase > thetaEdges1(j) & thetaPhase < thetaEdges2(j) & ...
                vel > 5; 
            var1 = scaledRelAngle(idx);
            var2 = freqISI(idx);
           [N,c{1}, c{2}] =histcounts2(var1, var2, scaleRelEdges, freqEdges);
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
            normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
            c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
%             imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)        
        D_scaledRel(:,j,:) = imresize(filteredIM,20); 
    end
end

x = c1Interp(1:end-1); 
y = (thetaEdges1+thetaEdges2)/2;
z = c2Interp(1:end-1); 

[X_rel,Y_rel,Z_rel]=meshgrid(x,y,z);


% gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];

% for i = 1:size(gainWin,1) 
n = 10; 
freqEdges = linspace(freqLim(1),freqLim(2), n); 
thetaEdges = linspace(0-90,360+90, n); 
scaleRelEdges1 = linspace(-0.3,1.0, 180); 
scaleRelEdges2 = linspace(0.0,1.3, 180); 
D_theta = [];

for i = binNum:binNum    
    scalePosLim = [0 1];
    for j = 1:length(scaleRelEdges1)
            idx =  (spikeGain) >= gainWin(i,1) & (spikeGain) <= gainWin(i,2) & ...
                scaledRelAngle < scaleRelEdges2(j) & scaledRelAngle > scaleRelEdges1(j) & ...
                freqISI > freqLim(1) & freqISI < freqLim(2) & ...
                thetaPhase > 0 & thetaPhase < 360 & ...
                vel > 5; 
            var1 = thetaPhase(idx);
            var2 = freqISI(idx);
           [N,c{1}, c{2}] =histcounts2(var1, var2, thetaEdges, freqEdges);
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
            normMatDiff(isnan(normMatDiff)) = 0; 
            normMatDiff(isinf(normMatDiff)) = 0; 
            normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
            c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
%             imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)  
%             clf; imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet);  set(gca,'YDir','normal');   
%             pause(0.1);
        D_theta(j,:,:) = (imresize((filteredIM),20)); 
    end
end

x = (scaleRelEdges1+scaleRelEdges2)/2;
y = c1Interp(1:end-1); 
z = c2Interp(1:end-1); 

[X_th,Y_th,Z_th]=meshgrid(x,y,z);


% gainWin = [0 0.6;0.6 0.99;1 1; 1.01 1.4; 1.4 2];

% for i = 1:size(gainWin,1) 
n = 10; 
freqEdges1 = linspace(-20,350, 180); 
freqEdges2 = linspace(20,390, 180); 
thetaEdges = linspace(0,360, n); 
scaleRelEdges = linspace(0,1, n); 
% scaleRelEdges2 = linspace(0.3,1, 180); 
D_isi = [];

for i =binNum:binNum  
    scalePosLim = [0 1];
    for j = 1:length(freqEdges1)
            idx =  (spikeGain) >= gainWin(i,1) & (spikeGain) <= gainWin(i,2) & ...
                scaledRelAngle < 1 & scaledRelAngle > 0 & ...
                freqISI > freqEdges1(j) & freqISI < freqEdges2(j) & ...
                thetaPhase > 0 & thetaPhase < 360 & ...
                vel > 5; 
            var1 = scaledRelAngle(idx);
            var2 = thetaPhase(idx);
           [N,c{1}, c{2}] =histcounts2(var1, var2, scaleRelEdges, thetaEdges);
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
            normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
            c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
%             imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)        
        D_isi(:,:,j) = imresize(filteredIM,20); 
    end
end


x = c1Interp(1:end-1); 
y = c2Interp(1:end-1); 
z = (freqEdges1+freqEdges2)/2;

[X_isi,Y_isi,Z_isi]=meshgrid(x,y,z);

%
n = 10; 
progressEdges1 = linspace(0-0.3,0.7+0.3, 180); 
progressEdges2 = linspace(0.3-0.3,1+0.3, 180);  
thetaEdges = linspace(0,360, n); 
freqEdges = linspace(freqLim(1),freqLim(2), n); 
% scaleRelEdges2 = linspace(0.3,1, 180); 
D_progress = [];

for i =binNum:binNum  
    scalePosLim = [0 1];
    for j = 1:length(progressEdges1)
            idx =  (spikeGain) >= gainWin(i,1) & (spikeGain) <= gainWin(i,2) & ...
                scaledRelAngle < 1 & scaledRelAngle > 0 & ...
                freqISI > freqLim(1) & freqISI < freqLim(2) & ...
                thetaPhase > 0 & thetaPhase < 360 & ...
                vel > 5 & progress > progressEdges1(j) & progress < progressEdges2(j) ; 
            var1 = thetaPhase(idx);
            var2 = freqISI(idx);
           [N,c{1}, c{2}] =histcounts2(var1, var2, thetaEdges, freqEdges);
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
            normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
            c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
%             clf; imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet);set(gca,'YDir','normal');    
%             pause(0.1); 
        D_progress(j,:,:) = imresize(filteredIM,20); 
    end
end


x = (progressEdges1+progressEdges2)/2;
y = c1Interp(1:end-1); 
z = c2Interp(1:end-1); 

[X_progress,Y_progress,Z_progress]=meshgrid(x,y,z);
%%

figure(40); 
% subplot(1,5,5)
isovalue = 0.82;
surf1 = isosurface(X_th,Y_th,Z_th,D_theta,isovalue);
p1 = patch(surf1);
% isonormals(x,y,z,D,p1);
set(p1,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface

 xlabel('Scaled Rel Pos');xlim([0 1]); zlabel('ISI Freq');zlim(freqLim); ylabel('Theta phase');ylim([0 360]); 
camlight; lighting gouraud; grid on; 
%%
figure(50); 
isovalue = 0.8;
surf1 = isosurface(X_rel,Y_rel,Z_rel,D_scaledRel,isovalue);
p1 = patch(surf1);
% isonormals(x,y,z,D,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface


 xlabel('Scaled Rel Pos');xlim([0 1]); zlabel('ISI Freq');zlim(freqLim); ylabel('Theta phase');ylim([0 360]); 
camlight; lighting gouraud; grid on; 
%%
figure(60); 
isovalue = 0.75;
surf1 = isosurface(X_isi,Y_isi,Z_isi,D_isi,isovalue);
p1 = patch(surf1);
% isonormals(x,y,z,D,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface


 xlabel('Scaled Rel Pos');xlim([0 1]); zlabel('ISI Freq');zlim([20 350]); ylabel('Theta phase');ylim([0 360]); 
camlight; lighting gouraud; grid on; 
%%
figure(3); 
isovalue = 0.8;
surf1 = isosurface(X_rel,Y_rel,Z_rel,D_scaledRel.*D_theta,isovalue);
p1 = patch(surf1);
% isonormals(x,y,z,D,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface


 xlabel('Scaled Rel Pos');xlim([0 1]); zlabel('ISI Freq');zlim([20 350]); ylabel('Theta phase');ylim([0 360]); 
camlight; lighting gouraud; grid on; 

%%
figure(70); 
% subplot(1,5,5)
isovalue = 0.82;
surf1 = isosurface(X_progress,Y_progress,Z_progress,D_progress,isovalue);
p1 = patch(surf1);
% isonormals(x,y,z,D,p1);
set(p1,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.5); % set the color, mesh and transparency level of the surface

 xlabel('Progress');xlim([0 1]); zlabel('ISI Freq');zlim(freqLim); ylabel('Theta phase');ylim([0 360]); 
camlight; lighting gouraud; grid on;