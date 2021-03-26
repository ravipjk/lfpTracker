clear; clc; 
load rippleStruct.mat
%%

idx =   rippleStruct.vel > 5 & rippleStruct.vel < 50 & (diff(rippleStruct.timestamps, [],2) < 0.05)  & (diff(rippleStruct.timestamps, [],2) > 0.03);

figure();
% var1= rippleStruct.gain(idx);
var1 = rippleStruct.vel(idx);

% var2= rippleStruct.theta(idx);
% % var1 = rippleStruct.vel(idx);
var2 = diff(rippleStruct.timestamps(idx,:), [],2);
% var1 = rippleStruct.gain(idx);

%             var1 = spikeGain(idx);
%             var1 = [var1; var1+1];

           
%             var2 = [var2; var2];

%            [N,c] =hist3([var1, var2],'CdataMode','auto','Nbins', [20 20]);
           [N,c{1}, c{2}] =histcounts2(var1, var2, [10 10]);
%            subplot(2,5,i+5); 
            normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
        %     normMatDiff = normMatDiff./repmat(sum(normMatDiff,2), 1,size(N,1));
%             normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
             c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
            imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
             set(gca,'YDir','normal')