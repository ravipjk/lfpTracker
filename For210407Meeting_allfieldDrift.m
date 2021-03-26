clear; clc; 
load('allFieldStruct.mat')
%%

unwrapRelPos_subtrFirst =  allFieldStruct.unwrapRelAng_subtrFirst; 
unwrapEncPos_subtrFirst =  allFieldStruct.unwrapEncAng_subtrFirst; 
gain =  allFieldStruct.gain; 
wrapRelAng =  allFieldStruct.wrapRelAng; 
fieldSize_relAng =  allFieldStruct.fieldSize_relAng; 
stability = allFieldStruct.stability; 
%%
% figure();
diffFieldRaw = [0; diff(wrapRelAng)];
diffField = diffFieldRaw./([0; diff(unwrapEncPos_subtrFirst)]);
% diffField = diffFieldRaw;

diffGain =  [0; diff(gain)];
diffGainIdx = diffGain ~= 0;
diffFieldIdx = abs(diffField) < 5 & abs(diffField) > 0 & stability > 0.8; 
plotIdx = diffGainIdx & diffFieldIdx ;

var1 = gain(plotIdx); 
var2 = diffField(plotIdx); 
[N,c{1}, c{2}] =histcounts2(var1, var2, [5 9]);


subplot(1,2,1); 

normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
normMatDiff(isnan(normMatDiff)) = 0; 
%             normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
            filteredIM_ramp = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
             c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
            c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
            imagesc(c1Interp, c2Interp,imresize(filteredIM_ramp,20)); colormap(jet)
            
xlabel('Gain'); ylabel({'Lap-by-lap drift rate dsitribution', '(degrees drift in LM frame per lap in lab frame'})
 set(gca,'YDir','normal');
 title(['Ramp']); 
 colormap jet;  ylim([-4 4]);  
hold on; plot(xlim, [0 0], 'w--', 'Linewidth',2);  plot([1 1],ylim ,'w--', 'Linewidth',2); hold off; 

% figure();
diffGain =  [0; diff(gain)];
diffGainIdx = diffGain == 0;
% diffFieldIdx = abs(diffField) <30 & abs(diffField) > 0; 
plotIdx = diffGainIdx & diffFieldIdx;


var1 = gain(plotIdx); 
var2 = diffField(plotIdx); 
[N,c{1}, c{2}] =histcounts2(var1, var2, [5 9]);

subplot(1,2,2); 

normMatConst = N'./repmat(sum(N,2), 1,size(N,2))';
normMatConst(isnan(normMatConst)) = 0; 
%     normMatConst = normMatConst./repmat(max(normMatConst,[],2), 1,size(N,1));
    filteredIM_const = imgaussfilt(imgaussfilt(normMatConst,0.6,'FilterDomain','auto'));
     c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
    c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
    imagesc(c1Interp, c2Interp,imresize(filteredIM_const,20)); colormap(jet)
            
xlabel('Gain'); ylabel({'Lap-by-lap drift rate dsitribution', '(degrees drift in LM frame per lap in lab frame'})
 set(gca,'YDir','normal');
             title(['Constant']); 
             colormap jet; ylim([-4 4]);  
hold on; plot(xlim, [0 0], 'w--', 'Linewidth',2);  plot([1 1],ylim, 'w--', 'Linewidth',2); hold off; 

% diffDistri = (imgaussfilt(imresize(filteredIM_const,20))./imgaussfilt(imresize(filteredIM_ramp,20))-1); 
% diffDistriPos = diffDistri;
% diffDistriPos(diffDistriPos<0) = 0;
% diffDistriNeg = diffDistri;
% diffDistriNeg(diffDistriNeg>0) = 0;
% diffDistriNeg = abs(diffDistriNeg);
% 
% figure();
% % subplot(2,1,1);
% imagesc(c1Interp, c2Interp, diffDistriPos); set(gca,'YDir','normal');
% xlabel('Gain'); ylabel('Relative change in distribution'); title('Areas where const has relatively higher activity')
% colorbar
% figure();
% imagesc(c1Interp, c2Interp, diffDistriNeg); set(gca,'YDir','normal');
% xlabel('Gain'); ylabel('Relative change in distribution'); title('Areas where const has relatively lower activity')
% colorbar
%%
diffFieldRaw = [0; diff(wrapRelAng)];
% diffField = diffFieldRaw;
diffField = diffFieldRaw./([0; diff(unwrapEncPos_subtrFirst)]);

diffGain =  [0; diff(gain)];
diffGainIdx = diffGain < 0;
diffFieldIdx = abs(diffField) < 10 & abs(diffField) > 0;
% diffFieldIdx = stability > 0.7;
plotIdx = diffFieldIdx;

var2 = diffField(plotIdx); 
% var1 = [var1; var1]; 
var1 = wrapRelAng(plotIdx); 
% var2 = [var2; var2+360]; 
% var2 = fieldSize_relAng(plotIdx); 


[N,c{1}, c{2}] =histcounts2(var1, var2, [9 50]);
figure();
normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
normMatDiff(isnan(normMatDiff)) = 0; 
% normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
 c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
            
xlabel('Gain'); ylabel('Lap-by-lap drift dsitribution')
 set(gca,'YDir','normal');
             title(['Ramp']); 
             colormap jet;
%%
diffFieldRaw = [0; diff(wrapRelAng)];
diffField = diffFieldRaw;

diffFieldSize = [0;diff(fieldSize_relAng)];

diffGain =  [0; diff(gain)];
diffGainIdx = diffGain ~= 0;
diffFieldIdx = abs(diffField) < 10 & abs(diffField) > 0; 
plotIdx = diffFieldIdx & abs(diffFieldSize) < 10 & gain < 1;

var2 = diffField(plotIdx); 
% var1 = [var1; var1]; 
var1 = diffField([0;plotIdx(1:end-1)]); 
% var2 = [var2; var2+360]; 
% var2 = fieldSize_relAng(plotIdx); 

idx = abs(var1) < 20 & abs(var2) < 20;

[N,c{1}, c{2}] =histcounts2(var1(idx), var2(idx), [10 10]);
figure();
normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
normMatDiff(isnan(normMatDiff)) = 0; 
% normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
 c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
            
xlabel('Gain'); ylabel('Lap-by-lap drift dsitribution')
 set(gca,'YDir','normal');
             title(['Ramp']); 
             colormap jet;
%%
diffFieldRaw = [0; diff(wrapRelAng)];
diffField = diffFieldRaw;

diffFieldSize = [0;diff(fieldSize_relAng)];

diffGain =  [0; diff(gain)];
diffGainIdx = diffGain ~= 0;
diffFieldIdx = abs(diffField) < 10 & abs(diffField) >= 0 & fieldSize_relAng < 50; 
plotIdx = diffFieldIdx & abs(diffFieldSize) < 10 ;

var2 = diffFieldSize(1:end-1); 
% var1 = [var1; var1]; 
var1 = diffFieldSize(2:end); 
% var2 = [var2; var2+360]; 
% var2 = fieldSize_relAng(plotIdx); 

idx = abs(var1) < 20 & abs(var2) < 20;

[N,c{1}, c{2}] =histcounts2(var1(idx), var2(idx), [10 10]);
figure();
normMatDiff = N'./repmat(sum(N,2), 1,size(N,2))';
normMatDiff(isnan(normMatDiff)) = 0; 
% normMatDiff = normMatDiff./repmat(max(normMatDiff,[],2), 1,size(N,1));
filteredIM = imgaussfilt(imgaussfilt(normMatDiff,0.6,'FilterDomain','auto'));
 c1Interp = interp1(1:length(c{1}), c{1}, 1:1/20:length(c{1}));
c2Interp = interp1(1:length(c{2}), c{2}, 1:1/20:length(c{2}));
imagesc(c1Interp, c2Interp,imresize(filteredIM,20)); colormap(jet)
            
xlabel('Gain'); ylabel('Lap-by-lap drift dsitribution')
 set(gca,'YDir','normal');
             title(['Ramp']); 
             colormap jet;