clear; clc; 
load('F:\dome_data\Rat638\160617_Rat638-15\analyzed\m1_v2_specGain.mat')
load('F:\dome_data\Rat638\160617_Rat638-15\analyzed\m1_rosdata.mat')
load('F:\dome_data\Rat638\160617_Rat638-15\analyzed\m1_hipp_hipp_clustExtended.mat')


%%

figure(); 
try
        landMsgIdx  = find(strcmp({rosdata.domeVisMsgs.type},'landmarks'));
        b           = [rosdata.domeVisMsgs.visible];
        landVis     = b(landMsgIdx);
        landOffTime = rosdata.domeVisTimes(landMsgIdx(find(landVis == 0,1)));
    catch
        idx         = find(strcmp({rosdata.domeEvMsgs.name},'landmarks'));
        landOffTime = rosdata.domeEvTimes(idx(find([rosdata.domeEvMsgs(idx).  value] == 0,1)));
end

if ~isempty(landOffTime)
    [~,landOffIdx]  = min(abs(rosdata.encTimes-landOffTime));
    landOffEnc      = rosdata.encAngle(landOffIdx);
else
    landOffEnc = inf;
end
    
k = 1; 
listClust = [2 10]; 
% listClust = [1:11]; 
legendCell = {}; 
for iLoop = 1:length(listClust)
    i = listClust(iLoop); 
   uniqueTags = unique(clustExtended(i).fieldTag); 
   uniqueTags(isnan(uniqueTags)) = [];
   
   if size(uniqueTags,2) == 1
       i
       fieldVec = [[clustExtended(i).fields.encAngle1]' [clustExtended(i).fields.encAngle2]'];
       clustGainEnc = mean(fieldVec(2:end,:),2); 
       clustGain = 360./mean(diff(fieldVec),2); 
       TF = isoutlier(clustGain,'movmedian',5); 
       outlierReplacementGain = interp1(specGain.encAngle,specGain.popGain,clustGainEnc(TF));
       ratioGain = outlierReplacementGain./clustGain(TF); 
       ratioInt = round(ratioGain);
       ratioRem = abs(ratioGain - ratioInt);
       ratioRemThresh = (ratioRem < 0.1 | ratioRem > 0.9) &  ismember(ratioInt,[2 3 4]);
       clustGain(TF) = clustGain(TF).*ratioInt; 
       clustGain = filloutliers(clustGain,'center');  
       plot(clustGainEnc/360, (clustGain));hold on; 
       
%        fieldVec = [clustExtended(i).fields.encAnglePeak]';
%        plot(fieldVec(2:end), filloutliers(360./diff(fieldVec),'linear'));

%         clustGainEnc = mean(fieldVec(2:end,:),2); 
        
%         clustGain = 360./(diff(mean(fieldVec,2))); 
%         TF = isoutlier(clustGain); 
%         outlierReplacementGain = interp1(specGain.encAngle,specGain.popGain,clustGainEnc(TF)); 
%         clustGain(TF) = outlierReplacementGain; 
% 
%         idx = max(find(specGain.encAngle<clustGainEnc(2))); 
% 
%         combinedEnc = [specGain.encAngle(1:idx); clustGainEnc];
%         combinedGain = [specGain.popGain(1:idx);clustGain];
%         
% %         [combinedEnc, idx] = sort([specGain.encAngle; clustGainEnc]);
% %         combinedGain = [specGain.popGain;clustGain];
% %         combinedGain = movmean(combinedGain(idx), 10); 
% 
%         hippGainEnc_combined = interp1(combinedEnc,combinedGain,rosdata.encAngle,'linear','extrap');
% %         rosHippAngleMon_combined = cumtrapz(rosdata.encAngle,abs(hippGainEnc_combined)) +rosdata.encAngle(1);
%         rosHippAngleMon_combined = rosdata.encAngle(1)+ cumsum([0 diff(rosdata.encAngle)].*hippGainEnc_combined);
%         
%         hippGainEnc = interp1(specGain.encAngle,specGain.popGain,rosdata.encAngle,'linear','extrap');
%         rosHippAngleMon = cumtrapz(rosdata.encAngle,abs(hippGainEnc)) +rosdata.encAngle(1);
%         rosHippAngleMon = rosdata.encAngle(1)+ cumsum([0 diff(rosdata.encAngle)].*hippGainEnc);
% % 
% %         hippAngle_combined      = interp1(rosdata.encTimes, rosHippAngleMon_combined, clustExtended(i).ts,'linear','extrap'); 
%         hippAngle      = interp1(rosdata.encTimes, rosHippAngleMon, clustExtended(i).ts,'linear','extrap'); 
% % 
%         velIdx = clustExtended(i).vel > 5; 
% %         subplot(2,1,1);
%         plot(clustExtended(i).encAngle(velIdx)/360, mod(hippAngle(velIdx)+180, 360), '.');
%         ylim([0 360]);
%         subplot(2,1,2); 
%         plot(clustExtended(i).encAngle(velIdx)/360, mod(hippAngle_combined(velIdx)+180, 360), '.');
%         ylim([0 360]);
%         pause();
        legendCell{k} =  clustExtended(i).name; 
        k = k+1; 
   end
   
end
legendCell{k} = 'Spectral gain'; 
plot(specGain.encAngle/360,specGain.popGain, 'k--', 'LineWidth', 2); xlabel('Laps in lab frame'); ylabel('Estimated Gain'); 
plot([landOffEnc/360 landOffEnc/360], ylim, 'k');
legendCell{k+1} = 'Landmark Off'; 
legend(legendCell);
%%
figure()
velIdx = clustExtended(i).vel > 5; 
subplot(3,1,1)
plot(clustExtended(i).encAngle(velIdx)/360, mod(clustExtended(i).encAngle(velIdx)+180, 360), '.');hold on; 
plot([landOffEnc/360 landOffEnc/360], ylim, 'k'); ylim([0 360]);xlim([0 135]); xlabel('Laps in lab frame'); ylabel('Angle in Lab Frame');
subplot(3,1,2)
gainVec = 1-rosdata.gain; 
gainVec(rosdata.gainTimes > landOffTime) = gainVec(max(find(rosdata.gainTimes < landOffTime)));
gainEnc = interp1(rosdata.gainTimes,gainVec,rosdata.encTimes,'linear','extrap');
% rosRelAngleMon = cumtrapz(rosdata.encAngle,abs(1-gainEnc));
rosRelAngleMon = rosdata.encAngle(1)+ cumsum([0 diff(rosdata.encAngle)].*gainEnc);
relAngle      = interp1(rosdata.encTimes, rosRelAngleMon, clustExtended(i).ts,'linear','extrap'); 
plot(clustExtended(i).encAngle(velIdx)/360, mod(relAngle(velIdx)+180, 360), '.'); hold on; 
plot([landOffEnc/360 landOffEnc/360], ylim, 'k'); ylim([0 360]);xlim([0 135]); xlabel('Laps in lab frame');ylabel('Angle in Landmark Frame');
subplot(3,1,3)
hippGainEnc = interp1(specGain.encAngle,specGain.popGain,rosdata.encAngle,'linear','extrap');
        rosHippAngleMon = rosdata.encAngle(1)+ cumsum([0 diff(rosdata.encAngle)].*hippGainEnc);
         hippAngle      = interp1(rosdata.encTimes, rosHippAngleMon, clustExtended(i).ts,'linear','extrap'); 
plot(clustExtended(i).encAngle(velIdx)/360, mod(hippAngle(velIdx)+180, 360), '.'); hold on; 
 ylim([0 360]);  xlim([0 135]); xlabel('Laps in lab frame');ylabel('Angle in Hippocampal Frame');plot([landOffEnc/360 landOffEnc/360], ylim, 'k');
 %%
load('C:\Users\Ravi\Dropbox\Dome\useThisSpecGains\Rat692\160919_Rat692-04\analyzed\m1_specGain.mat')
load('F:\dome_data\Rat692\160919_Rat692-04\analyzed\m1_rosdata.mat')
load('F:\dome_data\Rat692\160919_Rat692-04\analyzed\m1_clustExtended.mat')
specGain.popGain = 1-nanmedian([specGain.clGain.gain],2);
figure(1); clf;
try
        landMsgIdx  = find(strcmp({rosdata.domeVisMsgs.type},'landmarks'));
        b           = [rosdata.domeVisMsgs.visible];
        landVis     = b(landMsgIdx);
        landOffTime = rosdata.domeVisTimes(landMsgIdx(find(landVis == 0,1)));
    catch
        idx         = find(strcmp({rosdata.domeEvMsgs.name},'landmarks'));
        landOffTime = rosdata.domeEvTimes(idx(find([rosdata.domeEvMsgs(idx).  value] == 0,1)));
end

if ~isempty(landOffTime)
    [~,landOffIdx]  = min(abs(rosdata.encTimes-landOffTime));
    landOffEnc      = rosdata.encAngle(landOffIdx);
else
    landOffEnc = inf;
end
    
k = 1; 
% listClust = [2 10]; 
listClust = [1 3 11 12]; %6 7 10 8 9  14 13
legendCell = {}; 
for iLoop = 1:length(listClust)
    i = listClust(iLoop); 
   uniqueTags = unique(clustExtended(i).fieldTag); 
   uniqueTags(isnan(uniqueTags)) = [];
   
    for fLoop = 1:length(uniqueTags)
        idx = clustExtended(i).fieldTag == uniqueTags(fLoop); 
        if sum(idx) < 10
            continue; 
        end
       fieldVec = [[clustExtended(i).fields.encAngle1]' [clustExtended(i).fields.encAngle2]'];
       
       fieldVec =fieldVec(idx,:); 
       clustGainEnc = mean(fieldVec(2:end,:),2); 
       clustGain = 360./mean(diff(fieldVec),2); 
       TF = isoutlier(clustGain,'movmedian',7); 
       outlierReplacementGain = interp1(specGain.encAngle,specGain.popGain,clustGainEnc(TF));
       ratioGain = outlierReplacementGain./clustGain(TF); 
       ratioInt = round(ratioGain);
       ratioRem = abs(ratioGain - ratioInt);
       ratioRemThresh = (ratioRem < 0.1 | ratioRem > 0.9) &  ismember(ratioInt,[2 3 4]);
       clustGain(TF) = clustGain(TF).*ratioInt; 
%        clustGain = filloutliers(clustGain,'center');  
       plot(clustGainEnc/360, (clustGain));hold on; 
       
%        fieldVec = [clustExtended(i).fields.encAnglePeak]';
%        plot(fieldVec(2:end), filloutliers(360./diff(fieldVec),'linear'));

%         clustGainEnc = mean(fieldVec(2:end,:),2); 
        
%         clustGain = 360./(diff(mean(fieldVec,2))); 
%         TF = isoutlier(clustGain); 
%         outlierReplacementGain = interp1(specGain.encAngle,specGain.popGain,clustGainEnc(TF)); 
%         clustGain(TF) = outlierReplacementGain; 
% 
%         idx = max(find(specGain.encAngle<clustGainEnc(2))); 
% 
%         combinedEnc = [specGain.encAngle(1:idx); clustGainEnc];
%         combinedGain = [specGain.popGain(1:idx);clustGain];
%         
% %         [combinedEnc, idx] = sort([specGain.encAngle; clustGainEnc]);
% %         combinedGain = [specGain.popGain;clustGain];
% %         combinedGain = movmean(combinedGain(idx), 10); 
% 
%         hippGainEnc_combined = interp1(combinedEnc,combinedGain,rosdata.encAngle,'linear','extrap');
% %         rosHippAngleMon_combined = cumtrapz(rosdata.encAngle,abs(hippGainEnc_combined)) +rosdata.encAngle(1);
%         rosHippAngleMon_combined = rosdata.encAngle(1)+ cumsum([0 diff(rosdata.encAngle)].*hippGainEnc_combined);
%         
%         hippGainEnc = interp1(specGain.encAngle,specGain.popGain,rosdata.encAngle,'linear','extrap');
%         rosHippAngleMon = cumtrapz(rosdata.encAngle,abs(hippGainEnc)) +rosdata.encAngle(1);
%         rosHippAngleMon = rosdata.encAngle(1)+ cumsum([0 diff(rosdata.encAngle)].*hippGainEnc);
% % 
% %         hippAngle_combined      = interp1(rosdata.encTimes, rosHippAngleMon_combined, clustExtended(i).ts,'linear','extrap'); 
%         hippAngle      = interp1(rosdata.encTimes, rosHippAngleMon, clustExtended(i).ts,'linear','extrap'); 
% % 
%         velIdx = clustExtended(i).vel > 5; 
% %         subplot(2,1,1);
%         plot(clustExtended(i).encAngle(velIdx)/360, mod(hippAngle(velIdx)+180, 360), '.');
%         ylim([0 360]);
%         subplot(2,1,2); 
%         plot(clustExtended(i).encAngle(velIdx)/360, mod(hippAngle_combined(velIdx)+180, 360), '.');
%         ylim([0 360]);
%         pause();
        legendCell{k} =  clustExtended(i).name; 
        k = k+1; 
   end
   
end
legendCell{k} = 'Spectral gain'; 
plot(specGain.encAngle/360,specGain.popGain, 'k--', 'LineWidth', 2); 
xlabel('Laps in lab frame'); ylabel('Estimated Gain'); 
plot([landOffEnc/360 landOffEnc/360], ylim, 'k');
legendCell{k+1} = 'Landmark Off'; 
legend(legendCell);xlim([5 landOffEnc/360]); ylim([0.8 1.1]);