clear; clc; 
load('F:\dome_data\Rat638\160617_Rat638-15\analyzed\m1_v2_specGain.mat')
load('F:\dome_data\Rat638\160617_Rat638-15\analyzed\m1_rosdata.mat')
load('F:\dome_data\Rat638\160617_Rat638-15\analyzed\m1_hipp_hipp_clustExtended.mat')
%%
clear; clc;
addpath('F:\Documents\gitrepos\domeanalysis\');
addpackagepath('CircStat2012a');
%% Load all data
load expt.mat
% load expt_interneuron.mat

nExpt       = size(expt,2);
ratList     = [expt.rat];
dayList     = [expt.day];
epochList   = {expt.epoch};

warning ('off','MATLAB:dispatcher:nameConflict');

methodsProcessedData = struct;
nExpt = length(expt);


nExpt = length(expt);
for e = 1:nExpt
    endIdx = find(expt(e).specGain.encAngle>=(expt(e).rosdata.landOffAngle-6*360),1);
    if isempty(endIdx)
        endIdx = length(expt(e).specGain.encAngle);
    end 
    encAngle = expt(e).specGain.encAngle(1:endIdx)-expt(e).rosdata.encAngle(1);
    filtExptGain = expt(e).specGain.filtExptGain(1:endIdx)';
    filtPopGain = expt(e).specGain.filtPopGain(1:endIdx);
    
    periodicity = filtPopGain./filtExptGain;
    gainError = filtPopGain - filtExptGain;
    meanPeriodicity = nanmean(periodicity(filtExptGain>0.2));
    meanGainError = nanmean(gainError(filtExptGain>0.2));
    autoLmControl = meanPeriodicity>0.9 & meanPeriodicity<1.1;
    
    expt(e).meanPeriodicity = meanPeriodicity;
    expt(e).meanGainError = meanGainError;
    expt(e).autoLmControl = autoLmControl;
end
%%

for e = 1:nExpt
    X = []
    clc; 
    e
    disp(['Final Gain: ' num2str(expt(e).finalGain)]); 
    disp(['AutoLM: ' num2str(expt(e).autoLmControl)]); 
    
    if (expt(e).specGain.encAngle(end)-expt(e).rosdata.landOffAngle)>=6*360
        j=j+1;
        idx = find(expt(e).specGain.encAngle>=(expt(e).rosdata.landOffAngle+6*360),1);
        finalGain(e) = expt(e).finalGain; 
        recalibGain(e) = nanmedian(expt(e).specGain.filtPopGain(idx: idx));
        disp(['Recal Gain: ' num2str(recalibGain)]); 
    end
end   
%%
if 1
% for e = 1:1    
%     if expt(e).autoLmControl
        
        ratNum = expt(e).rat;
%         if ratNum ~= ratId(5)
%             continue;
%         end
        dayNum = expt(e).day;
        epoch = expt(e).epoch; 
        
        landOffAngle = expt(e).rosdata.landOffAngle;
        landOffTime = expt(e).rosdata.landOffTime;
        landOffIdx  = expt(e).rosdata.landOffIdx;
            
        [expFolders,epochs] = processArgs(ratNum, dayNum, epoch);
        
        % Load extended cluster file for epoch
        load(fullfile(expFolders{1},'analyzed',[epoch '_clustExtended.mat']));
%         gammaExtFile = fullfile(expFolders{1},'analyzed',[epoch '_gamma_trM.mat']);
        load(fullfile(expFolders{1},'analyzed',[epoch '_rosdata.mat']));
        figure(2); clf; 
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
        listClust = [2 10]; 
        legendCell = {}; 
        freqEdges = [20 75 150 350]; 
        phaseEdges = [0 180 360]; 
        for i = 1:length(clustExtended)

           uniqueFieldIDs  = unique(clustExtended(i).fieldTag);
            uniqueFieldIDs  = uniqueFieldIDs(~isnan(uniqueFieldIDs));
            clStability = [];
            for iField = 1:length(uniqueFieldIDs)  
                landOffidx      = (clustExtended(i).fieldTag == uniqueFieldIDs(iField)) & [clustExtended((i)).fields.encAngle2] <=landOffEnc;
        %         landOffidx      = (clustExtended(i).fieldTag == uniqueFieldIDs(iField));
                if sum(landOffidx) <5
                                continue;   
                end
        %         field_relRange      = [mod([clustExtended((i)).fields.relAngle1]',360) mod([clustExtended((i)).fields.relAngle2]',360)];
                field_encRange      = [[clustExtended((i)).fields.encAngle1]' [clustExtended((i)).fields.encAngle2]'];
                field_encRange      = field_encRange(landOffidx,:);
                field_relRange      = [[clustExtended((i)).fields.relAngle1]' [clustExtended((i)).fields.relAngle2]'];
                field_relRange      = field_relRange(landOffidx,:);

                medianFieldLimit = wrapTo360([rad2deg(circ_mean(deg2rad(field_relRange(:,1)))) rad2deg(circ_mean(deg2rad(field_relRange(:,2))))]);
                medianFieldMid = wrapTo360(rad2deg(circ_mean(deg2rad(medianFieldLimit(1)), deg2rad(medianFieldLimit(2)))));
                medianFieldSize = abs(rad2deg(circ_dist(deg2rad(medianFieldLimit(1)), deg2rad(medianFieldLimit(2)))));

               fieldVec_times = [[clustExtended(i).fields.ts1]' [clustExtended(i).fields.ts2]'];
               fieldVec_times = fieldVec_times(landOffidx,:); 

               for j = 1:length(fieldVec_times)
                   idx = clustExtended(i).ts >= fieldVec_times(j,1) & clustExtended(i).ts <= fieldVec_times(j,2); 
                   clustTs = clustExtended(i).ts(idx); 
                   spikeFiringFreq = [0; 1./(diff(clustTs)/1e6)]; 
                   scaledRelAngle = wrapTo360(rad2deg(circ_dist(deg2rad(mod(clustExtended((i)).relAngle(idx),360)), deg2rad(medianFieldMid)*ones(size(clustExtended((i)).relAngle(idx))))))/medianFieldSize; 
%                    scaledRelAngle = (clustExtended((i)).relAngle(idx)-field_relRange(j,1))./(field_relRange(j,2)-field_relRange(j,1));
                   selFreq = scaledRelAngle >=0 & scaledRelAngle<=1;
        %            figure(1); 
                   plot3(mean(clustExtended((i)).relAngle(idx))*ones(size(spikeFiringFreq(selFreq)))/360, scaledRelAngle(selFreq),movmean(spikeFiringFreq(selFreq),3), 'k'); hold on; grid on; 
        %            figure(2); 
                   N = histcounts(spikeFiringFreq(selFreq), freqEdges); 
%                    N = histcounts(rad2deg(clustExtended(i).thetaPhase(selFreq))+180, phaseEdges); 
%                    plot3(mean(clustExtended((i)).encAngle(idx))*ones(size(N))/360, freqEdges(1:end-1),(N), 'k'); hold on; grid on; 
                    binCen =  movmean(freqEdges,2);
%                     binCen =  movmean(phaseEdges,2);
                    X = [X; mean(clustExtended((i)).encAngle(idx))*ones(size(N))'/360, binCen(2:end)' (N./(sum(N)))'];
% %                    plot3(mean(clustExtended((i)).encAngle(idx))*ones(size(N))/360, binCen(2:end),(N./(sum(N))), 'k'); hold on; grid on; 
               end
           end
        end
%         clf; 
%         scatter3(X(:,1), X(:,2),X(:,3), 'k'); hold on; grid on; 
        pause();
%         plot3(landOffEnc*ones(size(ylim))/360, ylim,zlim, '--k', 'LineWidth', 3);
end