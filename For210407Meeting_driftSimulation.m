%%c
% Simulation as a proof of concept to test how CA3 predictive sequences and
% MEC sensory update interact in the hyperbolic world of gain
%%
clear; clc;
%%
trackRadii = 30;
nCells = 100; 
% fieldLength_angle = pi/6+ pi/20*randn(nCells,1);
fieldLength_angle = pi/6*ones(nCells,1);
fieldLength_angle = fieldLength_angle*trackRadii; 

fieldLocs_angle = 360*rand(nCells,1);
fieldLocs_xy = trackRadii*[cosd(fieldLocs_angle) sind(fieldLocs_angle)]; 

viscircles([0 0], trackRadii, 'LineWidth', 20, 'Color', 0.5*[1 1 1]); 
viscircles(fieldLocs_xy,fieldLength_angle/2); axis equal; grid on; 



%% Generate experiment protocol
laps_Epoch1     = 360*4; 
laps_Epoch2     = 360*26; 
laps_Epoch3     = 360*52; 
totalLaps       = [laps_Epoch1 laps_Epoch2 laps_Epoch3]; 
labIncr         = 0.1; 
labAngleVec     = 0:labIncr:sum(totalLaps); 
labAngleVec = labAngleVec';
% finalGainValue  = 2*rand(1); 
finalGainValue  = 0.1; 
gainVec_ep1     = ones(sum(labAngleVec < laps_Epoch1),1); 
gainVec_ep2     = linspace(1,finalGainValue,sum(labAngleVec >= laps_Epoch1 & labAngleVec < (laps_Epoch1+laps_Epoch2)))';  
gainVec_ep3     = finalGainValue*ones(sum(labAngleVec >= sum(totalLaps(1:2)) & labAngleVec <= sum(totalLaps)),1);
gainVec         = [gainVec_ep1;gainVec_ep2;gainVec_ep3]; 
relAngleVec     = labAngleVec(1) + cumsum(gainVec.*[0;diff(labAngleVec)]); 

%% Main simulation
dt = 1/14; 
rat_actualPos      = labAngleVec(1); 
rat_predPos        = rat_actualPos; 
rat_relPos         = rat_actualPos; 
mapLoc             = rat_actualPos; 

rat_vel             = 30; 
predictionOffset    = rat_vel*dt; 
predictionLength    = 1; 
fullCorrectionThresh = 0.3; 
piGain = 1; 
gainAtMoment = 1; 
weightMat = [1 0];  
i = 1;
learningRate = 0.0005; 
extinctionRate = 0.0001; 

while rat_actualPos(i) <= labAngleVec(end)
    clc; 
    (rat_actualPos(i) - labAngleVec(end))/360
    abs(rat_predPos(i) - rat_relPos(i))
%     abs(rat_predPos(i) - rat_relPos(i))
    
   
%     [~,labVecIdx]       = min(abs(rat_actualPos(i+1) - labAngleVec));
%     rat_relPos(i+1)     = relAngleVec(labVecIdx); 
    diffPos = rad2deg(circ_dist(deg2rad(rat_relPos(i)),deg2rad(rat_predPos(i)))); 
    if diffPos < 5
        piGain(i+1) = piGain(i) +  learningRate*diffPos/(piGain(i)*rat_vel*dt); 
    else
        piGain(i+1) = piGain(i); 
    end
    if abs(rat_predPos(i) - rat_relPos(i)) <= fullCorrectionThresh
        mapLoc(i+1) = rat_relPos(i); 
        piGain(i+1) = piGain(i);
        weightMat(i+1,:) = weightMat(i,:); 
    elseif abs(rat_predPos(i) - rat_relPos(i)) > fullCorrectionThresh && abs(rat_predPos(i) - rat_relPos(i)) < predictionLength
        1
%         mapLoc(i+1) = sum(weightMat(i,:).*[rat_relPos(i) rat_predPos(i)]);  
        mapLoc(i+1) = rat_relPos(i); 
        weightMat(i+1,:) = weightMat(i,:); 
    elseif abs(rat_predPos(i) - rat_relPos(i)) > predictionLength
        2
%         mapLoc(i+1) = sum(weightMat(i,:).*[rat_relPos(i) rat_predPos(i)]);
        mapLoc(i+1) = rat_relPos(i); 
%         piGain(i+1) = piGain(i); 
        weightMat(i+1,:) = weightMat(i,:) + extinctionRate*diffPos/(piGain(i)*rat_vel*dt)*[-1 1]; 
        if weightMat(i+1,1) < 0
            weightMat(i+1,:) = [0 1];
        end
    end
    rat_actualPos(i+1)  = rat_actualPos(i)+ rat_vel*dt; 
    rat_predPos(i+1)    = mapLoc(i+1)+ piGain(i)*rat_vel*dt; 
    gainAtMoment(i+1)        = interp1(labAngleVec,gainVec,rat_actualPos(i+1));
    rat_relPos(i+1)    = rat_relPos(i)+ gainAtMoment(i+1)*rat_vel*dt; 
    i = i+1; 
%     if i == 100
%         break; 
%     end
end

figure(10);clf; 
plot(rat_relPos); hold on; plot(rat_predPos)
figure(11); clf; plot(gainAtMoment);hold on; plot(piGain);
%% Hippocampal map control by cue

%% Remapping and extinction/birth rule

%% Binding rule

%% Learning rule
