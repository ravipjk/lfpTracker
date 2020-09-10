clear; clc;


exptFolder      =   uigetdir('F:\dome_data\Rat692\161011_Rat692-16\Neuralynx','Select Neuralynx Data Folder');

lfp_makeLFP(exptFolder);
% lfp_thetaProcess(ratNum,dayNum,epochName);
% lfp_gammaProcess(ratNum,dayNum,epochName);
