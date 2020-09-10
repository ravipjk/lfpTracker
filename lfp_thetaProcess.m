 function lfp_thetaProcess(varargin)
    addpath('..');
    [expFolders,epochs] = processArgs(varargin{:});
    for k = 1:length(expFolders)
        for e = 1:length(epochs{k})
             [~,dataset]     = fileparts(expFolders{k});
%             try
                disp (['lfp_thetaProcess: Attempting to process '  dataset '-' epochs{k}{e}]);
                lfp_thetaProcessFunc(expFolders{k},epochs{k}{e},varargin);
%             catch ME
%                
%                 disp (['Failed to process '  dataset '-' epochs{k}{e} ' due to error: ' newline]);
%                 disp(ME.message)
%             end
        end
    end
end

function lfp_thetaProcessFunc(expFolder,epoch,varargin)
    load(fullfile(expFolder,'analyzed',[epoch,'_lfp.mat']));
   
    [~,dataset]     = fileparts(expFolder);
    load('F:\Documents\gitrepos\domeanalysis\ravi\fitObject.mat');

    %% Filter signal for theta and gamma bands
    thetafreq       = [3 30];
    numChannels     = length(lfp.chData);
    lfp_chData      = lfp.chData;
    lfp_timestamps      = lfp.timestamps; 
    for chLoop = 1:length(lfp.chData)
        chData(chLoop).channelName = '';
        chData(chLoop).channelNumber = [];

        chData(chLoop).recreated.time = [];
        chData(chLoop).recreated.data = []; 

        chData(chLoop).peaks.time = [];
        chData(chLoop).peaks.amp = []; 
        chData(chLoop).peaks.freq = []; 

        chData(chLoop).troughs.time = [];
        chData(chLoop).troughs.amp = []; 
        chData(chLoop).troughs.freq = []; 
        
        chData(chLoop).sym.time = [];
        chData(chLoop).sym.data = []; 
        
        chData(chLoop).phase.time = [];
        chData(chLoop).phase.data = [];
    end
    
    for chLoop = 1:numChannels
        clc; 
        disp (['lfp_thetaProcess: Attempting to process '  dataset '-' epoch]);
        TTName = strsplit(lfp_chData(chLoop).channel,'.');
        TTNum = strsplit(TTName{1},'CSC');
        TTNum = str2double(TTNum{2});
        disp(['On CSC : ' num2str(chLoop) ' of ' num2str(numChannels)])
        x               = lfp_chData(chLoop).lfp;  
        
        %% Find 1/f fit 
%         fitObj = lfp_findPowerFit(x);
        fitObj = lfp_chData(chLoop).fitObj; 
        if fitObj.a < 100 || fitObj.b > 0
            chData(chLoop).channelName = TTName{1};
            chData(chLoop).channelNumber = TTNum; 

            chData(chLoop).recreated.time = [];
            chData(chLoop).recreated.data = []; 

            chData(chLoop).peaks.time = [];
            chData(chLoop).peaks.amp = []; 
            chData(chLoop).peaks.freq = []; 

            chData(chLoop).troughs.time = [];
            chData(chLoop).troughs.amp = []; 
            chData(chLoop).troughs.freq = []; 

            chData(chLoop).sym.time = [];
            chData(chLoop).sym.data = []; 

            chData(chLoop).phase.time = [];
            chData(chLoop).phase.data = [];
            chData(chLoop).fitObj     = fitObj; 
            chData(chLoop).fitObj_skip     = 1; 
            continue; 
        end
        %%
        x_thetafilt     = bz_Filter(x,'passband',thetafreq,'filter','fir1');
        lfp_rec_theta   = [];
        tVec            = []; 
        ns              = 2^14; % number of samples
        Fs              = 1/mean(diff(lfp_timestamps));

        numLoops = floor(length(lfp_timestamps)/ns);
        for i = 1:numLoops
%             clc
%             disp([num2str(i) ' of ' num2str(numLoops)])
            tVec        = [tVec;lfp_timestamps((i-1)*ns+1:ns*i)];
            in_theta    = x_thetafilt((i-1)*ns+1:ns*i);

            [cfs,f]     = cwt(in_theta,Fs,'amor');

            idx         = f < 120; 
            f           = f(idx); 
            cfs         = cfs(idx,:);
            p_Fit       = feval(fitObj,f)/1.2;
            p_Fit(f < 25)       = 5;

            multipFact  = repmat(p_Fit, 1,size(abs(cfs),2));
            cwtMat_invF_comp    = abs(cfs)./multipFact;

            threshIdx   = cwtMat_invF_comp < 400;

            filt_cwt    = cfs;
            filt_cwt(threshIdx) = 0;

            lfp_rec_theta       = [lfp_rec_theta icwt(filt_cwt)];
        end

        %% Find peaks and troughs
        lfp_rec_flip = -lfp_rec_theta;
        [pks, locs] =findpeaks(lfp_rec_theta,tVec,'MinPeakProminence',10,'MinPeakHeight',300,'MinPeakDistance',1/12);
        [pks_flip, locs_flip] =findpeaks(lfp_rec_flip,tVec,'MinPeakProminence',10,'MinPeakHeight',300,'MinPeakDistance',1/12);

        freq_pks = 1./diff([0;locs]);
        time_pks = locs;
        amp_pks = pks;
        
        freq_troughs = 1./diff([0;locs_flip]);
        time_troughs = locs_flip;
        amp_troughs = pks_flip;

        %% Compute Assymetry index
        troughPeak_Mat = [];
        locsCopy = locs;
        k = 1; 
        for i = 1:length(locs_flip)
            locsCopyIdx = find((locsCopy-locs_flip(i))<0 & (locsCopy-locs_flip(i)) > -0.5, 1, 'last' );
            if isempty(locsCopy)
                break;
            elseif isempty(locsCopyIdx)
                continue;
            else
                troughPeak_Mat(k,1) = locs_flip(i); % Trough
                troughPeak_Mat(k,2) = locsCopy(locsCopyIdx); % Peak before the trough
                locsCopy(1:locsCopyIdx) = [];
                k = k+1; 
            end
        end
        symTime  = troughPeak_Mat(1:end-1,1);
        symIndex = abs(diff(troughPeak_Mat(1:end-1,:),1,2))./diff(troughPeak_Mat(:,1)); 
        %% Compute Phase
        % Defining Trough to be zero phase
        phaseVec_pk_trough = zeros(2*size(troughPeak_Mat,1),2);
        tempVar = (fliplr(troughPeak_Mat))';
        phaseVec_pk_trough(:,1) = tempVar(:);
        phaseVec_pk_trough(:,2) = 180:180:180*size(phaseVec_pk_trough,1);
        
        phaseVec = interp1(phaseVec_pk_trough(:,1), phaseVec_pk_trough(:,2), tVec);
        
        %% Format into structure
        chData(chLoop).channelName = TTName{1};
        chData(chLoop).channelNumber = TTNum; 
        
        chData(chLoop).recreated.time = tVec*1e6;
        chData(chLoop).recreated.data = lfp_rec_theta; 
        
        chData(chLoop).peaks.time = time_pks*1e6;
        chData(chLoop).peaks.amp = amp_pks; 
        chData(chLoop).peaks.freq = freq_pks; 
        
        chData(chLoop).troughs.time = time_troughs*1e6;
        chData(chLoop).troughs.amp = amp_troughs; 
        chData(chLoop).troughs.freq = freq_troughs; 
        
        chData(chLoop).sym.time = symTime*1e6;
        chData(chLoop).sym.data = symIndex; 
        
        chData(chLoop).phase.time = tVec*1e6;
        chData(chLoop).phase.data = phaseVec;
        chData(chLoop).fitObj     = fitObj; 
        chData(chLoop).fitObj_skip     = 0; 
    end
        theta.chData    = chData;
        theta.dateOfProc    = date;
        theta.dataset   = dataset;
        theta.folder    = expFolder;
        theta.rat       = varargin{1}{1};
        theta.day       = varargin{1}{2};
        theta.epoch     = varargin{1}{3};
    %% Save data
    fprintf('Saving theta data\n');
    save(fullfile(expFolder,'analyzed',[epoch '_theta_trM.mat']),'theta','-v7.3');
end