 function lfp_gammaProcess(varargin)
    addpath('..');
    [expFolders,epochs] = processArgs(varargin{:});
    for k = 1:length(expFolders)
        for e = 1:length(epochs{k})
             [~,dataset]     = fileparts(expFolders{k});
%             try
                disp (['lfp_gammaProcess: Attempting to process '  dataset '-' epochs{k}{e}]);
                lfp_gammaProcessFunc(expFolders{k},epochs{k}{e},varargin);
%             catch ME
%                
%                 disp (['Failed to process '  dataset '-' epochs{k}{e} ' due to error: ' newline]);
%                 disp(ME.message)
%             end
        end
    end
end

function lfp_gammaProcessFunc(expFolder,epoch,varargin)
    load(fullfile(expFolder,'analyzed',[epoch,'_lfp.mat']));
    [~,dataset]     = fileparts(expFolder);
    load('F:\Documents\gitrepos\domeanalysis\ravi\fitObject.mat');

    %% Filter signal for gamma and gamma bands
    gammafreq       = [25 150];
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
    end
    %%
    for chLoop = 1:numChannels
        clc; 
        disp (['lfp_gammaProcess: Attempting to process '  dataset '-' epoch]);
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

            chData(chLoop).fitObj     = fitObj; 
            chData(chLoop).fitObj_skip     = 1; 
            continue; 
        end
        %%
        x_gammafilt     = bz_Filter(x,'passband',gammafreq,'filter','fir1');
        lfp_rec_gamma   = [];
        tVec            = []; 
        ns              = 2^14; % number of samples
        Fs              = 1/mean(diff(lfp_timestamps));
        numLoops        = floor(length(lfp_timestamps)/ns);
        
        for i = 1:numLoops
%             clc
%             disp([num2str(i) ' of ' num2str(numLoops)])
            tVec        = [tVec;lfp_timestamps((i-1)*ns+1:ns*i)];
            in_gamma    = x_gammafilt((i-1)*ns+1:ns*i);

            [cfs,f]     = cwt(in_gamma,Fs,'amor');

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

            lfp_rec_gamma       = [lfp_rec_gamma icwt(filt_cwt)];
        end

        %% Find peaks and troughs
        lfp_rec_flip    = -lfp_rec_gamma;
        [pks, locs]     = findpeaks(lfp_rec_gamma,tVec,'MinPeakProminence',10,'MinPeakHeight',300,'MinPeakDistance',1/100);
        [pks_flip, locs_flip] = findpeaks(lfp_rec_flip,tVec,'MinPeakProminence',10,'MinPeakHeight',300,'MinPeakDistance',1/100);

        freq_pks = 1./diff([0;locs]);
        time_pks = locs;
        amp_pks = pks;
        
        freq_troughs = 1./diff([0;locs_flip]);
        time_troughs = locs_flip;
        amp_troughs = pks_flip;
        
        %% Format into structure
        chData(chLoop).channelName = TTName{1};
        chData(chLoop).channelNumber = TTNum;
        
        chData(chLoop).recreated.time = tVec*1e6;
        chData(chLoop).recreated.data = lfp_rec_gamma; 
        
        chData(chLoop).peaks.time = time_pks*1e6;
        chData(chLoop).peaks.amp = amp_pks; 
        chData(chLoop).peaks.freq = freq_pks; 
        
        chData(chLoop).troughs.time = time_troughs*1e6;
        chData(chLoop).troughs.amp = amp_troughs; 
        chData(chLoop).troughs.freq = freq_troughs; 
        
        chData(chLoop).fitObj     = fitObj; 
        chData(chLoop).fitObj_skip     = 0; 
        
    end
        gamma.chData    = chData;
        gamma.dateOfProc    = date;
        gamma.dataset   = dataset;
        gamma.folder    = expFolder;
        gamma.rat       = varargin{1}{1};
        gamma.day       = varargin{1}{2};
        gamma.epoch     = varargin{1}{3};
    %% Save data
    fprintf('Saving gamma data\n');
    save(fullfile(expFolder,'analyzed',[epoch '_gamma_trM.mat']),'gamma','-v7.3');
end