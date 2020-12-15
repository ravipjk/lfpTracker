function kLab_makeLFPFromNLX(varargin)

    if isempty(varargin)
        expFolder = uigetdir();
    else 
        expFolder = varargin{1};
    end
    
    fprintf('\nData folder %s\n',expFolder);

    if ~exist(expFolder,'dir')
        error('Experiment folder does not exist')
    end
    x = strsplit(expFolder,'Neuralynx');
    if ~exist(fullfile(x{1},'extractedLFP'),'dir')
        mkdir(fullfile(x{1},'extractedLFP'));
    end

    % For each epoch, load csc files, and compute theta params
    if ismac || isunix
        addpath(['.', filesep, 'packages', filesep, 'Nlx2Mat_relDec15/binaries']);
    else
        addpath(['.', filesep, 'packages', filesep, 'MatlabImportExport_v6.0.0']);
    end
    
    ncsFiles    = dir(fullfile(expFolder,'Neuralynx','*.ncs'));
    ncsFiles    = sort({ncsFiles.name});
    nCh         = length(ncsFiles);

    cscData = cell(nCh,1);
    
        
    if ismac || isunix
        [ts,freq,Header] = Nlx2MatCSC_v3(fullfile(expFolder,'Neuralynx',ncsFiles{1}),[1 0 1 0 0],1,1,[]);
    else
        [ts,freq,Header] = Nlx2MatCSC(fullfile(expFolder,'Neuralynx',ncsFiles{1}), [1 0 1 0 0],1, 1,[]);
    end

    if length(unique(diff(ts)))>1
        warning('Data is not uniformly sampled')
    end

    freq    = mean(freq);
    startTs = ts(1);
    stopTs  = ts(end);
    

    sampleFreq = 1250;
    d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',sampleFreq);
    
    
    for k=1:nCh
        disp(['On CSC : ' num2str(k) ' of ' num2str(nCh)])
        if ismac || isunix
            samples = Nlx2MatCSC_v3(fullfile(x{1},'Neuralynx',ncsFiles{k}),[0 0 0 0 1],0,1,[]);
        else
            samples = Nlx2MatCSC(fullfile(x{1},'Neuralynx',ncsFiles{k}),[0 0 0 0 1],0,1,[]); 
        end
       
        cscData = resample(samples(:),sampleFreq,freq);
        cscData_60Notch = filtfilt(d,cscData);
        TTName = strsplit(ncsFiles{k},'.');
        TTNum = strsplit(TTName{1},'CSC');
        TTNum = str2double(TTNum{2});
        lfp.chData(k).data = cscData_60Notch; 
        lfp.chData(k).name = TTName{1};
        lfp.chData(k).chNum = TTNum;
        
        [fitObj,y, f] = lfp_findPowerFit(cscData_60Notch);
        lfp.chData(k).powerSpec.fitObj = fitObj;
        lfp.chData(k).powerSpec.yMean = y;
        lfp.chData(k).powerSpec.f = f;
    end
    timestamps = startTs:(1e6/sampleFreq):(startTs + 1e6/sampleFreq*(size(cscData,1)-1));


    fprintf('Read %d CSC channels\n',nCh);

    lfp.timestamps      = timestamps;
    lfp.channels        = nCh;
    lfp.samplingRate    = sampleFreq;
    lfp.dateOfProc      = date;
    lfp.header          = Header;
    lfp.folder          = expFolder;
    
    % Save epochname_theta file
    save(fullfile(expFolder,'extractedLFP','subSampLFP.mat'),'lfp','-v7.3');
    fprintf('Saved file %s\n','subSampLFP.mat');
end
