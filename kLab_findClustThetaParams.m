
function kLab_findClustThetaParams(varargin)   
    if isempty(varargin)
        expFolder = uigetdir();
    else 
        expFolder = varargin{1};
    end

    fprintf('\nData folder %s\n',expFolder);    
    if ~exist(fullfile(expFolder,'analyzed'),'dir')
        mkdir(fullfile(expFolder,'analyzed'));
    end

    % Load cl-files
    addpath(['.', filesep, 'packages', filesep, 'Enhanced_rdir']);
    clFiles = rdir(fullfile(expFolder,'**','*.*'), 'regexp(name,''[^OLD]cl-maze\d+\w?\.\d+'')', true);
    clFiles = {clFiles.name};
    nClust = length(clFiles);

    if nClust == 0
        error('No cluster files found from specified maze')
    end

    [clustTs] = deal(cell(nClust,1));
    for c = 1:nClust
        C = importdata(fullfile(expFolder,clFiles{c}));
        clustTs{c} = C.data(20:18:end);
    end
    fprintf('Loaded %d cluster files\n',nClust);


    thetaFile = fullfile(expFolder,'analyzed','theta_trM.mat');
    isTheta = exist(thetaFile,'file');
    if isTheta
        load(thetaFile);
        fprintf('Loaded theta file\n')
    end
    c = 0;
    for j = 1:nClust
        nameparts = strsplit(clFiles{j},filesep);
        name = [nameparts{2},'/',nameparts{3}];
        r = regexp(name,'TT(\d+)[/,\\]cl-maze(\d+).(\d+)','tokens');
        if ~isempty(r)
            c = c+1;
            clust(c).name = name;
            clust(c).ttnum = str2double(r{1}{1});
            clust(c).mazenum = str2double(r{1}{2});
            clust(c).clustnum = str2double(r{1}{3});
            ts = clustTs{j};
            clust(c).ts = ts;
            if isTheta
                cscIdx = find([theta.chData.ttNumber] == clust(c).ttnum);
                clust(c).thetaPhase = angle(interp1(theta.chData(cscIdx).phase.time,theta.chData(cscIdx).phase.wrappedData,clust(c).ts));
                clust(c).thetaTroughAmp = interp1(theta.chData(cscIdx).troughs.time,theta.chData(cscIdx).troughs.amp,clust(c).ts);
                clust(c).unwrappedThetaPhase = interp1(theta.chData(cscIdx).phase.time,theta.chData(cscIdx).phase.unwrappeddata,clust(c).ts);
                clust(c).thetaCycleIdx = floor(clust(c).unwrappedThetaPhase/(2*pi));
                clust(c).thetaTroughsFreq = interp1(theta.chData(cscIdx).troughs.time,theta.chData(cscIdx).troughs.freq,clust(c).ts);
                clust(c).symIndex = interp1(theta.chData(cscIdx).sym.time,theta.chData(cscIdx).sym.data,clust(c).ts);
            end

            % Metadata
            clust(c).dateOfProc    = date;
            clust(c).folder    = expFolder;
        end
    end
    clTheta = clust;
    fprintf('Saving cluster theta data\n');
    save(fullfile(expFolder,'analyzed','clTheta.mat'),'clTheta');
end