function [expFolders,epochs] = processArgs(varargin)
%
% In general,
% processArgs(expFolderSpecification,epochSpecification)
%
% **************************
% expFolderSpecification can be one argument or two arguments
%
% To prompt folder picking dialog
% processArgs()
% processArgs([])
%
% Directly specify full path of experiment folders
% processArgs('somepath/somewhere/Rat123/170806_Rat123-01');
% processArgs({'somepath/somewhere/Rat123/170806_Rat123-01','somepath/somewhere/Rat123/170807_Rat123-02'});
%
% Specify rats, and process all days
% processArgs(123)
% processArgs([123,256]);
%
% Specify rats, but pick days
% processArgs(123,[])
% processArgs([123,256],[]);
%
% Specify rats and days
% processArgs(123,1:10)
% processArgs([123,256],{1:5,2:7});
%
% **************************
% epochSpecification is one argument
%
% Choose all epochs from expFolderSpecification ...
% processArgs(...)
% 
% Pick epochs for each experiment folder
% processArgs(...,[])
%
% Specify epochs, else use regex to match pattern
% processArgs(...,'m1');
% processArgs(...,'m1*');
% processArgs(...,{'m1','m2*','m1b'});
%
% Manu S. Madhav
% 06-Aug-2017

%% Get data folder

addpackagepath('GetFullPath');

try
    load(fullfile('.','dataFolder.mat'),'dataFolder');
    if ~isdir(dataFolder)
        error('dataFolder is not a valid directory');
    end
catch
    q = questdlg('Select new dome data folder?','Error finding dome data folder');
    if strcmp(q,'Yes')
        dataFolder = uigetdir('.','Select Dome Data Folder');
        
        if ~isempty(dataFolder)
            save(GetFullPath(fullfile('.','dataFolder.mat'),'dataFolder'));
        else
            error('Invalid data folder');
        end        
    else
        error('Error finding dome data folder');
    end
end


%% Get experiment folders
args = varargin;
nArgs = length(args);

if nArgs==0
    expFolders = pickExpFolders();
    delArgs = 0;
else
    if isempty(args{1})
        expFolders = pickExpFolders();
        delArgs = 1;
    else
        validateattributes(args{1},{'char','cell','numeric'},{'vector'});
        if iscell(args{1})
            cellfun(@(x) validateattributes(x,{'char'},{'vector'}),args{1})
            expFolders = args{1};
            idx = cellfun(@(x) isdir(x),expFolders);
            expFolders = expFolders(idx);
            delArgs = 1;
        elseif ischar(args{1})
            expFolders = args(1);
            if ~isdir(expFolders{1})
                expFolders = [];
            end
            delArgs = 1;
        else
            rats = args{1};
            ratFolders = [];
            for k = 1:length(rats);
                ratFolder = fullfile(dataFolder,sprintf('Rat%d',rats(k)));
                d = dir(ratFolder);
                d = d(~cellfun('isempty',strfind({d.name},'_Rat')));
                [d.folder] = deal(ratFolder);
                
                ratFolders = [ratFolders;d];
            end
            
            if nArgs>1
                days = args{2};
                
                if isempty(days)
                    [idx,ok] = listdlg('PromptString','Select folders from rat','SelectionMode','multiple','ListString',{ratFolders.name});
                    if ok
                        ratFolders = ratFolders(idx);
                    end
                else
                    validateattributes(days,{'cell','numeric'},{'vector'});
                    if iscell(days)
                        if numel(days)~=length(rats)
                            error('Days not provided for all rats')
                        else
                            cellfun(@(x) validateattributes(x,{'numeric'},{'vector'}),days)
                            
                            selIdx = [];
                            for k = 1:length(days)
                                for j = 1:length(days{k})
                                    idx = find(~cellfun('isempty',regexp({ratFolders.name},sprintf('_Rat%d-%02d',rats(k),days{k}(j)))),1);
                                    selIdx = [selIdx,idx];
                                end
                            end
                            
                            ratFolders = ratFolders(selIdx);
                        end
                    else
                        if length(rats)~=1
                            error('Multiple rats provided with single day vector');
                        else
                            selIdx = [];
                            for j = 1:length(days)
                                idx = find(~cellfun('isempty',regexp({ratFolders.name},sprintf('_Rat%d-%02d',rats,days(j)))),1);
                                selIdx = [selIdx,idx];
                            end
                            
                            ratFolders = ratFolders(selIdx);
                        end
                    end
                end
                delArgs = 2;
            else
                delArgs = 1;
            end
            
            expFolders = arrayfun(@(x) fullfile(x.folder,x.name), ratFolders,'UniformOutput',false);
        end
    end
end

if delArgs
    args(1:delArgs) = [];
end

%% Get epochs
nArgs = length(args);

% Epochs can be either omitted, empty, char, cell array of chars, cell array of cell array of
% chars.

epochs = cell(size(expFolders));

if nArgs==0
    % No epochs specified - Find all the epochs in each folder!
    for k = 1:length(expFolders)
        try
            ep = findEpochs(expFolders{k});
            epochs{k} = {ep.name};
        catch
            epochs{k} = [];
        end
    end
else
    epochArg = args{1};
    
    if isempty(epochArg)
        % If empty, choose epochs for each experiment folder
        for k = 1:length(expFolders)
            ep = findEpochs(expFolders{k});
            [~,expName] = fileparts(expFolders{k});
            [idx,ok] = listdlg('PromptString',expName,'SelectionMode','multiple','ListString',{ep.name});
            if ok
                ep = ep(idx);
            end
            
            epochs{k} = {ep.name};
        end
    elseif ischar(epochArg)
        % If a single string is the epoch argument, unless an exact match is found, treat it as a regex
        % pattern to find epochs
        for k = 1:length(expFolders)
            ep = findEpochs(expFolders{k});
            idx = find(strcmp({ep.name},epochArg),1);
            if isempty(idx)
                idx = ~cellfun('isempty',regexp({ep.name},epochArg));
            end
            epochs{k} = {ep(idx).name};
        end
        
    elseif iscell(epochArg)
        % If the epoch argument is a cell array, each element of the cell
        % array is to be treated as the list of epochs
        if length(epochArg)~=length(expFolders)
            error('Length of epoch cell array must be same as number of experiment folders');
        else
            for k = 1:length(expFolders)
                ep = findEpochs(expFolders{k});
                
                if ischar(epochArg{k})
                    idx = find(strcmp({ep.name},epochArg{k}),1);
                    if isempty(idx)
                        idx = ~cellfun('isempty',regexp({ep.name},epochArg{k}));
                    end
                    epochs{k} = {ep(idx).name};
                elseif iscell(epochArg{k})
                    selIdx = [];
                    for j = 1:length(epochArg{k})
                        idx = find(strcmp({ep.name},epochArg{k}{j}),1);
                        if isempty(idx)
                            idx = find(~cellfun('isempty',regexp({ep.name},epochArg{k}{j})));
                        end
                        
                        selIdx = [selIdx,idx];
                    end
                    
                    epochs{k} = {ep(selIdx).name};
                end
            end
        end
    end
end

%% Final step: Eliminate folders for which there are no relevant epochs found
% delIdx = cellfun('isempty',epochs);
% expFolders(delIdx) = [];
% epochs(delIdx) = [];

