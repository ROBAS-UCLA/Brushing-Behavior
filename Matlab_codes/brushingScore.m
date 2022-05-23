clear; %close all; clc;

%% Global Final Variables
path_orig = '~/Desktop/research/Robas/Data/RobasPhase2';% CHANGE 1 : add the path to at home study data
 
addpath(genpath('~/Desktop/research/Robas'));

%% Global Variables
SampRate = 25;
SamplePeriod = 1/SampRate;
 
%% Intitialization

participant_names = {};

for i = 1:12
    participant_names = [participant_names; strcat('robas_', num2str(i))];
end


ids = [53:60]; %run 47 later there  was an error in 20180929 %also in robas52 20181112
        
classes = 1:16; %unique(YSegTrain);

% Days_labeled = {[], [], [], [], [], [], [], [], [], [], [], [1]};

% robas17 Day 20190716, robas18 Day 20180822, robas24 Day 20180916, robas35 Day 20180822, robas41 Day 20180916,
% robas52 Day 20181020, robas52 Day 20181112 pressure will have an error that suggests that the
% function is not well written and has bugs

for f = 1:length(ids)
    
    path_part = strcat(path_orig, '/robas_', num2str(ids(f)));
%     path_part = path_part{1};
    
    addpath(genpath(path_part));

    path_sessions = dir(path_part);
    
    for g = 1:numel(path_sessions)
        name_session_curr = path_sessions(g).name;
        
        if length(name_session_curr) <= 5 || contains(name_session_curr(end-4:end), '.')
            continue
        end
        if strcmp(name_session_curr(1:2), '20') && strcmp(name_session_curr(end-2:end), 'zip')               
            unzip(name_session_curr, path_part);        
            data_dir_add = name_session_curr(1:end-4);         
        elseif strcmp(name_session_curr(1:2), '20') 
            data_dir_add = name_session_curr;
        end    

        path_session_curr = strcat(path_part, '/', data_dir_add);
    %     for g = 1:length(Days_labeled{ids(f)})
        %% Data Loading
%         [accOralB, gyrOralB, timeOralB] = pathToDataNewNewestPhase2(path_part, Days_labeled{ids(f)}(g), ids(f)); 
        [noAcc, timeOralB, accOralB, pressureOralB] = pathToDataNewNewestPhase2(path_session_curr);          
         
    end
end

for f = 1
    for g= 1
        if noAcc == 1 
            continue
        end
        
        AccStat = LowPass(accOralB, SampRate, 3, 1);
        accN = normr(AccStat);               
        
        features = accN;
                
        predTransLTD = TransDetect(accN);
        predTransLTD = [1, predTransLTD, size(accN,1)];
        
%         brushStartSample = 1;
%         predTransLTD = TransDetect(accN(startIndex : endIndex, :), mag(startIndex : endIndex, :));
%         predTransLTD = predTransLTD + startIndex - 1;
%         predTransLTD = [lineInd(1), predTransLTD, lineInd(end)];

%         figure; plot(accN)
%         title(strcat('Acc', 'P', num2str(ids(f)), 'Day', name_session_curr));
%         vline_new(predTransLTD,'k', '.');
        
        netName = 'LSTMDataSegWO1';
        load(netName);

        %% make segments
        XSegTest = {};
        for k = 2:length(predTransLTD)
            XSegTest = [XSegTest; {features(predTransLTD(k-1):predTransLTD(k)-1, :)'}];         
        end
        
        %% Split segments
        [XSegTestSplit, segNoTest] = splitSegments(XSegTest, 'Test');
        [YSegPred1Split, scoresSplit] = classify(eval([netName]),XSegTestSplit);
        

        confidenceScores = zeros(length(XSegTest), size(scoresSplit,2));
        %     YSegPred1 = zeros(length(XSegTest),1);
        for i = 1:length(XSegTest)
        %         YSegPred1(i) = majorityvote(YSegPred1Split(segNoTest==i)); 
            confidenceScores(i,:) = mean(scoresSplit(segNoTest==i,:),1);
        end
        
        

        [~,idxScoresSorted] = sort(confidenceScores, 2, 'descend');
        YSegPred1 = classes(idxScoresSorted(:,1));
%         varargout{1} = categorical(YSegPred1);
%         varargout{2} = scores;
        
        YSegPredSession = [];
        for k = 1:length(XSegTest)
            YSegPredSession  = [YSegPredSession; repmat(YSegPred1(k), length(XSegTest{k}), 1)];
        end
        YSegPredSession = [YSegPredSession;YSegPredSession(end)]; %the last sample doesn't get assigned with translines
        
        sessionBrushingScore = brushingScoreFunction(YSegPredSession, pressureOralB)
        
        fid = fopen(strcat(path_session_curr, 'Score.txt'),'wt');
        fprintf(fid, strcat(num2str(sessionBrushingScore), '\n'));
        fclose(fid);
        
    end
end

