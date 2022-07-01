% ===============================================================================
% Matlab Code 
% ===============================================================================
% 
% Copyright (C) 2020, Mahmoud Essalat
% (mahmoudessalat[at]ucla[dot]edu)
% 
% 
% These programs and documents are distributed without any warranty,
% express or implied.  As the programs were written for research
% purposes only, they have not been tested to the degree that would be
% advisable in any important application.  All use of these programs is
% entirely at the user's own risk.
% 
% ===============================================================================
 
 
 
 
clear; %close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
  
%% Global Final Variables
path_orig = '~/Desktop/research/Robas/Data/AT-HOME STUDY - wo videos';% CHANGE 1 : add the path to at home study data

addpath(genpath('~/Desktop/research/Robas'));
%% Global Variables
SampRate = 25;
SamplePeriod = 1/SampRate;
classes = 1:16;
%% Intitialization

participant_names = {'Isabel', 'Sumukh', 'Apurva', 'Quit', 'Rachel', ...
    'Jorge', 'Dayanara', 'Mariana', 'Shivam', 'Nandini', 'Christian', ...
    'Veronica', 'June', 'Maria', 'David'};    

testMode = 0;
 
ids = [1:15]; 

Days_labeled = {[23,24,26:32,33], [19:28], [1:10], [], [1:3,5,6,8:11,15], [1:10], [15:24], [13:15, 17:22, 23], [22:31], [10,11,13:15,18:22], ...
    [7:14,16,17], [2,6:9,31,32,34,36,37], [18:27], [8,11,14,17:19,22,25:27], []};
 
% 10 Good Sessions with 0.5 sec resolution for Results Section

%P1 [23:37] % check this 33-37 added (9/15)
%P2 [18:28] % 11 sessions
%P3 [1:10]
%P4 []
%P5 [1:3,5,6,8:11,15]
%P6 [1:10]
%P7 [15:24]
%P8 [13:27] % check this 23-27 added (9/15)
%P9 [22:31]
%P10 [10,11,13:15,18:22]
%P11 [7:14,16,17] 
%P12 [2,6:9,31,32,34,36,37]
%P13 [18:27]
%P14 [8,11,14,17:19,22,25:27] 
%P15 []


%% problems
% P2 day 18, P8 day 16 basic AI should check the labels

% P1 day 23 it should be on 20210217 from the file name but in desiredTimes it shows that it is seems {'2021-01-03T20:25:46'}. It seems there is a shift in desiredTimes

% P8 day 13 labels seem to have a shift
% P5 day 5 its out of bounds, labels should be shifted
% P6 day 8 labels should be checked again , P8 day 22 labels have a little bit of shift 
% P12 Day 6 labels seem to have a shift
% P12 Day 8 labels seem to need a small shift
% P12 Day 37 there is one sample difference between YSegTruthSession and YSegPredSession
% P5 Day 8 there is two sample difference between YSegTruthSession and YSegPredSession

% P2 day 21 the accuracy is 4% ! P3 day 1 the accuracy is 8%! P3 day 2 the accuracy is 10%! P3 day 3 the accuracy is 9%! P3 day 4 the accuracy is 16%! 
% P3 day 5 the accuracy is 12%! P3 day 6 the accuracy is 13%! P3 day 8 the accuracy is 11%! P3 day 9 the accuracy is 14%! P3 day 10 the accuracy is 9%!
% P6 day 7 the accuracy is 13%! P6 day 8 the accuracy is 9%! P6 day 10 the accuracy is 18%! P9 day 23 the accuracy is 2%! P9 day 24 the accuracy is 4%!
% P9 day 25 the accuracy is 12%! P9 day 26 the accuracy is 10%!  P9 day 27 the accuracy is 11%! P9 Day 30 the accuracy is 10%! P9 Day 31 the accuracy is 10%!
% P10 idk if i missed any low accuracy!
% P11 day 7 the accuracy is 4%! P11 day 8 the accuracy is 10%! P11 day 16 the accuracy is 8%!P11 day 17 the accuracy is 10%!
% P14 day 8 the accuracy is 10%! P14 day 14 the accuracy is 5%! P14 day 17 the accuracy is 6%! P14 day 18 the accuracy is 4%!
% P14 day 19 the accuracy is 3%! P14 day 22 the accuracy is 6%! P14 day 27 the accuracy is 10%!
%% Getting data  

for f = 1:length(ids)
    
    path_part = strcat(path_orig, '/Participant', {' '}, num2str(ids(f)),  {' - '}, participant_names(ids(f)));
    path_part = path_part{1};
    
    for g = 1:length(Days_labeled{ids(f)})
        %% Data Loading
        [timeOralB, accOralB, gyrOralB, pressureOralB]  = pathToDataNewNewest(path_part, Days_labeled{ids(f)}(g), ids(f));         

        [lineInd, region_names, lineIndStart, lineIndEnd] = DaytolineIndandRegAtHome(Days_labeled{ids(f)}(g),ids(f));
     

        %%%% For at home study videos with time labels(in second) uncomment these three lines       
        brushStartTimeVid = lineIndStart(1);    
        lineIndStart = floor((lineIndStart-brushStartTimeVid)*SampRate)+1;
        lineIndEnd = floor((lineIndEnd-brushStartTimeVid)*SampRate)+1;
        %%%%
    
        lineIndStartNew = [lineIndStart, lineIndEnd(end)];
        lineIndEndNew = [lineIndStart(1), lineIndEnd];
        lineIndNew = floor(0.5*(lineIndStartNew+lineIndEndNew));
        region_names_new = [region_names, {'Rest'}];
    
        %%lineInd modification
        shiftLineInd = size(accOralB,1)-lineIndEndNew(end)
%         lineIndStartNew = lineIndStartNew + shiftLineInd;
%         lineIndEndNew = lineIndEndNew + shiftLineInd;
%         lineIndNew = lineIndNew + shiftLineInd;
%         lineIndStartNew(1) = 1;
%         lineIndEndNew(1) = 1;
%         lineIndNew(1) = 1;
        
        %%lineInd modification
        lineIndStartNew(end) = size(accOralB,1);
        lineIndEndNew(end) = size(accOralB,1);
        lineIndNew(end) = size(accOralB,1);
        
        figure;plot(accOralB);
        xlabel('sample');
        ylabel('g');
        title('Accelerometer OralB Scaled');
        legend('X', 'Y', 'Z');   
        vline_new(lineIndStartNew, 'k', region_names_new);               
%         vline_new(lineIndEndNew, 'b', '.');

        
%         figure;plot(gyrOralB);
%         xlabel('sample');
%         ylabel('g');
%         title('Gyroscope OralB Scaled');
%         legend('X', 'Y', 'Z');   
%         vline_new(lineIndStartNew, 'k', region_names_new);

%         vline_new(lineIndEndNew, 'k', '.');

    
%         figure;plot(gyrOralB);
%         xlabel('sample');
%         ylabel('rad/s');
%         title('Gyroscope OralB Scaled');
%         legend('X', 'Y', 'Z'); 
%         vline(lineIndStartNew, 'k', region_names_new);
%     
%         figure;plot(accOralB);
%         xlabel('sample');
%         ylabel('g');
%         title('Accelerometer OralB Scaled Avg');
%         legend('X', 'Y', 'Z');   
%         vline(lineIndNew, 'k', region_names_new);
  
        
       %% LPF accelerometer
               
%       accOralBStat = LowPass(accOralB, 25, 3, 1);
%       accOralB = accOralBStat;

        %% For generating train data mode sample
        [EA] = eulersMadgAccGyr(accOralB, gyrOralB, SamplePeriod);        
        
        %% generate lineInd
        startIndex = lineIndStartNew(1);      
        endIndex = lineIndEndNew(end);
        predTransLTD = TransDetect(accOralB(startIndex : endIndex, :));
        predTransLTD = predTransLTD + startIndex - 1;
        predTransLTD = [lineIndStartNew(1), predTransLTD, lineIndEndNew(end)];
        
        
        savingName = strcat('AtHome', 'P', num2str(ids(f)), 'Day', num2str(Days_labeled{ids(f)}(g)));
        EAunwrapped = labelsamples(accOralB, gyrOralB, [], [], EA, lineIndNew, region_names_new, savingName, lineIndStartNew, lineIndEndNew, predTransLTD, pressureOralB);
        
        %% For generating train data mode Chunk features
        
        %  this part uses dtw features and so load XTrain, so make sure XTrain is for all is
        %  saved for leaveOneOut = 0
        featureExtractChunk(accOralB, gyrOralB, [], [], EAunwrapped, lineIndNew, ...
            region_names_new, savingName, ids(f), Days_labeled{ids(f)}(g), lineIndStartNew, lineIndEndNew, predTransLTD)
        
        %% For generating train data mode Chunk lstm       
            
        labelChunkLSTM(accOralB, gyrOralB, [], [], EAunwrapped, lineIndNew, ...
            region_names_new, savingName, ids(f), Days_labeled{ids(f)}(g), lineIndStartNew, lineIndEndNew, predTransLTD)
 
        %% For generating train data Gyr Peaks Features
%         featureExtractGyrPeaks(gyr, lineInd, ...
%             region_names_new, ids(f), Days_labeled{ids(f)}(g), lineIndStart, lineIndEnd)
 
        
        
        YSegSessionWithZeros = zeros(size(accOralB,1),1);
        for k = 1:length(lineIndStartNew)-1
            YSegSessionWithZeros(lineIndStartNew(k):lineIndEndNew(k+1)-1)  = regionNamestoNo(region_names(k));
        end
%         YSegPredSession = [YSegPredSession;YSegPredSession(end)]; %the last sample doesn't get assigned with translines
        
        sessionBrushingScore = brushingScoreFunction(YSegSessionWithZeros, pressureOralB)
        
        AccStat = LowPass(accOralB, SampRate, 3, 1);
        accN = normr(AccStat);               
              
        features = [accN,EA];

        netName = 'AtHomeLSTMDataSegWO1';
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
        
%         figure;plot(EA)               
%         vline_new(predTransLTD, 'k', [regionNostoNames(YSegPred1), {'Rest'}]);

        YSegPredSession = [];
        for k = 1:length(XSegTest)
            YSegPredSession  = [YSegPredSession; repmat(YSegPred1(k), length(XSegTest{k}), 1)];
        end
        YSegPredSession = [YSegPredSession;YSegPredSession(end)]; %the last sample doesn't get assigned with translines
        
        YSegTruthSession = [];
        for k = 1:length(lineIndNew)-1
            YSegTruthSession  = [YSegTruthSession; repmat(regionNamestoNo(region_names_new{k}), lineIndNew(k+1)-lineIndNew(k), 1)];
        end
        YSegTruthSession = [YSegTruthSession;regionNamestoNo(region_names_new{end})]; %the last sample doesn't get assigned with translines
        
        accuracy = nnz(YSegTruthSession == YSegPredSession)/length(YSegTruthSession)
        
        sessionBrushingScorePred = brushingScoreFunction(YSegPredSession, pressureOralB)
        
    end
    close all;
end

return;
for i = 1:10
    for g = 1:10
        
        %% Plot 
    
   
        pauseTime = 0.5;
        windowTime = 10;
        slideTime = 5;

        windowLen = windowTime*SampRate;
        slideLen = slideTime*SampRate;

        brushStartSample = 1;
        brushEndSample = size(accOralB,1);
%         [brushStartSample, brushEndSample] = StartEndBrushSampAHS(Days_labeled{ids(f)}(g), ids(f));
%         brushEndSample = brushStartSample+3;
%         brushStartSample = 3;brushStartTime*SampRate;%CHANGE 3 : starting sample of his brushing session % SAVE 1
%         brushEndSample = 5;brushEndTime*SampRate;%CHANGE 4 : end sample of his brushing session % SAVE 2
    
        accClipped = accOralB(brushStartSample:brushEndSample,:);
        gyrClipped = gyrOralB(brushStartSample:brushEndSample,:);
    
%     if plotNo == 1
%         desSig = accClipped;
 
        figData = figure('NumberTitle', 'off', 'Name', ...
            strcat('Data for Participant ',...
            num2str(ids(f)), ' Session ', num2str(Days_labeled{ids(f)}(g))), ...
            'Position', [-20 0 600 800]);
        ax(1) = subplot(2,1,1);
        xlabel('sample');
        ylabel('g');
        title(strcat('Data for Participant ',...
            num2str(ids(f)), ' Session ', num2str(Days_labeled{ids(f)}(g)),...
            ' Accelerometer'));
        ylim([-1.3 1.3]);  
        txtPosAcc = -1.15;
%     elseif plotNo == 2
%         desSig = magClipped;
        
%         figMag = figure('NumberTitle', 'off', 'Name', 'Magnetometer', ...
%             'Position', [900 220 600 500]);
%         
       

        ax(2) = subplot(2,1,2);
        xlabel('sample');
        ylabel('degree/s');
        title('Gyroscope');
        ylim([-300 300]);  
        txtPosGyr = -250;
%     end
        No_Wind = max(0,floor((size(accClipped,1) - windowLen)/slideLen))+1;

        subplot(ax(1));
        hold on;
        plot(accClipped(:,1), 'r');
        plot(accClipped(:,2), 'g');
        plot(accClipped(:,3), 'b');
%     legend('X', 'Y', 'Z');
     
    
%     legend('X', 'Y', 'Z');
    
    subplot(ax(2));
    hold on;
    plot(gyrClipped(:,1), 'r');
    plot(gyrClipped(:,2), 'g');
    plot(gyrClipped(:,3), 'b');
%     legend('X', 'Y', 'Z');
    
    txtexistAcc = 0;
    txtexistGyr = 0;
    
    if size(accClipped,1) < windowLen
        for i = 1:size(accClipped,1)
            subplot(ax(1));
            plot(i, accClipped(i,1), 'r*');
            plot(i, accClipped(i,2), 'g*');
            plot(i, accClipped(i,3), 'b*');
            txt = num2str(i);
            if txtexistAcc
                delete(txthandAcc);
            end
            txthandAcc = text(i,txtPosAcc, txt);
            txtexistAcc = 1;
            set(0,'DefaultLegendAutoUpdate','off')
            drawnow;
            
            subplot(ax(2));
            plot(i, gyrClipped(i,1), 'r*');
            plot(i, gyrClipped(i,2), 'g*');
            plot(i, gyrClipped(i,3), 'b*'); 
            if txtexistGyr
                delete(txthandGyr);
            end
            txthandGyr = text(i,txtPosGyr, txt);
            txtexistGyr = 1;
            set(0,'DefaultLegendAutoUpdate','off')
            drawnow;
                     
            pause(pauseTime);  
        end
    else
        subplot(ax(1));
        xlim([1 windowLen]);
        subplot(ax(2));
        xlim([1 windowLen]);
        for i = 1:floor(windowLen*3/4)
            subplot(ax(1));
            plot(i, accClipped(i,1), 'r*');
            plot(i, accClipped(i,2), 'g*');
            plot(i, accClipped(i,3), 'b*');
            set(0,'DefaultLegendAutoUpdate','off')
            txt = num2str(i);
            if txtexistAcc
                delete(txthandAcc);
            end
            txthandAcc = text(i,txtPosAcc,txt);
            txtexistAcc = 1;
            drawnow;
            
            subplot(ax(2));
            plot(i, gyrClipped(i,1), 'r*');
            plot(i, gyrClipped(i,2), 'g*');
            plot(i, gyrClipped(i,3), 'b*');
            set(0,'DefaultLegendAutoUpdate','off')
            if txtexistGyr
                delete(txthandGyr);
            end
            txthandGyr = text(i,txtPosGyr,txt);
            txtexistGyr = 1;
            drawnow;
            
            pause(pauseTime);
            
            
        end
    end
 
curCurrPos = i;
 
subplot(ax(1));
hAcc = animatedline('Marker', '*', 'color', 'r');  
kAcc = animatedline('Marker', '*', 'color', 'g');  
lAcc = animatedline('Marker', '*', 'color', 'b'); 
 
subplot(ax(2));            
hGyr = animatedline('Marker', '*', 'color', 'r');  
kGyr = animatedline('Marker', '*', 'color', 'g');  
lGyr = animatedline('Marker', '*', 'color', 'b'); 
 
for i = 1:No_Wind-1
%     tic
    currInd = i*slideLen:i*slideLen+windowLen;
    subplot(ax(1));         
    xlim([currInd(1) currInd(end)]);
    subplot(ax(2));           
    xlim([currInd(1) currInd(end)]);
    for j = curCurrPos+1:curCurrPos+slideLen  %currInd(1):currInd(floor(mean(windowLen, slideLen)))
            subplot(ax(1));            
            addpoints(hAcc,j,accClipped(j,1));
            addpoints(kAcc,j,accClipped(j,2));
            addpoints(lAcc,j,accClipped(j,3));
            txt = num2str(j);
            if txtexistAcc
                delete(txthandAcc);
            end
            txthandAcc = text(j,txtPosAcc,txt);
            txtexistAcc = 1;
            drawnow;
            
            subplot(ax(2));            
            addpoints(hGyr,j,gyrClipped(j,1));
            addpoints(kGyr,j,gyrClipped(j,2));
            addpoints(lGyr,j,gyrClipped(j,3));
            if txtexistGyr
                delete(txthandGyr);
            end
            txthandGyr = text(j,txtPosGyr,txt);
            txtexistGyr = 1;
            drawnow;
            
            pause(pauseTime);           
    end
    curCurrPos = curCurrPos+slideLen;
%     tend = toc;
end
 
% if ~exist('currInd','var')
%     currInd = i;
% end
    
for j = curCurrPos+1:size(accClipped,1)
    subplot(ax(1));         
    xlim([curCurrPos-slideLen+1 size(accClipped,1)]);    
    subplot(ax(2));         
    xlim([curCurrPos-slideLen+1 size(accClipped,1)]);
    
    subplot(ax(1));            
    addpoints(hAcc,j,accClipped(j,1));
    addpoints(kAcc,j,accClipped(j,2));
    addpoints(lAcc,j,accClipped(j,3));
    txt = num2str(j);
    if txtexistAcc              
        delete(txthandAcc);           
    end
    txthandAcc = text(j,txtPosAcc,txt);
    txtexistAcc = 1;
    drawnow;
       
    subplot(ax(2));            
    addpoints(hGyr,j,gyrClipped(j,1));
    addpoints(kGyr,j,gyrClipped(j,2));
    addpoints(lGyr,j,gyrClipped(j,3));
    if txtexistGyr                
        delete(txthandGyr);        
    end
    txthandGyr = text(j,txtPosGyr,txt);
    txtexistGyr = 1;
    drawnow;
    
    pause(pauseTime);
    
end
 
    close all;
    end
end
   
