% ===============================================================================
% Matlab Code 
% ===============================================================================
% 
% Copyright (C) 2022, Mahmoud Essalat
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
rng('default')

addpath('~/Desktop/research/Robas/Utils');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
%% Global Final Variables
path_orig = '~/Desktop/research/Robas/Data/BrushingDataSamples/';% CHANGE 1 : add the path to at home study data
path_save = '~/Desktop/research/Robas/images for dental paper/Paper Images/';% CHANGE 1 : add the path to at home study data

%% Intitialization
daysRange = [123,29,28,23,53,9,63,45,47,31,54,67,121];
mean(daysRange)
std(daysRange)

participant_names = {'Isabel', 'Sumukh', 'Apurva', 'Quit', 'Rachel', ...
    'Jorge', 'Dayanara', 'Mariana', 'Shivam', 'Nandini', 'Christian', ...
    'Veronica', 'June', 'Maria', 'David'};    

region_names_new = {'ManRO', 'ManRB', 'ManRL', 'ManAB', 'ManAL', ...
    'ManLO', 'ManLB', 'ManLL', 'MaxRO', 'MaxRB', 'MaxRL', 'MaxAB', ...
    'MaxAL', 'MaxLO', 'MaxLB', 'MaxLL'};

sampRate = 25;

ids = [1:3,5:13]; 

Days_labeled = {[23,24,26:32, 33], [19:28], [1:10], [], [1:3,5,6,8:11,15], [1:10], [15:24], [13:15, 17:22, 23], [22:31], [10,11,13:15,18:22], ...
    [7:14,16,17], [2,6:9,31,32,34,36,37], [18:27], [8,11,14,17:19,22,25:27], []};
 

%% Initialization
brushingSampleCountsAll = cell(length(ids),1);
pressureSampleCountsAll = cell(length(ids),1);

sessionDurationAll = cell(length(ids),1);
sessionDurationNonEffAll = cell(length(ids),1);
sessionPressDurationAll = cell(length(ids),1);

Allparticipants = [];

tableRegionBrushingEachReg = table;
tableRegionBrushingAllReg = table;

tableRegionPressureEachReg = table;
tableRegionPressureAllReg = table;

%% Getting data  

for f = 1:length(ids)
    if isempty(Days_labeled{ids(f)})
        continue;
    end
    for g = 1:length(Days_labeled{ids(f)})
        %% Data Loading
          
        path_data = strcat(path_orig, 'AtHome', 'P', num2str(ids(f)),  'Day', num2str(Days_labeled{ids(f)}(g)), '.csv');
      
        data = importfileAtHome(path_data);
        
        %% Session duration
        Allparticipants = [Allparticipants; f];
        
        %% Region duration
        labels = table2array(data(table2array(data(:,2))==1,1));
        brushingSampleCounts = zeros(length(region_names_new),1);
        for i = 1:length(region_names_new)
            brushingSampleCounts(i) = nnz(labels == region_names_new(i));
            tableRegionBrushingEachReg = [tableRegionBrushingEachReg; table(f, g, brushingSampleCounts(i), {region_names_new{i}(1:3)}, ...
            {region_names_new{i}(4)}, {region_names_new{i}(5)}, 'VariableNames', {'Participant', 'Session', 'sampleCounts', 'MaxOrNot', 'Side', 'OLB'})];
        %'VariableTypes', {'categorical', 'categorical', 'double', 'categorical', 'categorical', 'categorical'}
        end
        brushingSampleCountsAll{f} = [brushingSampleCountsAll{f}, brushingSampleCounts];
        
        indicesOfSwitchInRegions = find(diff(table2array(data(:,2)))==-1);
        tablePrevRegs = table2cell(data(indicesOfSwitchInRegions,1));
        cellPrevRegs = cellstr(cat(1, tablePrevRegs{:}));
        indicesOfSwitchOutRegions = find(diff(table2array(data(:,2)))==1);
        if length(indicesOfSwitchOutRegions) < length(indicesOfSwitchInRegions)
            indicesOfSwitchOutRegions = [indicesOfSwitchOutRegions; size(data,1)];
        end        
        tableNextRegs = table2cell(data(indicesOfSwitchOutRegions-1,1));
        cellNextRegs = cellstr(cat(1, tableNextRegs{:}));
            
        indicesOfSwitchSextants = indicesOfSwitchInRegions(find(regionNamestoNoQuad(cellPrevRegs) ~= regionNamestoNoQuad(cellNextRegs)));
        array  = [f, g, length(indicesOfSwitchSextants), brushingSampleCounts'];
        tableRegionBrushingAllReg = [tableRegionBrushingAllReg; array2table(array, ...
             'VariableNames', [{'Participant'}, {'Session'}, {'Number of brush transitions'}, region_names_new(:)'])]; %'VariableTypes', VarTypes
        
              
        sessionDurationAll{f} = [sessionDurationAll{f}; sum(brushingSampleCounts)]; 
                
        sessionDurationNonEffAll{f} = [sessionDurationNonEffAll{f}; size(data,1)];
         
        %% Pressure duration
        pressure = table2array(data(table2array(data(:,2))==1,3));        
        pressureSampleCounts = zeros(length(region_names_new),1);
        for i = 1:length(region_names_new)
            pressureSampleCounts(i) = nnz(labels == region_names_new(i) & pressure == 1);
            tableRegionPressureEachReg = [tableRegionPressureEachReg; table(f, g, pressureSampleCounts(i), {region_names_new{i}(1:3)}, ...
            {region_names_new{i}(4)}, {region_names_new{i}(5)}, 'VariableNames', {'Participant', 'Session', 'sampleCounts', 'MaxOrNot', 'Side', 'OLB'})];
        end
        pressureSampleCountsAll{f} = [pressureSampleCountsAll{f}, pressureSampleCounts];
        array  = [f, g, pressureSampleCounts'];
        tableRegionPressureAllReg = [tableRegionPressureAllReg; array2table(array, ...
             'VariableNames', [{'Participant'}, {'Session'}, region_names_new(:)'])]; %'VariableTypes', VarTypes
        
        
        sessionPressDurationAll{f} = [sessionPressDurationAll{f}; sum(pressureSampleCounts)];
  
    end    
end

%% Add Sextant columns
tableRegionBrushingAllReg.Sextant1 = sum(tableRegionBrushingAllReg{:,4:6},2);
tableRegionBrushingAllReg.Sextant2 = sum(tableRegionBrushingAllReg{:,7:8},2);
tableRegionBrushingAllReg.Sextant3 = sum(tableRegionBrushingAllReg{:,9:11},2);
tableRegionBrushingAllReg.Sextant4 = sum(tableRegionBrushingAllReg{:,12:14},2);
tableRegionBrushingAllReg.Sextant5 = sum(tableRegionBrushingAllReg{:,15:16},2);
tableRegionBrushingAllReg.Sextant6 = sum(tableRegionBrushingAllReg{:,17:19},2);

%% Add Sextant columns
tableRegionPressureAllReg.Sextant1 = sum(tableRegionPressureAllReg{:,3:5},2);
tableRegionPressureAllReg.Sextant2 = sum(tableRegionPressureAllReg{:,6:7},2);
tableRegionPressureAllReg.Sextant3 = sum(tableRegionPressureAllReg{:,8:10},2);
tableRegionPressureAllReg.Sextant4 = sum(tableRegionPressureAllReg{:,11:13},2);
tableRegionPressureAllReg.Sextant5 = sum(tableRegionPressureAllReg{:,14:15},2);
tableRegionPressureAllReg.Sextant6 = sum(tableRegionPressureAllReg{:,16:18},2);

%% Brushing Scores 16 regions
BrushingRegsArrays = table2array(tableRegionBrushingAllReg(:,4:end-6));
PressureRegsArrays = table2array(tableRegionPressureAllReg(:,3:end-6));

TV = sum(abs(BrushingRegsArrays./sum(BrushingRegsArrays,2) - 1/16),2);

BrushWOPress = BrushingRegsArrays - PressureRegsArrays;

alpha = 5;
brushingScore = max(100*(sum(min(BrushWOPress./sampRate, 120/16),2))/120 - alpha*TV,0);
%%% Figure 11
figure;boxplot(brushingScore, Allparticipants);
ylim([0, 100])
xlabel('Participant number');
ylabel('Score');
title('Brushing score (out of 100)');
%%%%

mean(brushingScore)
std(brushingScore)


%% Brushing Scores sextants
BrushingRegsArraysSext = table2array(tableRegionBrushingAllReg(:,end-5:end));
PressureRegsArraysSext = table2array(tableRegionPressureAllReg(:,end-5:end));

TVSext = sum(abs(BrushingRegsArraysSext./sum(BrushingRegsArraysSext,2) - 1/6),2);

BrushWOPressSext = BrushingRegsArraysSext - PressureRegsArraysSext;

alpha = 5;
brushingScoreSext = max(100*(sum(min(BrushWOPressSext./sampRate, 120/6),2))/120 - alpha*TV,0);
%%% Figure 11
figure;boxplot(brushingScoreSext, Allparticipants);
ylim([0, 100])
xlabel('Participant number');
ylabel('Score');
title('Brushing score (out of 100)');

mean(brushingScoreSext)
std(brushingScoreSext)

figure;boxplot(100*(sum(min(BrushWOPressSext./sampRate, 120/6),2))/120, Allparticipants);
ylim([0, 100])
xlabel('Participant number');
ylabel('Score');
title('Pressure compensated coverage score (out of 100)');


figure;boxplot(100*(sum(min(BrushingRegsArraysSext./sampRate, 120/6),2))/120, Allparticipants);
ylim([0, 100])
xlabel('Participant number');
ylabel('Score');
title('Coverage score (out of 100)');


figure;boxplot(alpha*TVSext, Allparticipants);
ylim([0, 5])
xlabel('Participant number');
ylabel('Score');
title('Deviation from uniform brushing score');


figure;boxplot(sum(PressureRegsArrays,2)/25, Allparticipants);
% ylim([0, 8])
xlabel('Participant number');
ylabel('Duration in seconds');
title('Total Pressure duration in brushing sessions');


figure;boxplot(repmat([1:6], size(BrushingRegsArraysSext,1),1), BrushingRegsArraysSext/25);
ylim([0, 100])
xlabel('Sextant');
ylabel('Duration in seconds');
title('Duration of brushing of different sextants');


%% Save for Doug
writetable(tableRegionBrushingAllReg, 'Duration of brushing of each dental region.csv')
writetable(tableRegionPressureAllReg, 'Duration of brushing with excessive pressure on each dental region.csv')

[part1, sess1] = meshgrid([1:length(ids)], [1:length(Days_labeled{ids(f)})]);
allPart = part1(:); allSess = sess1(:);
tableNonEffDur = table(allPart, allSess, cell2mat(sessionDurationNonEffAll), 'VariableNames', {'Participant', 'Session', 'NonEff Brushing Duration'});
writetable(tableNonEffDur, 'Non effective brushing duration.csv');


%%
AllsessionDuration = cat(1,sessionDurationAll{:})/sampRate;
AllsessionDurationNonEff = cat(1,sessionDurationNonEffAll{:})/sampRate;
AllparticipantsCat = categorical(Allparticipants);

%%% Figure A.1
f1 = figure;histogram(tableRegionBrushingAllReg.MaxRO, 20, 'Normalization', 'pdf')
xlabel('Duration (# of samples)'); ylabel('Probability');
title('Distribution of brushing duration spent on MaxRO (# of samples)');
exportgraphics(f1,strcat(path_save, 'barchart.png'),'Resolution',300)
%%%%%

figure;histogram(floor(AllsessionDuration*sampRate), 20, 'Normalization', 'pdf')
xlabel('Duration (# of samples)'); ylabel('Probability');
title('Distribution of brushing time (# of samples)');


%% Prameeter of interest : Session Duration


%% population level


%% Subject-level and session-level
%%%% Figure 4
figure;boxplot(AllsessionDuration, Allparticipants);
ylim([0,130])
% ylim([0,170])
xlabel('Participant number');
ylabel('Duration (in seconds)');
% title('Effective durations of brushing sessions (in seconds)');
%%%%%%%%


%% Modeling
%%%% Appendix section 3
sessionDurationSampTable = table(floor(AllsessionDuration*sampRate), AllparticipantsCat, 'VariableNames', {'Duration', 'Participants'});
glmeSessDurSamp = fitglme(sessionDurationSampTable,...
		'Duration ~ 1 + (1|Participants)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaSessDurSamp,betanamesSessDurSamp,statsbetaSessDurSamp] = fixedEffects(glmeSessDurSamp);
[BSessDurSamp,BNamesSessDurSamp, statsSessDurSamp] = randomEffects(glmeSessDurSamp);
meanBrushingDuration  = exp(betaSessDurSamp)/sampRate
meanStandardDeviationDuration = sqrt(exp(betaSessDurSamp)/sampRate)

interPartVar= std(exp(betaSessDurSamp+BSessDurSamp))/sampRate

stdParts = zeros(length(ids),1);
for i = 1:length(ids)
    varParts(i) = var(sessionDurationSampTable.Duration(sessionDurationSampTable.Participants == categorical(i)));
        
    stdParts(i) = std(sessionDurationSampTable.Duration(sessionDurationSampTable.Participants == categorical(i)));
end
intraPartVar = sqrt(mean(varParts))/25
intraPartVar2 = (mean(stdParts))/25

% mean(exp(normrnd(0,sqrt(cell2mat(covarianceParameters(glmeSessDurSamp))),[1,1000])))
%%%%%%

%%%% Appendix section 3
sessionDurationNonEffSampTable = table(floor(AllsessionDurationNonEff*sampRate), AllparticipantsCat, 'VariableNames', {'Duration', 'Participants'});
glmeSessDurNonEffSamp = fitglme(sessionDurationNonEffSampTable,...
		'Duration ~ 1 + (1|Participants)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaSessDurNonEffSamp,betanamesSessDurNonEffSamp,statsbetaSessDurNonEffSamp] = fixedEffects(glmeSessDurNonEffSamp);
[BSessDurNonEffSamp,BNamesSessDurNonEffSamp, statsSessDurNonEffSamp] = randomEffects(glmeSessDurNonEffSamp);
meanBrushingDurationNonEff  = exp(betaSessDurNonEffSamp)/sampRate
meanStandardDeviationDurationNonEff = sqrt(exp(betaSessDurNonEffSamp)/sampRate)

interPartVarNonEff = std(exp(betaSessDurNonEffSamp+BSessDurNonEffSamp))/sampRate

stdPartsNonEff = zeros(length(ids),1);
for i = 1:length(ids)
    stdPartsNonEff(i) = var(sessionDurationNonEffSampTable.Duration(sessionDurationNonEffSampTable.Participants == categorical(i)));
end
intraPartNonEffVar = sqrt(mean(stdPartsNonEff))/25


%% Parameter of interest: Duration of brushing with excessive pressure in each session


%% population level


%% Subject-level and session-level
%%%% Figure 9
AllsessionPressDuration = cat(1,sessionPressDurationAll{:})/sampRate;
figure;boxplot(AllsessionPressDuration, Allparticipants);
ylim([0,7])
xlabel('Participant number');
ylabel('Duration (in seconds)');
% title('Durations of brushing with excessive pressure in brushing sessions (in seconds)');
%%%%%%%%

figure;histogram(tableRegionPressureAllReg.MaxRO, 20, 'Normalization', 'pdf')
xlabel('Duration (# of samples)'); ylabel('Probability');
title('Distribution of brushing duration with extra pressure spent on MaxRO (# of samples)');

figure;histogram(floor(AllsessionPressDuration*sampRate), 20, 'Normalization', 'pdf')
xlabel('Duration (# of samples)'); ylabel('Probability');
title('Distribution of brushing time with excessive pressure (# of samples)');


%% Session level


%% Modeling
P1OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='1')==2;
P2OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='2')==2;
P3OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='3')==2;
P4OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='4')==2;
P5OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='5')==2;
P6OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='6')==2;
P7OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='7')==2;
P8OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='8')==2;
P9OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='9')==2;
P10OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='10')==2;
P11OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='11')==2;
P12OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant)=='12')==2;

resPressDummySub = EMzeropoisson_mat([tableRegionPressureEachReg.sampleCounts P1OrNotValuesPress ...
    P2OrNotValuesPress P3OrNotValuesPress P4OrNotValuesPress P5OrNotValuesPress P6OrNotValuesPress ...
    P7OrNotValuesPress P8OrNotValuesPress P9OrNotValuesPress P10OrNotValuesPress P11OrNotValuesPress P12OrNotValuesPress], ...
    0.5,'initialtau','stable','tol',0.001)
fittedParamsPressDummySub = resPressDummySub(end,:);
meanBrushingPressDurationZI = exp(fittedParamsPressDummySub(2))/sampRate
meanStandardDeviationPressDurationZI = sqrt(exp(fittedParamsPressDummySub(2))/sampRate)

%%%% Appendix section 3
sessionPressDurationSampTable = table(floor(AllsessionPressDuration*sampRate), AllparticipantsCat, 'VariableNames', {'Duration', 'Participants'});
glmeSessPressDurSamp = fitglme(sessionPressDurationSampTable,...
		'Duration ~ 1 + (1|Participants)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaSessPressDurSamp,betanamesSessPressDurSamp] = fixedEffects(glmeSessPressDurSamp);
[BSessPressDurSamp,BNamesSessPressDurSamp, statsSessPressDurSamp] = randomEffects(glmeSessPressDurSamp);
meanBrushingPressDuration = exp(betaSessPressDurSamp)/sampRate
meanStandardDeviationPressDuration = sqrt(exp(betaSessPressDurSamp)/sampRate)
%%%%%%


%% Prameter of interest: Duration of brushing on each region

%% Population level
AllsessionsRegionDuration = [];
Allregions = categorical([]);

AllsessionsRegionDurationSum = [];
AllregionsSum = categorical([]);
%% Participant-level and Session-level
%%% Figure 6
for i = 1:length(ids)
    AllsessionsRegionDurationCurrPart = reshape(brushingSampleCountsAll{i}, ...
        size(brushingSampleCountsAll{i},1)*size(brushingSampleCountsAll{i},2),1);
    AllregionsCurrPart = repmat(categorical(region_names_new'),size(brushingSampleCountsAll{i},2),1);
    
%     figure;boxplot(AllsessionsRegionDurationCurrPart/sampRate, AllregionsCurrPart);
%     xlabel('Dental surfaces');
%     ylabel('Duration (in seconds)');
% %     ylim([-2,30])    
% %     title("Time spent by " + "P" + num2str(i) + " brushing each dental surfaces (in seconds)");
%     
%     %%% Maxillary vs. Mandibular
%     AllregionsCurrPartUD = repmat({'Mandibular'}, length(AllregionsCurrPart),1);
%     idxCurrPartUp = find(AllregionsCurrPart == 'MaxRO' | AllregionsCurrPart == 'MaxRB' | AllregionsCurrPart == 'MaxRL' | AllregionsCurrPart == 'MaxAB' | ...
%     AllregionsCurrPart == 'MaxAL' | AllregionsCurrPart == 'MaxLO' | AllregionsCurrPart == 'MaxLB' | AllregionsCurrPart == 'MaxLL');
%     AllregionsCurrPartUD(idxCurrPartUp) = repmat({'Maxillary'}, length(idxCurrPartUp),1);
%     figure;boxplot(AllsessionsRegionDurationCurrPart/sampRate, AllregionsCurrPartUD);
%     xlabel('Dental surfaces');
%     ylabel('Duration (in seconds)');
% %     ylim([-2,30])
% %     title("Time spent by " + "P" + num2str(i) + " brushing maxillary and mandibular dental surfaces (in seconds)");
%    
%     %%% Right vs Left vs Anterior
%     AllregionsCurrPartSide = repmat({'Right'}, length(AllregionsCurrPart),1);
%     idxCurrPartLeft = find(AllregionsCurrPart == 'ManLO' | AllregionsCurrPart == 'ManLB' | AllregionsCurrPart == 'ManLL' | AllregionsCurrPart == 'MaxLO' | ...
%     AllregionsCurrPart == 'MaxLB' | AllregionsCurrPart == 'MaxLL');
%     AllregionsCurrPartSide(idxCurrPartLeft) = repmat({'Left'}, length(idxCurrPartLeft),1);
%     idxCurrPartAnt = find(AllregionsCurrPart == 'ManAB' | AllregionsCurrPart == 'ManAL' | AllregionsCurrPart == 'MaxAB' | AllregionsCurrPart == 'MaxAL');
%     AllregionsCurrPartSide(idxCurrPartAnt) = repmat({'Anterior'}, length(idxCurrPartAnt),1);
%     figure;boxplot(AllsessionsRegionDurationCurrPart/sampRate, AllregionsCurrPartSide);
%     xlabel('Dental surfaces');
%     ylabel('Duration (in seconds)');
% %     ylim([-2,30])
% %     title("Time spent by " + "P" + num2str(i) + " brushing right, left, and anterior dental surfaces (in seconds)");
% 
%     %%% Occlusal vs Lingual vs Buccal
%     AllregionsOLB = repmat({'Occlusal'}, length(AllregionsCurrPart),1);
%     idxCurrPartLing = find(AllregionsCurrPart == 'ManRL' | AllregionsCurrPart == 'ManAL' | AllregionsCurrPart == 'ManLL' | AllregionsCurrPart == 'MaxRL' | ...
%     AllregionsCurrPart == 'MaxAL' | AllregionsCurrPart == 'MaxLL');
%     AllregionsOLB(idxCurrPartLing) = repmat({'Lingual'}, length(idxCurrPartLing),1);
%     idxCurrPartBucc = find(AllregionsCurrPart == 'ManRB' | AllregionsCurrPart == 'ManAB' | AllregionsCurrPart == 'ManLB' | AllregionsCurrPart == 'MaxRB' | ...
%     AllregionsCurrPart == 'MaxAB' | AllregionsCurrPart == 'MaxLB');
%     AllregionsOLB(idxCurrPartBucc) = repmat({'Buccal'}, length(idxCurrPartBucc),1);
%     figure;boxplot(AllsessionsRegionDurationCurrPart/sampRate, AllregionsOLB);
%     xlabel('Dental surfaces');
%     ylabel('Duration (in seconds)');
% %     ylim([-2,30])
% %     title("Time spent by " + "P" + num2str(i) + " brushing occlusal, lingual, and buccal dental surfaces (in seconds)");
% 

    AllsessionsRegionDuration = [AllsessionsRegionDuration; AllsessionsRegionDurationCurrPart];
    Allregions = [Allregions; AllregionsCurrPart];
    
    

    %%% Maxillary vs. Mandibular


    brushingSampCurrPart = brushingSampleCountsAll{i};
    
    ManRegs = {'ManRO', 'ManRB', 'ManRL', 'ManAB', 'ManAL', 'ManLO', 'ManLB', 'ManLL'};
    idxMan = regionNamestoNo(ManRegs);
    brushingSampMan = sum(brushingSampCurrPart(idxMan,:),1)';
    MaxRegs = {'MaxRO', 'MaxRB', 'MaxRL', 'MaxAB', 'MaxAL', 'MaxLO', 'MaxLB', 'MaxLL'};
    idxMax = regionNamestoNo(MaxRegs);
    brushingSampMax = sum(brushingSampCurrPart(idxMax,:),1)';
    
    AllsessionsRegionDurationManMaxCurrPart = [brushingSampMan; brushingSampMax];
    AllregionsManMaxCurrPart = [repmat({'Mandibular'}, length(brushingSampMan),1); repmat({'Maxillary'}, length(brushingSampMax),1)];   
    figure;boxplot(AllsessionsRegionDurationManMaxCurrPart/sampRate, AllregionsManMaxCurrPart);
    xlabel('Dental surfaces');
    ylabel('Duration (in seconds)');
    
    
       
    %%% Right vs Left vs Anterior

    RightRegs = {'ManRO', 'ManRB', 'ManRL', 'MaxRO', 'MaxRB', 'MaxRL'};
    idxRight = regionNamestoNo(RightRegs);
    brushingSampRight = sum(brushingSampCurrPart(idxRight,:),1)';
    AntRegs = {'ManAB', 'ManAL', 'MaxAB', 'MaxAL'};
    idxAnt = regionNamestoNo(AntRegs);
    brushingSampAnt = sum(brushingSampCurrPart(idxAnt,:),1)';
    LeftRegs = {'ManLO', 'ManLB', 'ManLL', 'MaxLO', 'MaxLB', 'MaxLL'};
    idxLeft = regionNamestoNo(LeftRegs);
    brushingSampLeft = sum(brushingSampCurrPart(idxLeft,:),1)';
    
    AllsessionsRegionDurationSideCurrPart = [brushingSampRight; brushingSampAnt; brushingSampLeft];
    AllregionsSideCurrPart = [repmat({'Right'}, length(brushingSampRight),1); ...
        repmat({'Anterior'}, length(brushingSampAnt),1); repmat({'Left'}, length(brushingSampLeft),1)];   
    figure;boxplot(AllsessionsRegionDurationSideCurrPart/sampRate, AllregionsSideCurrPart);
    xlabel('Dental surfaces');
    ylabel('Duration (in seconds)');
       
    %%% Occlusal vs Lingual vs Buccal

    OcclusalRegs = {'ManRO', 'ManLO', 'MaxRO', 'MaxLO'};
    idxOcc = regionNamestoNo(OcclusalRegs);
    brushingSampOcc = sum(brushingSampCurrPart(idxOcc,:),1)';
    LingualRegs = {'ManRL', 'ManAL', 'ManLL', 'MaxRL', 'MaxAL', 'MaxLL'};
    idxLin = regionNamestoNo(LingualRegs);
    brushingSampLin = sum(brushingSampCurrPart(idxLin,:),1)';
    BuccalRegs = {'ManRB', 'ManAB', 'ManLB', 'MaxRB', 'MaxAB', 'MaxLB'};
    idxBucc = regionNamestoNo(BuccalRegs);
    brushingSampBucc = sum(brushingSampCurrPart(idxBucc,:),1)';
    
    
    AllsessionsRegionDurationOLBCurrPart = [brushingSampOcc; brushingSampLin; brushingSampBucc];
    AllregionsOLBCurrPart = [repmat({'Occlusal'}, length(brushingSampOcc),1); ...
        repmat({'Lingual'}, length(brushingSampLin),1); repmat({'Buccal'}, length(brushingSampBucc),1)];   
    figure;boxplot(AllsessionsRegionDurationOLBCurrPart/sampRate, AllregionsOLBCurrPart);
    xlabel('Dental surfaces');
    ylabel('Duration (in seconds)');
    
    
    AllsessionsRegionDurationSum = [AllsessionsRegionDurationSum; AllsessionsRegionDurationManMaxCurrPart; ...
        AllsessionsRegionDurationSideCurrPart; AllsessionsRegionDurationOLBCurrPart];
    AllregionsSum = [AllregionsSum; AllregionsManMaxCurrPart; AllregionsSideCurrPart; AllregionsOLBCurrPart];
            
    
    close all
    
end
%%%%


brushingSampCurrPart = brushingSampleCountsAll{11};

OcclusalRegs = {'ManRO', 'ManLO', 'MaxRO', 'MaxLO'};
idxOcc = regionNamestoNo(OcclusalRegs);
brushingSampOcc = sum(brushingSampCurrPart(idxOcc,:),1)';
LingualRegs = {'ManRL', 'ManAL', 'ManLL', 'MaxRL', 'MaxAL', 'MaxLL'};
idxLin = regionNamestoNo(LingualRegs);
brushingSampLin = sum(brushingSampCurrPart(idxLin,:),1)';
BuccalRegs = {'ManRB', 'ManAB', 'ManLB', 'MaxRB', 'MaxAB', 'MaxLB'};
idxBucc = regionNamestoNo(BuccalRegs);
brushingSampBucc = sum(brushingSampCurrPart(idxBucc,:),1)';


AllsessionsRegionDurationOLBCurrPart = [brushingSampOcc; brushingSampLin; brushingSampBucc];
AllregionsOLBCurrPart = [repmat({'Occlusal'}, length(brushingSampOcc),1); ...
    repmat({'Lingual'}, length(brushingSampLin),1); repmat({'Buccal'}, length(brushingSampBucc),1)];   
figure; subplot(1,2,1); boxplot(AllsessionsRegionDurationOLBCurrPart/sampRate, AllregionsOLBCurrPart);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,80])

brushingSampCurrPart = brushingSampleCountsAll{9};

OcclusalRegs = {'ManRO', 'ManLO', 'MaxRO', 'MaxLO'};
idxOcc = regionNamestoNo(OcclusalRegs);
brushingSampOcc = sum(brushingSampCurrPart(idxOcc,:),1)';
LingualRegs = {'ManRL', 'ManAL', 'ManLL', 'MaxRL', 'MaxAL', 'MaxLL'};
idxLin = regionNamestoNo(LingualRegs);
brushingSampLin = sum(brushingSampCurrPart(idxLin,:),1)';
BuccalRegs = {'ManRB', 'ManAB', 'ManLB', 'MaxRB', 'MaxAB', 'MaxLB'};
idxBucc = regionNamestoNo(BuccalRegs);
brushingSampBucc = sum(brushingSampCurrPart(idxBucc,:),1)';


AllsessionsRegionDurationOLBCurrPart = [brushingSampOcc; brushingSampLin; brushingSampBucc];
AllregionsOLBCurrPart = [repmat({'Occlusal'}, length(brushingSampOcc),1); ...
    repmat({'Lingual'}, length(brushingSampLin),1); repmat({'Buccal'}, length(brushingSampBucc),1)];   
subplot(1,2,2); boxplot(AllsessionsRegionDurationOLBCurrPart/sampRate, AllregionsOLBCurrPart);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,80])



%%%%

brushingSampCurrPart = brushingSampleCountsAll{4};

RightRegs = {'ManRO', 'ManRB', 'ManRL', 'MaxRO', 'MaxRB', 'MaxRL'};
idxRight = regionNamestoNo(RightRegs);
brushingSampRight = sum(brushingSampCurrPart(idxRight,:),1)';
AntRegs = {'ManAB', 'ManAL', 'MaxAB', 'MaxAL'};
idxAnt = regionNamestoNo(AntRegs);
brushingSampAnt = sum(brushingSampCurrPart(idxAnt,:),1)';
LeftRegs = {'ManLO', 'ManLB', 'ManLL', 'MaxLO', 'MaxLB', 'MaxLL'};
idxLeft = regionNamestoNo(LeftRegs);
brushingSampLeft = sum(brushingSampCurrPart(idxLeft,:),1)';

AllsessionsRegionDurationSideCurrPart = [brushingSampRight; brushingSampAnt; brushingSampLeft];
AllregionsSideCurrPart = [repmat({'Right'}, length(brushingSampRight),1); ...
    repmat({'Anterior'}, length(brushingSampAnt),1); repmat({'Left'}, length(brushingSampLeft),1)];   
figure; subplot(1,2,1); boxplot(AllsessionsRegionDurationSideCurrPart/sampRate, AllregionsSideCurrPart);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,50])  


brushingSampCurrPart = brushingSampleCountsAll{7};

RightRegs = {'ManRO', 'ManRB', 'ManRL', 'MaxRO', 'MaxRB', 'MaxRL'};
idxRight = regionNamestoNo(RightRegs);
brushingSampRight = sum(brushingSampCurrPart(idxRight,:),1)';
AntRegs = {'ManAB', 'ManAL', 'MaxAB', 'MaxAL'};
idxAnt = regionNamestoNo(AntRegs);
brushingSampAnt = sum(brushingSampCurrPart(idxAnt,:),1)';
LeftRegs = {'ManLO', 'ManLB', 'ManLL', 'MaxLO', 'MaxLB', 'MaxLL'};
idxLeft = regionNamestoNo(LeftRegs);
brushingSampLeft = sum(brushingSampCurrPart(idxLeft,:),1)';

AllsessionsRegionDurationSideCurrPart = [brushingSampRight; brushingSampAnt; brushingSampLeft];
AllregionsSideCurrPart = [repmat({'Right'}, length(brushingSampRight),1); ...
    repmat({'Anterior'}, length(brushingSampAnt),1); repmat({'Left'}, length(brushingSampLeft),1)];   
subplot(1,2,2); boxplot(AllsessionsRegionDurationSideCurrPart/sampRate, AllregionsSideCurrPart);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,50])


%%%
brushingSampCurrPart = brushingSampleCountsAll{4};

ManRegs = {'ManRO', 'ManRB', 'ManRL', 'ManAB', 'ManAL', 'ManLO', 'ManLB', 'ManLL'};
idxMan = regionNamestoNo(ManRegs);
brushingSampMan = sum(brushingSampCurrPart(idxMan,:),1)';
MaxRegs = {'MaxRO', 'MaxRB', 'MaxRL', 'MaxAB', 'MaxAL', 'MaxLO', 'MaxLB', 'MaxLL'};
idxMax = regionNamestoNo(MaxRegs);
brushingSampMax = sum(brushingSampCurrPart(idxMax,:),1)';

AllsessionsRegionDurationManMaxCurrPart = [brushingSampMan; brushingSampMax];
AllregionsManMaxCurrPart = [repmat({'Mandibular'}, length(brushingSampMan),1); repmat({'Maxillary'}, length(brushingSampMax),1)];   
figure; subplot(1,2,1); boxplot(AllsessionsRegionDurationManMaxCurrPart/sampRate, AllregionsManMaxCurrPart);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,90])
  
brushingSampCurrPart = brushingSampleCountsAll{7};

ManRegs = {'ManRO', 'ManRB', 'ManRL', 'ManAB', 'ManAL', 'ManLO', 'ManLB', 'ManLL'};
idxMan = regionNamestoNo(ManRegs);
brushingSampMan = sum(brushingSampCurrPart(idxMan,:),1)';
MaxRegs = {'MaxRO', 'MaxRB', 'MaxRL', 'MaxAB', 'MaxAL', 'MaxLO', 'MaxLB', 'MaxLL'};
idxMax = regionNamestoNo(MaxRegs);
brushingSampMax = sum(brushingSampCurrPart(idxMax,:),1)';

AllsessionsRegionDurationManMaxCurrPart = [brushingSampMan; brushingSampMax];
AllregionsManMaxCurrPart = [repmat({'Mandibular'}, length(brushingSampMan),1); repmat({'Maxillary'}, length(brushingSampMax),1)];   
subplot(1,2,2); boxplot(AllsessionsRegionDurationManMaxCurrPart/sampRate, AllregionsManMaxCurrPart);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,90])    
 
    

%%%%%%
brushingSampCurrPart = brushingSampleCountsAll{2};

OcclusalRegs = {'ManRO', 'ManLO', 'MaxRO', 'MaxLO'};
idxOcc = regionNamestoNo(OcclusalRegs);
brushingSampOcc = sum(brushingSampCurrPart(idxOcc,:),1)';
LingualRegs = {'ManRL', 'ManAL', 'ManLL', 'MaxRL', 'MaxAL', 'MaxLL'};
idxLin = regionNamestoNo(LingualRegs);
brushingSampLin = sum(brushingSampCurrPart(idxLin,:),1)';
BuccalRegs = {'ManRB', 'ManAB', 'ManLB', 'MaxRB', 'MaxAB', 'MaxLB'};
idxBucc = regionNamestoNo(BuccalRegs);
brushingSampBucc = sum(brushingSampCurrPart(idxBucc,:),1)';


AllsessionsRegionDurationOLBCurrPart = [brushingSampOcc; brushingSampLin; brushingSampBucc];
AllregionsOLBCurrPart = [repmat({'Occlusal'}, length(brushingSampOcc),1); ...
    repmat({'Lingual'}, length(brushingSampLin),1); repmat({'Buccal'}, length(brushingSampBucc),1)];   
figure; boxplot(AllsessionsRegionDurationOLBCurrPart/sampRate, AllregionsOLBCurrPart);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
% ylim([-2,80])



%% Population-level continued
%%% all dental surfaces
figure;boxplot(AllsessionsRegionDuration/sampRate, Allregions);
xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,40])
% title("Population-level brushing duration spent on each dental surface (in seconds)");


% %%% Maxillary vs. Mandibular
% AllregionsUD = repmat({'Mandibular'}, length(Allregions),1);
% idxUp = find(Allregions == 'MaxRO' | Allregions == 'MaxRB' | Allregions == 'MaxRL' | Allregions == 'MaxAB' | ...
% Allregions == 'MaxAL' | Allregions == 'MaxLO' | Allregions == 'MaxLB' | Allregions == 'MaxLL');
% AllregionsUD(idxUp) = repmat({'Maxillary'}, length(idxUp),1);
% figure;boxplot(AllsessionsRegionDuration/sampRate, AllregionsUD);
% xlabel('Dental surfaces');
% ylabel('Duration (in seconds)');
% % ylim([-2,42])
% % title("Population-level time spent on maxillary and mandibular dental surfaces (in seconds)");
% 
% %%% Right vs Left vs Anterior
% AllregionsSide = repmat({'Right'}, length(Allregions),1);
% idxLeft = find(Allregions == 'ManLO' | Allregions == 'ManLB' | Allregions == 'ManLL' | Allregions == 'MaxLO' | ...
% Allregions == 'MaxLB' | Allregions == 'MaxLL');
% AllregionsSide(idxLeft) = repmat({'Left'}, length(idxLeft),1);
% idxAnt = find(Allregions == 'ManAB' | Allregions == 'ManAL' | Allregions == 'MaxAB' | Allregions == 'MaxAL');
% AllregionsSide(idxAnt) = repmat({'Anterior'}, length(idxAnt),1);
% figure;boxplot(AllsessionsRegionDuration/sampRate, AllregionsSide);
% xlabel('Dental surfaces');
% ylabel('Duration (in seconds)');
% % ylim([-2,42])
% % title("Population-level time spent on right, left, and anterior dental surfaces (in seconds)");
% 
% %%% Occlusal vs Lingual vs Buccal
% AllregionsOLB = repmat({'Occlusal'}, length(Allregions),1);
% idxLing = find(Allregions == 'ManRL' | Allregions == 'ManAL' | Allregions == 'ManLL' | Allregions == 'MaxRL' | ...
% Allregions == 'MaxAL' | Allregions == 'MaxLL');
% AllregionsOLB(idxLing) = repmat({'Lingual'}, length(idxLing),1);
% idxBucc = find(Allregions == 'ManRB' | Allregions == 'ManAB' | Allregions == 'ManLB' | Allregions == 'MaxRB' | ...
% Allregions == 'MaxAB' | Allregions == 'MaxLB');
% AllregionsOLB(idxBucc) = repmat({'Buccal'}, length(idxBucc),1);
% figure;boxplot(AllsessionsRegionDuration/sampRate, AllregionsOLB);
% xlabel('Dental surfaces');
% ylabel('Duration (in seconds)');
% % ylim([-2,42])
% % title("Population-level time spent on occlusal, lingual, and buccal dental surfaces (in seconds)");



%%% Mandibular vs Maxillary
figure;
% tiledlayout(1,3);
% ax1 = nexttile;
subplot(1,3,1);
idxSumManMax = find(AllregionsSum == 'Mandibular' | AllregionsSum == 'Maxillary');
boxplot(AllsessionsRegionDurationSum(idxSumManMax)/sampRate, reorderlevels(nominal(AllregionsSum(idxSumManMax)), {'Mandibular', 'Maxillary'}));
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,97])
% title("Population-level time spent on mandibular and maxillary dental surfaces (in seconds)");

%%% Right vs Left vs Anterior
% ax2 = nexttile;
subplot(1,3,2);
idxSumSide = find(AllregionsSum == 'Right' | AllregionsSum == 'Anterior' | AllregionsSum == 'Left' );
boxplot(AllsessionsRegionDurationSum(idxSumSide)/sampRate, reorderlevels(nominal(AllregionsSum(idxSumSide)), {'Right', 'Anterior', 'Left'}));
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,97])
% title("Population-level time spent on right, anterior, and left dental surfaces (in seconds)");

%%% Occlusal vs Lingual vs Buccal
% ax3 = nexttile;
subplot(1,3,3);
idxSumOLB = find(AllregionsSum == 'Occlusal' | AllregionsSum == 'Lingual' | AllregionsSum == 'Buccal');
boxplot(AllsessionsRegionDurationSum(idxSumOLB)/sampRate, reorderlevels(nominal(AllregionsSum(idxSumOLB)), {'Occlusal', 'Lingual', 'Buccal'}));
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
ylim([-2,97])
% title("Population-level time spent on occlusal, lingual, and buccal dental surfaces (in seconds)");

%% Prameter of interest: Duration of brushing with extra pressure on each region

%% Population level
AllsessionsRegionDurationPress = [];
AllregionsPress = categorical([]);

AllsessionsRegionDurationPressSum = [];
AllregionsPressSum = categorical([]);
%% Participant-level and Session-level
%%% Figure 6
for i = 1:length(ids)
    AllsessionsRegionDurationPressCurrPart = reshape(pressureSampleCountsAll{i}, ...
        size(pressureSampleCountsAll{i},1)*size(pressureSampleCountsAll{i},2),1);
    AllregionsPressCurrPart = repmat(categorical(region_names_new'),size(pressureSampleCountsAll{i},2),1);
    
%     figure;boxplot(AllsessionsRegionDurationPressCurrPart/sampRate, AllregionsPressCurrPart);
%     xlabel('Dental surfaces');
%     ylabel('Duration (in seconds)');
% %     title("Time spent by " + "P" + num2str(i) + " brushing with excessive pressure on each dental surfaces (in seconds)");
%     
%     %%% Maxillary vs. Mandibular
%     AllregionsPressCurrPartUD = repmat({'Mandibular'}, length(AllregionsPressCurrPart),1);
%     idxPressCurrPartUp = find(AllregionsPressCurrPart == 'MaxRO' | AllregionsPressCurrPart == 'MaxRB' | AllregionsPressCurrPart == 'MaxRL' | AllregionsPressCurrPart == 'MaxAB' | ...
%     AllregionsPressCurrPart == 'MaxAL' | AllregionsPressCurrPart == 'MaxLO' | AllregionsPressCurrPart == 'MaxLB' | AllregionsPressCurrPart == 'MaxLL');
%     AllregionsPressCurrPartUD(idxPressCurrPartUp) = repmat({'Maxillary'}, length(idxPressCurrPartUp),1);
%     figure;boxplot(AllsessionsRegionDurationPressCurrPart/sampRate, AllregionsPressCurrPartUD);
%     xlabel('Dental surfaces');
%     ylabel('Duration (in seconds)');
% %     ylim([-2,42])
% %     title("Time spent by " + "P" + num2str(i) + " brushing with excessive pressure on maxillary and mandibular dental surfaces (in seconds)");
%    
%     %%% Right vs Left vs Anterior
%     AllregionsPressCurrPartSide = repmat({'Right'}, length(AllregionsPressCurrPart),1);
%     idxPressCurrPartLeft = find(AllregionsPressCurrPart == 'ManLO' | AllregionsPressCurrPart == 'ManLB' | AllregionsPressCurrPart == 'ManLL' | AllregionsPressCurrPart == 'MaxLO' | ...
%     AllregionsPressCurrPart == 'MaxLB' | AllregionsPressCurrPart == 'MaxLL');
%     AllregionsPressCurrPartSide(idxPressCurrPartLeft) = repmat({'Left'}, length(idxPressCurrPartLeft),1);
%     idxPressCurrPartAnt = find(AllregionsPressCurrPart == 'ManAB' | AllregionsPressCurrPart == 'ManAL' | AllregionsPressCurrPart == 'MaxAB' | AllregionsPressCurrPart == 'MaxAL');
%     AllregionsPressCurrPartSide(idxPressCurrPartAnt) = repmat({'Anterior'}, length(idxPressCurrPartAnt),1);
%     figure;boxplot(AllsessionsRegionDurationPressCurrPart/sampRate, AllregionsPressCurrPartSide);
%     xlabel('Dental surfaces');
%     ylabel('Duration (in seconds)');
% %     ylim([-2,42])
% %     title("Time spent by " + "P" + num2str(i) + " brushing with excessive pressure on right, left, and anterior dental surfaces (in seconds)");
% 
%     %%% Occlusal vs Lingual vs Buccal
%     AllregionsPressOLB = repmat({'Occlusal'}, length(AllregionsPressCurrPart),1);
%     idxPressCurrPartLing = find(AllregionsPressCurrPart == 'ManRL' | AllregionsPressCurrPart == 'ManAL' | AllregionsPressCurrPart == 'ManLL' | AllregionsPressCurrPart == 'MaxRL' | ...
%     AllregionsPressCurrPart == 'MaxAL' | AllregionsPressCurrPart == 'MaxLL');
%     AllregionsPressOLB(idxPressCurrPartLing) = repmat({'Lingual'}, length(idxPressCurrPartLing),1);
%     idxPressCurrPartBucc = find(AllregionsPressCurrPart == 'ManRB' | AllregionsPressCurrPart == 'ManAB' | AllregionsPressCurrPart == 'ManLB' | AllregionsPressCurrPart == 'MaxRB' | ...
%     AllregionsPressCurrPart == 'MaxAB' | AllregionsPressCurrPart == 'MaxLB');
%     AllregionsPressOLB(idxPressCurrPartBucc) = repmat({'Buccal'}, length(idxPressCurrPartBucc),1);
%     figure;boxplot(AllsessionsRegionDurationPressCurrPart/sampRate, AllregionsPressOLB);
%     xlabel('Dental surfaces');
%     ylabel('Duration (in seconds)');
% %     ylim([-2,42])
% %     title("Time spent by " + "P" + num2str(i) + " brushing with excessive pressure occlusal, lingual, and buccal dental surfaces (in seconds)");


    AllsessionsRegionDurationPress = [AllsessionsRegionDurationPress; AllsessionsRegionDurationPressCurrPart];
    AllregionsPress = [AllregionsPress; AllregionsPressCurrPart];
    
    

    %%% Maxillary vs. Mandibular


    pressureSampCurrPart = pressureSampleCountsAll{i};
    
    ManRegs = {'ManRO', 'ManRB', 'ManRL', 'ManAB', 'ManAL', 'ManLO', 'ManLB', 'ManLL'};
    idxMan = regionNamestoNo(ManRegs);
    pressureSampMan = sum(pressureSampCurrPart(idxMan,:),1)';
    MaxRegs = {'MaxRO', 'MaxRB', 'MaxRL', 'MaxAB', 'MaxAL', 'MaxLO', 'MaxLB', 'MaxLL'};
    idxMax = regionNamestoNo(MaxRegs);
    pressureSampMax = sum(pressureSampCurrPart(idxMax,:),1)';
    
    AllsessionsRegionDurationManMaxPressCurrPart = [pressureSampMan; pressureSampMax];
    AllregionsManMaxPressCurrPart = [repmat({'Mandibular'}, length(pressureSampMan),1); repmat({'Maxillary'}, length(pressureSampMax),1)];   
    figure;boxplot(AllsessionsRegionDurationManMaxPressCurrPart/sampRate, AllregionsManMaxPressCurrPart, 'whisker', 1.5);
    xlabel('Dental surfaces');
    ylabel('Duration (in seconds)');
    
    
       
    %%% Right vs Left vs Anterior

    RightRegs = {'ManRO', 'ManRB', 'ManRL', 'MaxRO', 'MaxRB', 'MaxRL'};
    idxRight = regionNamestoNo(RightRegs);
    pressureSampRight = sum(pressureSampCurrPart(idxRight,:),1)';
    AntRegs = {'ManAB', 'ManAL', 'MaxAB', 'MaxAL'};
    idxAnt = regionNamestoNo(AntRegs);
    pressureSampAnt = sum(pressureSampCurrPart(idxAnt,:),1)';
    LeftRegs = {'ManLO', 'ManLB', 'ManLL', 'MaxLO', 'MaxLB', 'MaxLL'};
    idxLeft = regionNamestoNo(LeftRegs);
    pressureSampLeft = sum(pressureSampCurrPart(idxLeft,:),1)';
    
    AllsessionsRegionDurationSidePressCurrPart = [pressureSampRight; pressureSampAnt; pressureSampLeft];
    AllregionsSidePressCurrPart = [repmat({'Right'}, length(pressureSampRight),1); ...
        repmat({'Anterior'}, length(pressureSampAnt),1); repmat({'Left'}, length(pressureSampLeft),1)];   
    figure;boxplot(AllsessionsRegionDurationSidePressCurrPart/sampRate, AllregionsSidePressCurrPart, 'whisker', 1.5);
    xlabel('Dental surfaces');
    ylabel('Duration (in seconds)');
       
    %%% Occlusal vs Lingual vs Buccal

    OcclusalRegs = {'ManRO', 'ManLO', 'MaxRO', 'MaxLO'};
    idxOcc = regionNamestoNo(OcclusalRegs);
    pressureSampOcc = sum(pressureSampCurrPart(idxOcc,:),1)';
    LingualRegs = {'ManRL', 'ManAL', 'ManLL', 'MaxRL', 'MaxAL', 'MaxLL'};
    idxLin = regionNamestoNo(LingualRegs);
    pressureSampLin = sum(pressureSampCurrPart(idxLin,:),1)';
    BuccalRegs = {'ManRB', 'ManAB', 'ManLB', 'MaxRB', 'MaxAB', 'MaxLB'};
    idxBucc = regionNamestoNo(BuccalRegs);
    pressureSampBucc = sum(pressureSampCurrPart(idxBucc,:),1)';
    
    
    AllsessionsRegionDurationOLBPressCurrPart = [pressureSampOcc; pressureSampLin; pressureSampBucc];
    AllregionsOLBPressCurrPart = [repmat({'Occlusal'}, length(pressureSampOcc),1); ...
        repmat({'Lingual'}, length(pressureSampLin),1); repmat({'Buccal'}, length(pressureSampBucc),1)];   
    figure;boxplot(AllsessionsRegionDurationOLBPressCurrPart/sampRate, AllregionsOLBPressCurrPart, 'whisker', 1.5);
    xlabel('Dental surfaces');
    ylabel('Duration (in seconds)');
    
    
    AllsessionsRegionDurationPressSum = [AllsessionsRegionDurationPressSum; AllsessionsRegionDurationManMaxPressCurrPart; ...
        AllsessionsRegionDurationSidePressCurrPart; AllsessionsRegionDurationOLBPressCurrPart];
    AllregionsPressSum = [AllregionsPressSum; AllregionsManMaxPressCurrPart; AllregionsSidePressCurrPart; AllregionsOLBPressCurrPart];
    
    close all;
     
end
%%%%


AllsessionsRegionDurationPressCurrPart = reshape(pressureSampleCountsAll{6}, ...
    size(pressureSampleCountsAll{6},1)*size(pressureSampleCountsAll{6},2),1);
AllregionsPressCurrPart = repmat(categorical(region_names_new'),size(pressureSampleCountsAll{6},2),1);

figure; subplot(1,2,1); boxplot(AllsessionsRegionDurationPressCurrPart/sampRate, AllregionsPressCurrPart);
%     xlabel('Dental surfaces');
ylabel('Duration (in seconds)');


AllsessionsRegionDurationPressCurrPart = reshape(pressureSampleCountsAll{11}, ...
    size(pressureSampleCountsAll{11},1)*size(pressureSampleCountsAll{11},2),1);
AllregionsPressCurrPart = repmat(categorical(region_names_new'),size(pressureSampleCountsAll{11},2),1);

subplot(1,2,2); boxplot(AllsessionsRegionDurationPressCurrPart/sampRate, AllregionsPressCurrPart);
%     xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
   

%% pressure in each session

M = movsum(AllsessionsRegionDuration, [0,15]);
indxMult = floor(length(AllsessionsRegionDuration)/16)-1;
brushingSumSession = M(16*[0:indxMult] + 1);

MPress = movsum(AllsessionsRegionDurationPress, [0,15]);
indxMultPress = floor(length(AllsessionsRegionDurationPress)/16)-1;
pressureSumSession = MPress(16*[0:indxMultPress] + 1);

%% Population-level continued
%%% all dental surfaces
figure;boxplot(AllsessionsRegionDurationPress/sampRate, AllregionsPress, 'whisker', 1.5);
xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
% title("Population-level brushing duration with excessive pressure on each dental surface (in seconds)");


% %%% Maxillary vs. Mandibular
% AllregionsPressUD = repmat({'Mandibular'}, length(AllregionsPress),1);
% idxPressUp = find(AllregionsPress == 'MaxRO' | AllregionsPress == 'MaxRB' | AllregionsPress == 'MaxRL' | AllregionsPress == 'MaxAB' | ...
% AllregionsPress == 'MaxAL' | AllregionsPress == 'MaxLO' | AllregionsPress == 'MaxLB' | AllregionsPress == 'MaxLL');
% AllregionsPressUD(idxPressUp) = repmat({'Maxillary'}, length(idxPressUp),1);
% figure;boxplot(AllsessionsRegionDurationPress/sampRate, AllregionsPressUD);
% xlabel('Dental surfaces');
% ylabel('Duration (in seconds)');
% % ylim([-2,42])
% % title("Population-level time spent with excessive pressure on maxillary and mandibular dental surfaces (in seconds)");
% 
% %%% Right vs Left vs Anterior
% AllregionsPressSide = repmat({'Right'}, length(AllregionsPress),1);
% idxPressLeft = find(AllregionsPress == 'ManLO' | AllregionsPress == 'ManLB' | AllregionsPress == 'ManLL' | AllregionsPress == 'MaxLO' | ...
% AllregionsPress == 'MaxLB' | AllregionsPress == 'MaxLL');
% AllregionsPressSide(idxPressLeft) = repmat({'Left'}, length(idxPressLeft),1);
% idxPressAnt = find(AllregionsPress == 'ManAB' | AllregionsPress == 'ManAL' | AllregionsPress == 'MaxAB' | AllregionsPress == 'MaxAL');
% AllregionsPressSide(idxPressAnt) = repmat({'Anterior'}, length(idxPressAnt),1);
% figure;boxplot(AllsessionsRegionDurationPress/sampRate, AllregionsPressSide);
% xlabel('Dental surfaces');
% ylabel('Duration (in seconds)');
% % ylim([-2,42])
% % title("Population-level time spent with excessive pressure on right, left, and anterior dental surfaces (in seconds)");
% 
% %%% Occlusal vs Lingual vs Buccal
% AllregionsPressOLB = repmat({'Occlusal'}, length(AllregionsPress),1);
% idxPressLing = find(AllregionsPress == 'ManRL' | AllregionsPress == 'ManAL' | AllregionsPress == 'ManLL' | AllregionsPress == 'MaxRL' | ...
% AllregionsPress == 'MaxAL' | AllregionsPress == 'MaxLL');
% AllregionsPressOLB(idxPressLing) = repmat({'Lingual'}, length(idxPressLing),1);
% idxPressBucc = find(AllregionsPress == 'ManRB' | AllregionsPress == 'ManAB' | AllregionsPress == 'ManLB' | AllregionsPress == 'MaxRB' | ...
% AllregionsPress == 'MaxAB' | AllregionsPress == 'MaxLB');
% AllregionsPressOLB(idxPressBucc) = repmat({'Buccal'}, length(idxPressBucc),1);
% figure;boxplot(AllsessionsRegionDurationPress/sampRate, AllregionsPressOLB);
% xlabel('Dental surfaces');
% ylabel('Duration (in seconds)');
% % ylim([-2,42])
% % title("Population-level time spent with excessive pressure on occlusal, lingual, and buccal dental surfaces (in seconds)");


%%% Mandibular vs Maxillary
figure;
% tiledlayout(1,3);
% ax1 = nexttile;
subplot(1,3,1);
idxSumManMaxPress = find(AllregionsPressSum == 'Mandibular' | AllregionsPressSum == 'Maxillary');
boxplot(AllsessionsRegionDurationPressSum(idxSumManMaxPress)/sampRate, ...
    reorderlevels(nominal(AllregionsPressSum(idxSumManMaxPress)), {'Mandibular', 'Maxillary'}), 'whisker', 1.5);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
% ylim([-2,97])
% title("Population-level time spent with excessive pressure on mandibular and maxillary dental surfaces (in seconds)");

%%% Right vs Left vs Anterior
% ax2 = nexttile;
subplot(1,3,2);
idxSumSidePress = find(AllregionsPressSum == 'Right' | AllregionsPressSum == 'Anterior' | AllregionsPressSum == 'Left' );
boxplot(AllsessionsRegionDurationPressSum(idxSumSidePress)/sampRate, reorderlevels(nominal(AllregionsPressSum(idxSumSidePress)), ...
    {'Right', 'Anterior', 'Left'}), 'whisker', 1.5);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
% ylim([-2,97])
% title("Population-level time spentwith excessuve pressure on right, anterior, and left dental surfaces (in seconds)");

%%% Occlusal vs Lingual vs Buccal
% ax3 = nexttile;
subplot(1,3,3);
idxSumOLBPress = find(AllregionsPressSum == 'Occlusal' | AllregionsPressSum == 'Lingual' | AllregionsPressSum == 'Buccal');
boxplot(AllsessionsRegionDurationPressSum(idxSumOLBPress)/sampRate, reorderlevels(nominal(AllregionsPressSum(idxSumOLBPress)), ...
    {'Occlusal', 'Lingual', 'Buccal'}), 'whisker', 1.5);
% xlabel('Dental surfaces');
ylabel('Duration (in seconds)');
% ylim([-2,97])
% title("Population-level time spent with excessive pressure on occlusal, lingual, and buccal dental surfaces (in seconds)");

%% Modeling brushing region
FixedglmeBrushRegDur = fitglme(tableRegionBrushingEachReg([2:end, 1], :),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB ',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[FixedbetaBrush,FixedbetanamesBrush, statsFixedbetaBrush] = fixedEffects(FixedglmeBrushRegDur);
[FixedBBrushRegDur,FixedBNamesBrushRegDur] = randomEffects(FixedglmeBrushRegDur);


%%% appendix section 1
glmeBrushRegDur = fitglme(tableRegionBrushingEachReg([2:end, 1], :),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1|Participant) + (1|Participant:Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaBrush,betanamesBrush, statsbetaBrush] = fixedEffects(glmeBrushRegDur);
[BBrushRegDur,BNamesBrushRegDur, statsBrushRegDur] = randomEffects(glmeBrushRegDur);
%%%%

glmeBrushRegDur2 = fitglme(tableRegionBrushingEachReg([2:end, 1], :),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1 + MaxOrNot + Side + OLB | Participant) + (1|Participant:Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaBrush2, betanamesBrush2, statsbetaBrush2] = fixedEffects(glmeBrushRegDur2);
[BBrushRegDur2, BNamesBrushRegDur2, statsBrushRegDur2] = randomEffects(glmeBrushRegDur2);


glmeBrushRegDur3 = fitglme(tableRegionBrushingEachReg([2:end, 1], :),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (-1 + MaxOrNot|Participant) + (-1 + Side|Participant) +  (-1 + OLB|Participant) + (1|Participant) + (1|Participant:Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaBrush3,betanamesBrush3, statsbetaRegDur3] = fixedEffects(glmeBrushRegDur3);
[BBrushRegDur3,BNamesBrushRegDur3, statsBrushRegDur3] = randomEffects(glmeBrushRegDur3);
% E(exp(N(0,0.22)))
%%% Zero inflated Poisson regression

%%%%%%% Appendix section 1
MaxOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.MaxOrNot)=='Max')==2;
AntOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.Side)=='A')==2;
LeftOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.Side)=='L')==2;
% BuccalOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.OLB)=='B')==2;
LingualOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.OLB)=='L')==2;
OcclusalOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.OLB)=='O')==2;

resBrush = EMzeropoisson_mat([tableRegionBrushingEachReg.sampleCounts MaxOrNotValuesBrush ...
    AntOrNotValuesBrush LeftOrNotValuesBrush LingualOrNotValuesBrush OcclusalOrNotValuesBrush],0.5,'initialtau','stable','tol',0.001);
fittedParams = resBrush(end,:);

%%% finding p-values using permutation test
no_iter = 1000;
fittedParamsPermut = zeros(no_iter,7);

% for i = 1:no_iter
%     randomPerm = randperm(length(tableRegionBrushingEachReg.sampleCounts));
%     samplesPerm = tableRegionBrushingEachReg.sampleCounts;
%     samplesPerm = samplesPerm(randomPerm);
%     resBrush = EMzeropoisson_mat([samplesPerm MaxOrNotValuesBrush ...
%         AntOrNotValuesBrush LeftOrNotValuesBrush LingualOrNotValuesBrush OcclusalOrNotValuesBrush]);
%     fittedParamsPermut(i,:) = resBrush(end,:);
% end

variates = [MaxOrNotValuesBrush, AntOrNotValuesBrush, LeftOrNotValuesBrush, LingualOrNotValuesBrush, OcclusalOrNotValuesBrush];   

    
for i = 1:no_iter               
    for j = 1:size(variates,2)
        variatesCurr = variates(:,j);
        variatesPerm = variates;

        randomPerm = randperm(length(variatesCurr));
        variatesCurrPerm = variatesCurr(randomPerm);

        variatesPerm(:,j) = variatesCurrPerm;
        resBrush = EMzeropoisson_mat([tableRegionBrushingEachReg.sampleCounts variatesPerm]);
        fittedParamsPermut(i,j+2) = resBrush(end,j+2);
    end
end



pValuesBrushingNew = sum(abs(fittedParams) < abs(fittedParamsPermut),1)/no_iter

% pValuesBrushing = min([sum(fittedParams > fittedParamsPermut, 1); sum(fittedParams < fittedParamsPermut, 1)])/no_iter

for i = 1:size(fittedParamsPermut,2)
    figure; histogram(fittedParamsPermut(:,i),100); 
    xlabel('Values');
    ylabel('Magnitude');
    title("Distribution of coefficient #" + num2str(i));
end

%%%%%%%%

no_bootstraps = 1000;
fittedParamsBS = zeros(no_bootstraps,7);

    
no_participants = length(unique(tableRegionBrushingEachReg.Participant));

for n = 1:no_bootstraps

    bootsrapBrushEachReg = table;
    
    MaxOrNotValuesBrushBS = [];
    AntOrNotValuesBrushBS = [];
    LeftOrNotValuesBrushBS = [];
    LingualOrNotValuesBrushBS = [];
    OcclusalOrNotValuesBrushBS = [];
            

    for i = 1:no_participants
        randParticipant = randsample(no_participants, 1);        
           
        no_sessions = length(unique(tableRegionBrushingEachReg.Session(tableRegionBrushingEachReg.Participant==randParticipant)));
        idxBS = find(tableRegionBrushingEachReg.Participant == randParticipant & tableRegionBrushingEachReg.Session == randSession);
        
        for j = 1:no_sessions
            randSession = randsample(no_sessions, 1);     
            bootsrapBrushEachReg = [bootsrapBrushEachReg; tableRegionBrushingEachReg( ...
                idxBS, :)];   
            
            MaxOrNotValuesBrushBS = [MaxOrNotValuesBrushBS; MaxOrNotValuesBrush(idxBS)];
            AntOrNotValuesBrushBS = [AntOrNotValuesBrushBS; AntOrNotValuesBrush(idxBS)];
            LeftOrNotValuesBrushBS = [LeftOrNotValuesBrushBS; LeftOrNotValuesBrush(idxBS)];
            LingualOrNotValuesBrushBS = [LingualOrNotValuesBrushBS; LingualOrNotValuesBrush(idxBS)];
            OcclusalOrNotValuesBrushBS = [OcclusalOrNotValuesBrushBS; OcclusalOrNotValuesBrush(idxBS)];
            
        end
    end
    
    resBrushBS = EMzeropoisson_mat([bootsrapBrushEachReg.sampleCounts MaxOrNotValuesBrushBS ...
        AntOrNotValuesBrushBS LeftOrNotValuesBrushBS LingualOrNotValuesBrushBS OcclusalOrNotValuesBrushBS]);
    fittedParamsBS(n,:) = resBrushBS(end,:);
end


N = n;   %100  

CIBrushing = zeros(size(fittedParamsBS,2),2); 
SEBrushing = zeros(size(fittedParamsBS,2),1); 
for i = 1:size(fittedParamsBS,2)
%     figure; histogram(fittedParamsBS(1:N,i),100); 
%     xlabel('Values');
%     ylabel('Magnitude');
%     title("Distribution of coefficient #" + num2str(i));
    
    CIBrushing(i,:) = [quantile(fittedParamsBS(1:N,i), 0.025) , quantile(fittedParamsBS(1:N,i), 0.975)];
%     SEBrushing(i) = (CIBrushing(i,2) - CIBrushing(i,1))/3.92;
    SEBrushing(i) = std(fittedParamsBS(1:N,i));                            
end



% CI95 = tinv([0.025 0.975], N-1);                

% meanBootstrapInterc = mean(fittedParamsBS(1:N,2),1);                             
% ySEMInterc = std(fittedParamsBS(1:N,2))/sqrt(N);                            
% yCI95Interc = bsxfun(@times, ySEMInterc, CI95(:));            
% yCI95Interc+meanBootstrapInterc


% meanBootstrap = mean(fittedParamsBS(1:N,2),1);                             
% ySEM = std(fittedParamsBS(1:N,2))/sqrt(N);                            
% yCI95 = bsxfun(@times, ySEM, CI95(:));            
% figure
% plot(1:size(fittedParamsBS,2), meanBootstrap)                                   
% hold on
% plot(1:size(fittedParamsBS,2), yCI95+meanBootstrap)                       
% hold off
% grid 
% yCI95+meanBootstrap




%% Modeling pressure region
FixedglmeBrushRegDurPress = fitglme(tableRegionPressureEachReg([1:end],:),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB ',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[FixedbetaPress,FixedbetanamesPress, statsFixedbetaPress] = fixedEffects(FixedglmeBrushRegDurPress);
[FixedBBrushRegDurPress,FixedBNamesBrushRegDurPress] = randomEffects(FixedglmeBrushRegDurPress);


%%% appendix section 1
glmePressRegDur = fitglme(tableRegionPressureEachReg([1:end],:),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1|Participant) + (1|Participant:Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaPress,betanamesPress, statsbetaPress] = fixedEffects(glmePressRegDur);
[BPressRegDur,BNamesPressRegDur, statsPressRegDur] = randomEffects(glmePressRegDur);
%%%%

glmePressRegDur2 = fitglme(tableRegionPressureEachReg([1:end],:),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1 + MaxOrNot + Side + OLB | Participant) + (1|Participant:Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaPress2, betanamesPress2, statsbetaPress2] = fixedEffects(glmePressRegDur2);
[BPressRegDur2, BNamesPressRegDur2, statsPressRegDur2] = randomEffects(glmePressRegDur2);


glmePressRegDur3 = fitglme(tableRegionPressureEachReg([1:end],:),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (-1 + MaxOrNot|Participant) + (-1 + Side|Participant) +  (-1 + OLB|Participant) + (1|Participant) + (1|Participant:Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaPress3,betanamesPress3, statsbetaPress3] = fixedEffects(glmePressRegDur3);
[BPressRegDur3,BNamesPressRegDur3,  statsPressRegDur3] = randomEffects(glmePressRegDur3);


%%% Zero inflated Poisson regression

%%%%%%% Appendix section 1
MaxOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.MaxOrNot)=='Max')==2;
AntOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Side)=='A')==2;
LeftOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Side)=='L')==2;
BuccalOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.OLB)=='B')==2;
LingualOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.OLB)=='L')==2;
% OcclusalOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.OLB)=='O')==2;

resPress = EMzeropoisson_mat([tableRegionPressureEachReg.sampleCounts MaxOrNotValuesPress ...
    AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesPress LingualOrNotValuesPress],0.5,'initialtau','stable','tol',0.001);
fittedParamsPress = resPress(end,:);

%%% finding p-values using permutation test
no_iter = 1000;
fittedParamsPermutPress = zeros(no_iter,7);
% for i = 1:no_iter
%     randomPerm = randperm(length(tableRegionPressureEachReg.sampleCounts));
%     samplesPerm = tableRegionPressureEachReg.sampleCounts;
%     samplesPerm = samplesPerm(randomPerm);
%     resPress = EMzeropoisson_mat([samplesPerm MaxOrNotValuesPress ...
%         AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesPress LingualOrNotValuesPress]);
%     fittedParamsPermutPress(i,:) = resPress(end,:);
% end

variatesPress = [MaxOrNotValuesPress, AntOrNotValuesPress, LeftOrNotValuesPress, BuccalOrNotValuesPress, LingualOrNotValuesPress];   

   
for i = 1:no_iter               
    for j = 1:size(variatesPress,2)
        variatesCurr = variatesPress(:,j);
        variatesPerm = variatesPress;

        randomPerm = randperm(length(variatesCurr));
        variatesCurrPerm = variatesCurr(randomPerm);

        variatesPerm(:,j) = variatesCurrPerm;
        resPress = EMzeropoisson_mat([tableRegionBrushingEachReg.sampleCounts variatesPerm]);
        fittedParamsPermutPress(i,j+2) = resPress(end,j+2);
    end
end


% what porportion of the absolute values of the permuted esitmates are larger than the absoloute values of the actual estimates
pValuesPressNew = sum(abs(fittedParamsPress) < abs(fittedParamsPermutPress),1)/no_iter

% pValuesPress = min([sum(fittedParamsPress > fittedParamsPermutPress, 1); sum(fittedParamsPress < fittedParamsPermutPress, 1)])/no_iter

for i = 1:size(fittedParamsPermutPress,2)
    figure; histogram(fittedParamsPermutPress(:,i),100); 
    xlabel('Values');
    ylabel('Magnitude');
    title("Distribution of coefficient #" + num2str(i));
end

%%%%%%%%

no_bootstraps = 1000;
fittedParamsPressBS = zeros(no_bootstraps,7);

    
no_participants = length(unique(tableRegionPressureEachReg.Participant));



for n = 1:no_bootstraps

    bootsrapPressEachReg = table;
    
    MaxOrNotValuesPressBS = [];
    AntOrNotValuesPressBS = [];
    LeftOrNotValuesPressBS = [];
    BuccalOrNotValuesPressBS = [];
    LingualOrNotValuesPressBS = [];


    for i = 1:no_participants
        randParticipant = randsample(no_participants, 1);        
           
        no_sessions = length(unique(tableRegionPressureEachReg.Session(tableRegionPressureEachReg.Participant==randParticipant)));
        
        for j = 1:no_sessions
            randSession = randsample(no_sessions, 1);     
            idxBS = find(tableRegionPressureEachReg.Participant == randParticipant & tableRegionPressureEachReg.Session == randSession);
            bootsrapPressEachReg = [bootsrapPressEachReg; tableRegionPressureEachReg( ...
                idxBS, :)];
            MaxOrNotValuesPressBS = [MaxOrNotValuesPressBS;MaxOrNotValuesPress(idxBS)];
            AntOrNotValuesPressBS = [AntOrNotValuesPressBS;AntOrNotValuesPress(idxBS)];
            LeftOrNotValuesPressBS = [LeftOrNotValuesPressBS;LeftOrNotValuesPress(idxBS)];
            BuccalOrNotValuesPressBS = [BuccalOrNotValuesPressBS;BuccalOrNotValuesPress(idxBS)];
            LingualOrNotValuesPressBS = [LingualOrNotValuesPressBS;LingualOrNotValuesPress(idxBS)];
           
        end
    end
    
    resPressBS = EMzeropoisson_mat([bootsrapPressEachReg.sampleCounts MaxOrNotValuesPressBS ...
        AntOrNotValuesPressBS LeftOrNotValuesPressBS BuccalOrNotValuesPressBS LingualOrNotValuesPressBS]);
    fittedParamsPressBS(n,:) = resPressBS(end,:);
end

if n < no_bootstraps
    N = n-1;   %100  
else
    N = n;
end

CIPressure = zeros(size(fittedParamsPressBS,2),2); 
SEPressure = zeros(size(fittedParamsPressBS,2),1); 
for i = 1:size(fittedParamsPressBS,2)
%     figure; histogram(fittedParamsPressBS(1:N,i),100); 
%     xlabel('Values');
%     ylabel('Magnitude');
%     title("Distribution of coefficient #" + num2str(i));
    
    CIPressure(i,:) = [quantile(fittedParamsPressBS(1:N,i), 0.025) , quantile(fittedParamsPressBS(1:N,i), 0.975)];
%     SEPressure(i) = (CIPressure(i,2) - CIPressure(i,1))/3.92;
    SEPressure(i) = std(fittedParamsPressBS(1:N,i));                            

end


1;

% CI95 = tinv([0.025 0.975], N-1);                

% meanPressBootstrapInterc = mean(fittedParamsPressBS(1:N,2),1);                             
% yPressSEMInterc = std(fittedParamsPressBS(1:N,2))/sqrt(N);                            
% yPressCI95Interc = bsxfun(@times, yPressSEMInterc, CI95(:));            
% yPressCI95Interc+meanPressBootstrapInterc


% meanPressBootstrap = mean(fittedParamsPressBS(1:N,2),1);                             
% yPressSEM = std(fittedParamsPressBS(1:N,2))/sqrt(N);                            
% yPressCI95 = bsxfun(@times, yPressSEM, CI95(:));            
% figure
% plot(1:size(fittedParamsPressBS,2), meanPressBootstrap)                                   
% hold on
% plot(1:size(fittedParamsPressBS,2), yCI95+meanPressBootstrap)                       
% hold off
% grid 
% yPressCI95+meanPressBootstrap


