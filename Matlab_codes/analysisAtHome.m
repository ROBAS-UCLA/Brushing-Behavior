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
path_orig = '~/Desktop/research/Robas/Data/BrushingDataSamples/';% CHANGE 1 : add the path to at home study data

%% Intitialization

participant_names = {'Isabel', 'Sumukh', 'Apurva', 'Quit', 'Rachel', ...
    'Jorge', 'Dayanara', 'Mariana', 'Shivam', 'Nandini', 'Christian', ...
    'Veronica', 'June', 'Maria', 'David'};    

region_names_new = {'ManRO', 'ManRB', 'ManRL', 'ManAB', 'ManAL', ...
    'ManLO', 'ManLB', 'ManLL', 'MaxRO', 'MaxRB', 'MaxRL', 'MaxAB', ...
    'MaxAL', 'MaxLO', 'MaxLB', 'MaxLL'};
sampRate = 25;

ids = [1:3,5:14]; 

Days_labeled = {[23,24,26:32, 33], [19:28], [1:10], [], [1:3,5,6,8:11,15], [1:10], [15:24], [13:15, 17:22, 23], [22:31], [10,11,13:15,18:22], ...
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
 
% P2 day 21 the accuracy is 4% ! P3 day 1 the accuracy is 8%! P3 day 2 the accuracy is 10%! P3 day 3 the accuracy is 9%! P3 day 4 the accuracy is 16%! 
% P3 day 5 the accuracy is 12%! P3 day 6 the accuracy is 13%! P3 day 8 the accuracy is 11%! P3 day 9 the accuracy is 14%! P3 day 10 the accuracy is 9%!
% P6 day 7 the accuracy is 13%! P6 day 8 the accuracy is 9%! P6 day 10 the accuracy is 18%! P9 day 23 the accuracy is 2%! P9 day 24 the accuracy is 4%!
% P9 day 25 the accuracy is 12%! P9 day 26 the accuracy is 10%!  P9 day 27 the accuracy is 11%! P9 Day 30 the accuracy is 10%! P9 Day 31 the accuracy is 10%!
% P10 idk if i missed any low accuracy!
% P11 day 7 the accuracy is 4%! P11 day 8 the accuracy is 10%! P11 day 16 the accuracy is 8%!P11 day 17 the accuracy is 10%!
% P12 day  the accuracy is %!


%% Initialization
brushingSampleCountsAll = cell(length(ids),1);
pressureSampleCountsAll = cell(length(ids),1);
sessionDurationAll = cell(length(ids),1);
transitionMatrixAll = cell(length(ids),1);

subjectSumBrushingSampleCounts = zeros(length(region_names_new), length(ids));
subjectSumPressureSampleCounts = zeros(length(region_names_new), length(ids));
subjectSumTransitionMatrix = zeros(length(region_names_new), length(region_names_new), length(ids));


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
        sessionDurationAll{f} = [sessionDurationAll{f}; size(data,1)/sampRate];
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
        
%         doubleVars = cell(1, length(brushingSampleCounts));
%         doubleVars{:} = {'double'};
%         VarTypes = [{'categorical'}, {'categorical'}, doubleVars(:)];
        array  = [f, g, brushingSampleCounts'];
        tableRegionBrushingAllReg = [tableRegionBrushingAllReg; array2table(array, ...
             'VariableNames', [{'Participant'}, {'Session'}, region_names_new(:)'])]; %'VariableTypes', VarTypes
        
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
        
        %% Transition matrix
        transitionMatrix = zeros(length(region_names_new));
        for i = 2:length(labels)
            if labels(i) ~= labels(i-1)
                transitionMatrix(labels(i-1), labels(i)) = transitionMatrix(labels(i-1), labels(i)) + 1;
            end
        end
        transitionMatrixAll{f} = [transitionMatrixAll{f}, {transitionMatrix}];
%         figure;heatmap(region_names_new, region_names_new, transitionMatrix)
    end
    sumBrushingSampleCounts = sum(brushingSampleCountsAll{f},2);
    subjectSumBrushingSampleCounts(:,f) = sumBrushingSampleCounts;
    normalizedSumBrushingSampleCounts = sumBrushingSampleCounts'/sum(sumBrushingSampleCounts);
    nnz(isnan(normalizedSumBrushingSampleCounts))
    normalizedSumBrushingSampleCounts(isnan(normalizedSumBrushingSampleCounts)) = 0;
    
    figure;bar(categorical(region_names_new), normalizedSumBrushingSampleCounts);
    xlabel('Dental regions');
    ylabel('Mean proportion of brushing duration');
    title("Distribution (pmf) of time spent by "+ "P"+ num2str(f)+ " brushing each region");
    
    sumPressureSampleCounts = sum(pressureSampleCountsAll{f},2);
    subjectSumPressureSampleCounts(:,f) = sumPressureSampleCounts;
    normalizedSumPressureSampleCounts = sumPressureSampleCounts'/sum(sumPressureSampleCounts);
    normalizedSumPressureSampleCounts(isnan(normalizedSumPressureSampleCounts)) = 0;
    figure;bar(categorical(region_names_new), normalizedSumPressureSampleCounts);
    xlabel('Dental regions');
    ylabel('Mean proportion of brushing duration');
    title("Distribution (pmf) of time spent by "+ "P"+ num2str(f)+ " brushing each region with excessive pressure");
   
%     figure;histogram(pressureSampleCountsAll{f}(5,:), 10);
    sumTransitionMatrix=sum(cat(3,transitionMatrixAll{f}{:}),3);
    subjectSumTransitionMatrix(:,:,f) = sumTransitionMatrix;
    normalizedSumTransMat = sumTransitionMatrix./sum(sumTransitionMatrix,2);
    normalizedSumTransMat(isnan(normalizedSumTransMat))=0;
    
    
    figure;heatmap(region_names_new, region_names_new, normalizedSumTransMat)
    title(strcat('P', num2str(f), 'transition matrix'));
   
    prunedNormalizedSumTransMat = normalizedSumTransMat;
    prunedNormalizedSumTransMat(prunedNormalizedSumTransMat<0.3) = 0;
    
    G = digraph(prunedNormalizedSumTransMat, region_names_new);
    ModifiedG = rmnode(G,find(indegree(G) == 0 & outdegree(G) == 0));
    LWidths = 4*ModifiedG.Edges.Weight/max(ModifiedG.Edges.Weight);
    
    %%% figure 10
    figure;plot(ModifiedG,'Layout','force','EdgeLabel',ModifiedG.Edges.Weight,'LineWidth',LWidths, 'ArrowSize', 13, 'NodeFontSize', 9)
    title(strcat('P', num2str(f),' transition graph'));
    %%%
end


%% population level
%% Session Duration

%%%% Figure 4
AllsessionDuration = cat(1,sessionDurationAll{:});
AllparticipantsCat = categorical(Allparticipants);
figure;boxplot(AllsessionDuration, Allparticipants);
ylim([0,170])
xlabel('Participant number');
ylabel('Duration (in seconds)');
title('Brushing session duration (in seconds)');
% median(AllsessionDuration(1:9))
% mean(AllsessionDuration)
% std(AllsessionDuration)
%%%%%%%%

%% Save for Doug
brushingDurationsAllSess = table(AllparticipantsCat, AllsessionDuration);
writetable(brushingDurationsAllSess, 'Total duration of brushing session.csv')
writetable(tableRegionBrushingAllReg, 'Duration of brushing dental regions.csv')
%%

sessionDurationTable = table(AllsessionDuration, AllparticipantsCat, 'VariableNames', {'Duration', 'Participants'});
glmeSessDur = fitglme(sessionDurationTable,...
		'Duration ~ 1 + (1|Participants)',...
		'Distribution','Gamma','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaSessDur,betanamesSessDur] = fixedEffects(glmeSessDur);
[BSessDur,BNamesSessDur] = randomEffects(glmeSessDur);


%%%% Appendix section 3
sessionDurationSampTable = table(floor(AllsessionDuration*25), AllparticipantsCat, 'VariableNames', {'Duration', 'Participants'});
glmeSessDurSamp = fitglme(sessionDurationSampTable,...
		'Duration ~ 1 + (1|Participants)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaSessDurSamp,betanamesSessDurSamp] = fixedEffects(glmeSessDurSamp);
[BSessDurSamp,BNamesSessDurSamp, statsSessDurSamp] = randomEffects(glmeSessDurSamp);
exp(betaSessDurSamp)/25
%%%%%%

%% Brushing region duration

%%% Figure 6
for i = 1:length(ids)
    AllsessionsRegionDuration = reshape(brushingSampleCountsAll{i}, ...
        size(brushingSampleCountsAll{i},1)*size(brushingSampleCountsAll{i},2),1);
    AllregionsCat = repmat(categorical(region_names_new'),size(brushingSampleCountsAll{i},2),1);
    figure;boxplot(AllsessionsRegionDuration/sampRate, AllregionsCat);
    xlabel('Dental regions');
    ylabel('Duration (in seconds)');
    title("Time spent by " + "P" + num2str(i) + " on each dental region (in seconds)");
end
%%%%

%%%% Figure 5
populationSumBrushingSampleCounts = sum(subjectSumBrushingSampleCounts,2);
figure;bar(categorical(region_names_new), populationSumBrushingSampleCounts'/sum(populationSumBrushingSampleCounts));
xlabel('Dental regions');
ylabel('Mean proportion of brushing duration');
title('Population-level distribution (pmf) of time spent brushing on each region');
% alaki = populationSumBrushingSampleCounts'/sum(populationSumBrushingSampleCounts);
% sum(alaki(9:end))
% sum(alaki(1:8))
%%%%

FixedglmeBrushRegDur = fitglme(tableRegionBrushingEachReg,...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB ',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[FixedbetaBrush,FixedbetanamesBrush] = fixedEffects(FixedglmeBrushRegDur);
[FixedBBrushRegDur,FixedBNamesBrushRegDur] = randomEffects(FixedglmeBrushRegDur);

glmeBrushRegDur = fitglme(tableRegionBrushingEachReg,...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1|Participant) + (1|Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaBrush,betanamesBrush] = fixedEffects(glmeBrushRegDur);
[BBrushRegDur,BNamesBrushRegDur] = randomEffects(glmeBrushRegDur);

glmeBrushRegDurNewestOrd = fitglme(tableRegionBrushingEachReg([4:end],:),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1|Participant) + (1|Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaBrushNewestOrd,betanamesBrushNewestOrd] = fixedEffects(glmeBrushRegDurNewestOrd);
[BBrushRegDurNewestOrd,BNamesBrushRegDurNewestOrd] = randomEffects(glmeBrushRegDurNewestOrd);

%%% appendix section 1
glmeBrushRegDurInter = fitglme(tableRegionBrushingEachReg([2:end, 1],:),...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1|Participant) + (1|Participant:Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaBrushInter,betanamesBrushInter] = fixedEffects(glmeBrushRegDurInter);
[BBrushRegDurInter,BNamesBrushRegDurInter, statsBrushRegDurInter] = randomEffects(glmeBrushRegDurInter);
%%%%


NoCorrglmeBrushRegDur = fitglme(tableRegionBrushingEachReg,...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (-1 + MaxOrNot|Participant) + (-1 + Side|Participant) + (-1 + OLB|Participant) + (-1 + MaxOrNot|Session) + (-1 + Side|Session) + (-1 + OLB|Session) + (1|Participant) + (1|Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[NoCorrbetaBrush,NoCorrbetanamesBrush] = fixedEffects(NoCorrglmeBrushRegDur);
[NoCorrBBrushRegDur,NoCorrBNamesBrushRegDur] = randomEffects(NoCorrglmeBrushRegDur);


CorrglmeBrushRegDur = fitglme(tableRegionBrushingEachReg,...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (MaxOrNot|Participant) + (Side|Participant) +  (OLB|Participant) + (MaxOrNot|Session) + (Side|Session) +  (OLB|Session) + (1|Participant) + (1|Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[CorrbetaBrush,CorrbetanamesBrush] = fixedEffects(CorrglmeBrushRegDur);
[CorrBBrushRegDur,CorrBNamesBrushRegDur] = randomEffects(CorrglmeBrushRegDur);

%%% Figure A.1
figure;histogram(tableRegionBrushingAllReg.MaxRO, 20, 'Normalization', 'pdf')
xlabel('Duration (# of samples)'); ylabel('Probability');
title('Distribution of brushing duration spent on MaxRO (# of samples)');
%%%%%
%% Zero inflated Poisson regression Brush Region
% fit the model 

%%%%%%% Appendix section 1
MaxOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.MaxOrNot))==2;
AntOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.Side))==1;
LeftOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.Side))==2;
BuccalOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.OLB))==1;
LingualOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.OLB))==2;
OcclusalOrNotValuesBrush = grp2idx(categorical(tableRegionBrushingEachReg.OLB))==3;

resBrush = EMzeropoisson_mat([tableRegionBrushingEachReg.sampleCounts MaxOrNotValuesBrush ...
    AntOrNotValuesBrush LeftOrNotValuesBrush LingualOrNotValuesBrush OcclusalOrNotValuesBrush]);
fittedtauBrush = resBrush(end,1);
fittedinterceptBrush = resBrush(end,2);
fittedbetaBrush = resBrush(end,3:end);
EMzeropoisson_mat([tableRegionBrushingEachReg.sampleCounts MaxOrNotValuesBrush ...
    AntOrNotValuesBrush LeftOrNotValuesBrush BuccalOrNotValuesBrush LingualOrNotValuesBrush],0.5,'initialtau','stable','tol',0.001)




no_iter = 1000;
fittedtauBrushingAll = zeros(no_iter,1);
fittedinterceptBrushingAll = zeros(no_iter,1);
fittedbetaBrushingAll = zeros(no_iter,5);
for i = 1:no_iter
    randomPerm = randperm(length(tableRegionBrushingEachReg.sampleCounts));
    samplesPerm = tableRegionBrushingEachReg.sampleCounts;
    samplesPerm = samplesPerm(randomPerm);
    resBrush = EMzeropoisson_mat([samplesPerm MaxOrNotValuesBrush ...
        AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesBrush LingualOrNotValuesBrush]);
    fittedtauBrushingAll(i) = resBrush(end,1);
    fittedinterceptBrushingAll(i) = resBrush(end,2);
    fittedbetaBrushingAll(i,:) = resBrush(end,3:end);
end

nnz(fittedtauBrush > fittedtauBrushingAll)/size(fittedtauBrushingAll,1)

figure; histogram(fittedinterceptBrushingAll,100); 
xlabel('Values');
ylabel('Probability');
title('Distribution of intercept');
nnz(fittedinterceptBrush > fittedinterceptBrushingAll)/size(fittedinterceptBrushingAll,1)

for i = 1:size(fittedbetaBrushingAll,2)
%     figure; histogram(fittedbetaPressAll(:,i),100); 
    min([nnz(fittedbetaBrush(i) > ...
    fittedbetaBrushingAll(:,i))/size(fittedbetaBrushingAll,1), nnz(fittedbetaBrush(i) < ...
    fittedbetaBrushingAll(:,i))/size(fittedbetaBrushingAll,1)]) 
%     xlabel('Values');
%     ylabel('Magnitude');
%     title('Distribution of coefficient');
end


%%%%%%%%

%% Pressure region duration
% for i = 1:length(ids)
%     AllsessionsRegionDuration = reshape(brushingSampleCountsAll{i}, ...
%         size(brushingSampleCountsAll{i},1)*size(brushingSampleCountsAll{i},2),1);
%     AllregionsCat = repmat(categorical(region_names_new'),size(brushingSampleCountsAll{i},2),1);
%     figure;boxplot(AllsessionsRegionDuration/sampRate, AllregionsCat);
%     xlabel('Region names');
%     ylabel('Duration (in seconds)');
%     title(strcat('P', num2str(i), ' Time spent on each dental region (in seconds)'));
% end


%%%%% Figure 8
for i = 1:length(ids)
    AllsessionsRegionPressureDuration = reshape(pressureSampleCountsAll{i}, ...
        size(pressureSampleCountsAll{i},1)*size(pressureSampleCountsAll{i},2),1);
    AllregionsCat = repmat(categorical(region_names_new'),size(pressureSampleCountsAll{i},2),1);
    figure;boxplot(AllsessionsRegionPressureDuration/sampRate, AllregionsCat);
    xlabel('Dental regions');
    ylabel('Duration (in seconds)');
    title("Time spent (in seconds) by " + "P" + num2str(i) + " brushing each region with excessive pressure");
end
%%%%%

%%%% Figure 7
populationSumPressureSampleCounts = sum(subjectSumPressureSampleCounts,2);
figure;bar(categorical(region_names_new), populationSumPressureSampleCounts'/sum(populationSumPressureSampleCounts));
xlabel('Dental regions');
ylabel('Probability');
title('Population-level distribution (pmf) of the time spent brushing each region with excessive pressure');
%%%%%


FixedglmePressRegDur = fitglme(tableRegionPressureEachReg,...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB ',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[FixedbetaPress,FixedbetanamesPress] = fixedEffects(FixedglmePressRegDur);
[FixedBPressRegDur,FixedBNamesPressRegDur] = randomEffects(FixedglmePressRegDur);

glmePressRegDur = fitglme(tableRegionPressureEachReg,...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1|Participant) + (1|Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaPress,betanamesPress] = fixedEffects(glmePressRegDur);
[BPressRegDur,BNamesPressRegDur] = randomEffects(glmePressRegDur);


%%%% Appendix section 2
glmePressRegDur = fitglme(tableRegionPressureEachReg,...
		'sampleCounts ~ 1 + MaxOrNot + Side + OLB + (1|Participant) + (1|Participant:Session)',...
		'Distribution','Poisson','FitMethod','Laplace');%, 'Link','log', 'DummyVarCoding','effects'
[betaPress,betanamesPress] = fixedEffects(glmePressRegDur);
[BPressRegDur,BNamesPressRegDur] = randomEffects(glmePressRegDur);
%%%%

figure;histogram(tableRegionPressureAllReg.ManRB, 20)

%% Zero inflated Poisson regression Pressure
% fit the model 
MaxOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.MaxOrNot))==2;
AntOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Side))==1;
LeftOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Side))==2;
BuccalOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.OLB))==1;
LingualOrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.OLB))==2;
P1OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==1;
P2OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==2;
P3OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==3;
P4OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==4;
P5OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==5;
P6OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==6;
P7OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==7;
P8OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==8;
P9OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==9;
P10OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==10;
P11OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==11;
P12OrNotValuesPress = grp2idx(categorical(tableRegionPressureEachReg.Participant))==12;

%%% Appendix section 2
resPress = EMzeropoisson_mat([tableRegionPressureEachReg.sampleCounts MaxOrNotValuesPress ...
    AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesPress LingualOrNotValuesPress]);
fittedtauPress = resPress(end,1);
fittedinterceptPress = resPress(end,2);
fittedbetaPress = resPress(end,3:end);
EMzeropoisson_mat([tableRegionPressureEachReg.sampleCounts MaxOrNotValuesPress ...
    AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesPress LingualOrNotValuesPress],0.5,'initialtau','stable','tol',0.001)

rng(0)
no_iter = 1000;
fittedtauPressAll = zeros(no_iter,1);
fittedinterceptPressAll = zeros(no_iter,1);
fittedbetaPressAll = zeros(no_iter,5);
for i = 1:no_iter
    randomPerm = randperm(length(tableRegionPressureEachReg.sampleCounts));
    samplesPerm = tableRegionPressureEachReg.sampleCounts;
    samplesPerm = samplesPerm(randomPerm);
    resPress = EMzeropoisson_mat([samplesPerm MaxOrNotValuesPress ...
        AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesPress LingualOrNotValuesPress]);
    fittedtauPressAll(i) = resPress(end,1);
    fittedinterceptPressAll(i) = resPress(end,2);
    fittedbetaPressAll(i,:) = resPress(end,3:end);
end

figure; histogram(fittedinterceptPressAll,100); 
xlabel('Values');
ylabel('Probability');
title('Distribution of intercept');
nnz(fittedinterceptPress > fittedinterceptPressAll)/size(fittedinterceptPressAll,1)

for i = 1:size(fittedbetaPressAll,2)
%     figure; histogram(fittedbetaPressAll(:,i),100); 
    min([nnz(fittedbetaPress(i) > ...
    fittedbetaPressAll(:,i))/size(fittedbetaPressAll,1), nnz(fittedbetaPress(i) < ...
    fittedbetaPressAll(:,i))/size(fittedbetaPressAll,1)]) 
%     xlabel('Values');
%     ylabel('Magnitude');
%     title('Distribution of coefficient');
end

rng(0)
no_bootstraps = 900;
fittedtauPressBSAll = zeros(no_bootstraps,1);
fittedinterceptPressBSAll = zeros(no_bootstraps,1);
fittedbetaPressBSAll = zeros(no_bootstraps,5);

    
no_participants = length(unique(tableRegionPressureEachReg.Participant));

for n = 1:no_bootstraps
    rng(n)

    bootsrapPressEachReg = table;

    for i = 1:no_participants
        randParticipant = randsample(no_participants, 1);        
           
        no_sessions = length(unique(tableRegionPressureEachReg.Session(tableRegionPressureEachReg.Participant==randParticipant)));
        
        for j = 1:no_sessions
            randSession = randsample(no_sessions, 1);     
            bootsrapPressEachReg = [bootsrapPressEachReg; tableRegionPressureEachReg( ...
                tableRegionPressureEachReg.Participant == randParticipant & tableRegionPressureEachReg.Session == randSession, :)];
            
        end
    end
    
    resPressBS = EMzeropoisson_mat([bootsrapPressEachReg.sampleCounts MaxOrNotValuesPress ...
        AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesPress LingualOrNotValuesPress]);
    fittedtauPressBSAll(n) = resPressBS(end,1);
    fittedbetaPressBSAll(n,:) = resPressBS(end,3:end);
    fittedinterceptPressBSAll(n) = resPressBS(end,2);
        
end

figure; histogram(fittedinterceptPressBSAll(1:100),100); 
xlabel('Values');
ylabel('Probability');
title('Distribution of intercept');

for i = 1:size(fittedbetaPressBSAll,2)
    figure; histogram(fittedbetaPressBSAll(1:900,i),100); 
    xlabel('Values');
    ylabel('Magnitude');
    title('Distribution of coefficient');
end

N = 100; %no_bootstrap;       
CI95 = tinv([0.025 0.975], N-1);                

% meanBootstrapInterc = mean(fittedinterceptPressBSAll(1:N,:),1);                             
% ySEMInterc = std(fittedinterceptPressBSAll(1:N,:))/sqrt(N);                            
% yCI95Interc = bsxfun(@times, ySEMInterc, CI95(:));            
% yCI95Interc+meanBootstrapInterc
quantile(fittedinterceptPressBSAll(1:N,:), 0.0275)
quantile(fittedinterceptPressBSAll(1:N,:), 0.975)


meanBootstrap = mean(fittedbetaPressBSAll(1:N,:),1);                             
ySEM = std(fittedbetaPressBSAll(1:N,:))/sqrt(N);                            
yCI95 = bsxfun(@times, ySEM, CI95(:));            
figure
plot(1:size(fittedbetaPressBSAll,2), meanBootstrap)                                   
hold on
plot(1:size(fittedbetaPressBSAll,2), yCI95+meanBootstrap)                       
hold off
grid 
yCI95+meanBootstrap



%%%%%


resPressDummySub = EMzeropoisson_mat([tableRegionPressureEachReg.sampleCounts MaxOrNotValuesPress ...
    AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesPress LingualOrNotValuesPress P1OrNotValuesPress ...
    P2OrNotValuesPress P3OrNotValuesPress P4OrNotValuesPress P5OrNotValuesPress P6OrNotValuesPress ...
    P7OrNotValuesPress P8OrNotValuesPress P9OrNotValuesPress P10OrNotValuesPress P11OrNotValuesPress ...
    P12OrNotValuesPress]);
fittedtauPressDummySub = resPressDummySub(end,1);
fittedinterceptPressDummySub = resPressDummySub(end,2);
fittedbetaPressDummySub = resPressDummySub(end,3:end);
EMzeropoisson_mat([tableRegionPressureEachReg.sampleCounts MaxOrNotValuesPress ...
    AntOrNotValuesPress LeftOrNotValuesPress BuccalOrNotValuesPress LingualOrNotValuesPress P1OrNotValuesPress ...
    P2OrNotValuesPress P3OrNotValuesPress P4OrNotValuesPress P5OrNotValuesPress P6OrNotValuesPress ...
    P7OrNotValuesPress P8OrNotValuesPress P9OrNotValuesPress P10OrNotValuesPress P11OrNotValuesPress ...
    P12OrNotValuesPress],0.5,'initialtau','stable','tol',0.001)

%% Transition Matrix and graph

%%%% Figure 9
populationSumTransitionMatrix = sum(subjectSumTransitionMatrix,3);
populationNormalizedSumTransMat = populationSumTransitionMatrix./sum(populationSumTransitionMatrix,2);
figure;heatmap(region_names_new, region_names_new, populationNormalizedSumTransMat)
title('Population level transition matrix');

populationPrunedNormalizedSumTransMat = populationNormalizedSumTransMat;
populationPrunedNormalizedSumTransMat(populationPrunedNormalizedSumTransMat<0.2) = 0;
    
G = digraph(populationPrunedNormalizedSumTransMat, region_names_new);
ModifiedG = rmnode(G,find(indegree(G) == 0 & outdegree(G) == 0));
LWidths = 4*ModifiedG.Edges.Weight/max(ModifiedG.Edges.Weight);
figure;plot(simplify(ModifiedG),'Layout','force','EdgeLabel',ModifiedG.Edges.Weight,'LineWidth',LWidths, 'ArrowSize', 13, 'NodeFontSize', 9)
title('Population level transition graph');
%%%%%



transitionMatrixSubject=transitionMatrixAll{5}{8};
normalizedSumTransMat = transitionMatrixSubject./sum(transitionMatrixSubject,2);
normalizedSumTransMat(isnan(normalizedSumTransMat))=0;
figure;heatmap(region_names_new, region_names_new, normalizedSumTransMat)
title(strcat('P', num2str(1), ' transition matrix'));

prunedNormalizedSumTransMat = normalizedSumTransMat;
prunedNormalizedSumTransMat(prunedNormalizedSumTransMat<0.4) = 0;

G = digraph(prunedNormalizedSumTransMat, region_names_new);
ModifiedG = rmnode(G,find(indegree(G) == 0 & outdegree(G) == 0));
LWidths = 4*ModifiedG.Edges.Weight/max(ModifiedG.Edges.Weight);
%%% Figure 10
figure;plot(ModifiedG,'Layout','force','EdgeLabel',ModifiedG.Edges.Weight,'LineWidth',LWidths, 'ArrowSize', 13, 'NodeFontSize', 9)
title(strcat('P', num2str(1),' transition matrix'));
%%%%

%% Brushing Scores     
BrushingRegsArrays = table2array(tableRegionBrushingAllReg(:,3:end));
PressureRegsArrays = table2array(tableRegionPressureAllReg(:,3:end));

TV = sum(abs(BrushingRegsArrays./sum(BrushingRegsArrays,2) - 1/16),2);

BrushWOPress = BrushingRegsArrays - PressureRegsArrays;

alpha = 5;
brushingScore = max(100*(sum(min(BrushWOPress./sampRate, 7.5),2))/120 - alpha*TV,0);
%%% Figure 11
figure;boxplot(brushingScore, Allparticipants);
ylim([0, 100])
xlabel('Participant number');
ylabel('Score');
title('Brushing score (out of 100)');
%%%%

mean(brushingScore)
std(brushingScore)


1;
