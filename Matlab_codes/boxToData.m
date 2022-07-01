function [epochSess, acc, gyr, pressure] = boxToData(path_day, varargin)

addpath(genpath(path_day));

datefolder_dir = dir(path_day);
 

for i = 1:numel(datefolder_dir)
    datefolder_dir_curr = datefolder_dir(i).name;
    if length(datefolder_dir(i).name) <= 5
        continue;
    end
    if strcmp(datefolder_dir_curr(1:2), '20') && strcmp(datefolder_dir_curr(end-2:end), 'zip')               
        unzip(datefolder_dir_curr, path_day);        
        data_dir_add = datefolder_dir_curr(1:end-4);         
        break;
    elseif strcmp(datefolder_dir_curr(1:2), '20')
        data_dir_add = datefolder_dir_curr;
        break;
    end    
end

datefoldername = strcat(path_day, '/', data_dir_add);
data_dir = dir(datefoldername);

sensorSelected = {'ORALB_BRUSH_MOTION', 'ORALB_PRESSURE'};

for j = 1:length(sensorSelected)
for i = 1:numel(data_dir)
    data_dir_curr = data_dir(i).name;
    if length(data_dir(i).name) <= 10
        continue;
    end
    if strcmp(data_dir_curr(1:length(sensorSelected{j})), sensorSelected{j}) && strcmp(data_dir_curr(end-2:end), 'zip')               
        unzip(data_dir_curr, path_day);        
        data_dir_add = data_dir_curr(1:end-4);         
        break;
    elseif strcmp(data_dir_curr(1:length(sensorSelected{j})), sensorSelected{j})
        data_dir_add = data_dir_curr;
        break;
    end    
end


motion_dir = dir(strcat(datefoldername, '/', data_dir_add));
% addpath(strcat(datefoldername, '/', data_dir_add))

textData = {};

% gunzip(strcat(path_day, '/', data_dir_add, '/*.gz'))

for i = 1:numel(motion_dir)
    if length(motion_dir(i).name) <= 10 || motion_dir(i).bytes == 0
        continue;
    end
    
    if strcmp(motion_dir(i).name(end-1:end), 'gz') 
        gzData = gunzip(motion_dir(i).name);
        fid = fopen(gzData{1},'r');
    else
        fid = fopen(motion_dir(i).name,'r');
    end
    textCurr = textscan(fid,'%s');
    textCurr = textCurr{1};
    textData = [textData;textCurr];
    fclose(fid);
end


    
% Allocate Memeory
sensData = []; %zeros(size(textData,1),8); % stores split string data

% Seperate text data based on data
for i = 1:length(textData) 
   sensData = [sensData; str2double(strsplit(textData{i,1},','))];    
end % end for loop

% Find epoch data in and adjust the order so epochVec goes from max to min
epochVec =  sensData(:,1);   % Stores epoch time
[epochVecSorted, sortIn] = sort(epochVec, 'ascend');
sensDataSorted = sensData(sortIn, :);

% Convvert epochVec from epoch to date time and store in timeVec
% for i = 1:length(epochVec)
%   timeVec(i)= datetime(epochVec(i)/1000, 'convertfrom','posixtime', 'TimeZone', 'America/Los_Angeles');
% end % end for loop


% Find if there any gaps in time and store in skippedTime array
if strcmp(sensorSelected{j}, 'ORALB_BRUSH_MOTION')
timeTh = 1000;
markStartGap = zeros(size(epochVecSorted));
for i = 1:length(epochVecSorted)-1
    if abs(epochVecSorted(i+1)-epochVecSorted(i)) > timeTh
        markStartGap(i) = 1;
       % fprintf("Hello");
    end
end 

startGapIndices = find(markStartGap == 1);
endGapIndices = [1;startGapIndices+1];

%% Just for plotting and debugging
timeVec = datetime(epochVecSorted/1000, 'TimeZone', 'America/Chicago', 'convertfrom','posixtime');


% % plot the sensor data
% figure 
% scatter(timeVec,accX);
% % 
% % Label graph
% title('Acc X-Axis Data')
% xlabel('Time')
% ylabel('AccX')
% %$legend('x-axis');
gapsStartTime = timeVec(startGapIndices);
gapsEndTime = timeVec(endGapIndices);

%%
if max(diff(epochVecSorted)) > timeTh && nargin < 2
    error('give a start time');
end

disp(strcat(num2str(length(endGapIndices)),' different gaps in the data'));

if nargin > 1 
    if length(varargin{1}) > 8
        startTimeFull = datetime(varargin{1},'InputFormat','yyyy-MM-dd''T''HH:mm:ss', 'TimeZone', 'America/Chicago');

    else
        startTime = datetime(varargin{1},'InputFormat','HH:mm:ss', 'TimeZone', 'America/Chicago');
        startDate = timeVec(1);
        startDate.Hour = 0;
        startDate.Minute = 0;
        startDate.Second = 0;
        startTimeFull = startDate + timeofday(startTime);
    end
    
    
    
%     startTimeFormatted = datetime(startTimeFull,'InputFormat','yyyy-MM-dd HH:mm:ss', 'TimeZone', 'America/Los_Angeles');
    startTimeStampUnix = posixtime(startTimeFull);
    [extra,closeInd] = min(abs(startTimeStampUnix-epochVecSorted/1000));
    [extra2, minInd] = min(abs(endGapIndices - closeInd));   
    indStartSess = endGapIndices(minInd);
else
    indStartSess = 1;
end

if minInd == length(endGapIndices)
    indEndSess = size(sensDataSorted,1);
else    
    indEndSess = startGapIndices(minInd);
end


indices = indStartSess:indEndSess;

epochSess = epochVecSorted(indices);

time = timeVec(indices);

end

if strcmp(sensorSelected{j}, 'ORALB_BRUSH_MOTION')

    % Store the sensor data
    epochSessMotion = epochSess;
    accX =  sensDataSorted(indices,3);      % stores sensor data in x axis
    accY =  sensDataSorted(indices,4);      % stores sensor data in y axis 
    accZ =  sensDataSorted(indices,5);      % stores sensor data in z axis
    gyrX =  sensDataSorted(indices,6);      % stores sensor data in x axis
    gyrY =  sensDataSorted(indices,7);      % stores sensor data in y axis 
    gyrZ =  sensDataSorted(indices,8);      % stores sensor data in z axis
    acc = [accX,accY,accZ];
    gyr = [gyrX,gyrY,gyrZ];
    
elseif strcmp(sensorSelected{j}, 'ORALB_PRESSURE')
    
    pressureDataSorted =  [];
    epochVecSortedPressure = [];
    for l = 1:length(epochVecSorted)
        if epochVecSorted(l) >= epochSessMotion(1) && epochVecSorted(l) <= epochSessMotion(end)
            epochVecSortedPressure = [epochVecSortedPressure;epochVecSorted(l)];
            pressureDataSorted = [pressureDataSorted;sensDataSorted(l,3)];               
        end
    end
    
    pressureStruct = struct('Time',epochVecSortedPressure,...
                  'Value',pressureDataSorted);
              
%     timeVecPressEpoch = datetime(epochVecSortedPressure/1000, 'TimeZone', 'America/Chicago', 'convertfrom','posixtime');
%     figure; stairs(timeVecPressEpoch, pressureDataSorted)

    [Start_Time, End_Time] = timesPressureAtHome(pressureStruct);
%     Elapsed_Time = durationPressure(Start_Time, End_Time);

%     epochSessPressure = epochSess;
    starttimeVecPressEpoch = datetime(Start_Time/1000, 'TimeZone', 'America/Chicago', 'convertfrom','posixtime');
    endtimeVecPressEpoch = datetime(End_Time/1000, 'TimeZone', 'America/Chicago', 'convertfrom','posixtime');
        
    pressure = zeros(length(epochSessMotion),1);
    for k = 1:length(Start_Time)
        [extra2, startInd] = min(abs(Start_Time(k) - epochSessMotion));
        [extra3, endInd] = min(abs(End_Time(k) - epochSessMotion));
        
        pressure(startInd:endInd) = 1;
    end
end



%% Pressure



end

end
