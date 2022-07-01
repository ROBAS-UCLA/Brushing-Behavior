function [timeStamps, accModifiedScaled, gyrModifiedScaled, pressureModified] = pathToDataNewNewest(path_part, DayNum, id, varargin)
%if isempty(find(dirpath_part.name,  num2str(Days_labeling{ids(f)}(g))))

if contains(path_part, 'AT-HOME')
    NewOrNot = 0;
else
    NewOrNot = NewDays(DayNum, id);
end

if NewOrNot
    DayorSession = '/Session';
    path_day = strcat(path_part, DayorSession, {' '}, num2str(DayNum), {' '}, 'New');
    path_day = path_day{1};

    if ~exist(path_day, 'dir')
        DayorSession = '/Day';
        path_day = strcat(path_part, DayorSession, {' '}, num2str(DayNum), {' '}, 'New');
        path_day = path_day{1};
    end

else
    DayorSession = '/Session';
    path_day = strcat(path_part, DayorSession, {' '}, num2str(DayNum));
    path_day = path_day{1};

    if ~exist(path_day, 'dir')
        DayorSession = '/Day';
        path_day = strcat(path_part, DayorSession, {' '}, num2str(DayNum));
        path_day = path_day{1};
    end          

end

addpath(genpath(path_day));

% exist(path_day, 'dir')
% [rawAccBrush, rawQuaBrush, rawMagBrush, rawAccHead, rawQuaHead, rawMagHead] ...
%         = DaytorawUpdated(DayNum, id);

dTime = desiredTimes(id, DayNum);
[time, acc, gyr, pressure] = boxToData(path_day, dTime);
       
timeStamps=min(time):1000/25:max(time);


if length(timeStamps) > 10^5
    error('TimeStamps array is too big');
end
%% Data Interpolation
[uniqueTime, uniqueIndices] = unique(time);
accUniqueTime = acc(uniqueIndices,:);
gyrUniqueTime = gyr(uniqueIndices,:);
pressureUniqueTime = pressure(uniqueIndices);

accModified = interp1(uniqueTime,accUniqueTime,timeStamps,'spline');
gyrModified = interp1(uniqueTime,gyrUniqueTime,timeStamps,'spline');
pressureModified = [interp1(uniqueTime,pressureUniqueTime,timeStamps,'previous')]';


%% Axes polarity
polarityAcc = [1, -1, 1];
polarityGyr = [-1, 1, -1];

accModifiedPol = accModified./repmat(polarityAcc, size(accModified,1), 1);
gyrModifiedPol = gyrModified./repmat(polarityGyr, size(gyrModified,1), 1);

% figure;plot(accModifiedPol);
% xlabel('sample');
% ylabel('g');
% title('Accelerometer OralB Sign Corrected');
% legend('X', 'Y', 'Z');
% 
% figure;plot(gyrModifiedPol);
% xlabel('sample');
% ylabel('rad/s');
% title('Gyroscope OralB Sign Corrected');
% legend('X', 'Y', 'Z');
        

%% Axes scales
scalesAcc = 32.17405*[1, 1, 1];
scalesGyr = 0.16*[1, 1, 1];


accModifiedScaled = accModifiedPol./repmat(scalesAcc, size(accModifiedPol,1), 1);
gyrModifiedScaled = gyrModifiedPol./repmat(scalesGyr, size(gyrModifiedPol,1), 1);


%% Scale 


end
