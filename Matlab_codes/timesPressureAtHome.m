function [startVec,endVec] = timesPressure(pressure)
% Objective: Determine the start and end times of when the pressure was 
%            held. Start time will be determined by any instance of 1. 
%            End time will be determined by the first 0 that follows the 
%            1 and does not have the same time stamp as the 1
% Input Variables: pressure = a structure 
% Output Variables: startVec = vector containing the start time of the 
%                              pressures
%                   endVec = vector containing the end times of the
%                            pressurs 
%Function called: None

pressureValues = pressure.Value;
pressureValues = [0;pressureValues];
pressureTimes = pressure.Time;

startVec = pressureTimes(find(diff(pressureValues) == 1));
endVec = pressureTimes(find(diff(pressureValues) == -1)+1);

end % end of fucntion durationPressure

