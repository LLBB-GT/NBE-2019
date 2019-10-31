clear
clc

%%
%Pumping Pressure Measurement
%Pressure and intensity data will be imported from text files (tab delimited) 
%and will be analyzed to determine the measured pumping pressure.
%Browse for pressure and intensity files
[pFileName, pPathName] = uigetfile('*.txt','Select the pressure data');
[iFileName, iPathName] = uigetfile('*.txt','Select the intensity data');
%Full file names for the pressure and intensity files
pressureFile = strcat(pPathName,pFileName);
intensityFile = strcat(iPathName,iFileName);
%Input frame where pumping pressure measurement started to adjust time data
%for pressure measurement
prompt = 'Which frame did the pumping pressure measurement start on?';
dlg_title = 'Input Frame';
num_lines = 1;
PressureFrameStart = inputdlg(prompt,dlg_title,num_lines,{'0'});
PressureFrameStartNum = str2double(PressureFrameStart{1});
%Read in data from pressure file, and separate out adjusted time data and
%the pressure data
PressureWksht = tdfread(pressureFile);
PressureFieldNames = fieldnames(PressureWksht);
PressureTimeData = PressureWksht.(PressureFieldNames{1})+PressureFrameStartNum;
PressureData = PressureWksht.(PressureFieldNames{2});
%Read in data from intensity file, and separate out time data and the
%intensity data
IntensityWksht = tdfread(intensityFile);
IntensityFieldNames = fieldnames(IntensityWksht);
% IntensityTimeData = IntensityWksht.(IntensityFieldNames{1});
% IntensityData = IntensityWksht.(IntensityFieldNames{2});
IntensityData = IntensityWksht.(IntensityFieldNames{1});
IntensityTimeData=zeros(1,length(IntensityData))';
for k=2:length(IntensityData)
    IntensityTimeData(k) = IntensityTimeData(k-1)+1.15;
end
%Filter Intensity Data using a moving average filter (averaging over 15
%points)
timeIntervalFilter = 15;
FilterMat = ones(1,timeIntervalFilter)/timeIntervalFilter;
IntensityDataFiltered = filter(FilterMat, 1, IntensityData);
% Code to view filtered data vs actual data:
% figure
% hold on
% plot(IntensityTimeData,IntensityData)
% plot(IntensityTimeData,IntensityDataFiltered)
% legend('Intensity Data','Intensity Data Filtered')

%%
%Determining Pumping Pressure Value
count = 0;
%Find time point where pressure first goes above 80 mmHg (t1) and when it
%decreases below 60 mmHg after holding the pressure at 80 mmHg (t2)
for i=1:length(PressureData)
    if PressureData(i)>80 && count==0
        t1_P = floor(PressureTimeData(i));
        [Val1, t1] = min(abs(IntensityTimeData-t1_P)); 
        count = 1;
    elseif PressureData(i)<60 && count==1
        t2_P = floor(PressureTimeData(i));
        [Val2, t2] = min(abs(IntensityTimeData-t2_P));
        count = 2;
    elseif PressureData(i)<1 && count==2
        t3_P = floor(PressureTimeData(i));
        [Val3, t3] = min(abs(IntensityTimeData-t3_P));
        break;
    end
end
%Find min intensity after t1 and find index of value in
%IntensityDataFiltered
[minIntensity, minIntIndex_NA] = min(IntensityDataFiltered(t2:t3));
minIntIndex = length(IntensityDataFiltered) - length(IntensityDataFiltered(t2:end)) + minIntIndex_NA;
%Find max intensity after t2 and find index of value in
%IntensityDataFiltered
[maxIntensity, maxIntIndex_NA] = max(IntensityDataFiltered(t2:t3));
maxIntIndex = length(IntensityDataFiltered) - length(IntensityDataFiltered(t2:end)) + maxIntIndex_NA;
%Average min and max intensity to find the mid intensity
midIntensity = mean([minIntensity maxIntensity]);
%Find the index and time where the intensity data is closest to the mid
%intensity
[midIntVal, midIntIndex_NA] = min(abs(IntensityDataFiltered(minIntIndex:maxIntIndex)-midIntensity));
midIntIndex = length(IntensityDataFiltered) - length(IntensityDataFiltered(minIntIndex:end)) + midIntIndex_NA;
midIntTime = IntensityTimeData(midIntIndex);
%Find the index and pumping pressure where the intensity data is closest to
%the mid intensity
[pressPumpVal, pressPumpIndex] = min(abs(PressureTimeData-midIntTime));
PumpingPressure = PressureData(pressPumpIndex);
uiwait(msgbox(strcat('The Pumping Pressure is:',num2str(PumpingPressure),' mmHg')));

%%
%Emptying Rate
%Find how quickly the intensity decreases to its minimum value
[maxIntEmpRate, maxIntEmpRateIndex_NA] = max(IntensityDataFiltered(t1:t2));
maxIntEmpRateIndex = length(IntensityDataFiltered) - length(IntensityDataFiltered(t1:end)) + maxIntEmpRateIndex_NA;
minThreshold = 1.10*minIntensity;
for j=maxIntEmpRateIndex:t3
    if IntensityDataFiltered(j)<minThreshold
        threshIntIndex = j;
        break;
    end
end
EmpRate = round((maxIntEmpRate-minThreshold)./(IntensityTimeData(threshIntIndex)-IntensityTimeData(maxIntEmpRateIndex)),1);
uiwait(msgbox(strcat('The Emptying Rate is:',num2str(EmpRate),' IU/s')));

%%
%Plotting Data
figure
hold on
title('Pumping Pressure Measurement')
yyaxis right
plot(PressureTimeData,PressureData)
axis([0 max(max(PressureTimeData),max(IntensityTimeData)) -5 90])
xlabel('Time (s)')
ylabel('Pumping Pressure (mmHg)')
yyaxis left
plot(IntensityTimeData,IntensityDataFiltered)
%Line to denote the pumping pressure on the graph
line([PressureTimeData(pressPumpIndex) PressureTimeData(pressPumpIndex)],[min(IntensityDataFiltered)-500 max(IntensityDataFiltered)+500],'Marker','.','LineStyle','-','Color','k')
%Line for emptying rate
line([IntensityTimeData(maxIntEmpRateIndex) IntensityTimeData(threshIntIndex)],[maxIntEmpRate minThreshold],'Marker','.','LineStyle','-','Color','r')
axis([0 max(max(PressureTimeData),max(IntensityTimeData)) min(IntensityDataFiltered)-500 max(IntensityDataFiltered)+500])
ylabel('Average Intensity')