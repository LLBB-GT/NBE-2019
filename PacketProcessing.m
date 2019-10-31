clear
clc
threshold = 0.001; %change

%Packet Processing
%Intensity data will be imported from text files (tab delimited) 
%and will be analyzed to determine the functional metrics.
%Browse for intensity file
[FileName, PathName] = uigetfile('*.txt','Select the intensity data');
%Full file names for the intensity file
intensityFile = strcat(PathName,FileName);

% filepath = '\\130.207.40.141\public\Matthew Cribb\Single Vessel Ligation Experiments\10 Weeks Post Surgery - 2.21.17\Mouse 8\Mouse8_RightFunction\';
% filename = 'Right Vessel Average Intensity';
% extension = '.txt';
% 
% fileloc= strcat(filepath, filename, extension);
x= dlmread(intensityFile, '\t');

%%x2=Wounded3(9573:9921,1);

y2=smooth(x,3);
y2=y2(200:end);

[maxtab2, min2]=peakdet(y2,150); %Intensity change to be considered a peak


if maxtab2(1,1) < min2(1,1)
    maxtab2 = maxtab2(2:end,:);
end

if maxtab2(end,1) > min2(end,1)
    maxtab2 = maxtab2(1:end-1,:);
end

mintab2=min2(1,:);
for n=2:length(min2)-1
    index=min2(n);
    while y2(index) <= min2(n,2)*(1+threshold)
        index=index-1;
    end
    mintemp2(1,1)=index;mintemp2(1,2)=y2(index);
    index=min2(n);
    while y2(index) <= min2(n,2)*(1+threshold)
        index=index+1;
    end
    mintemp2(2,1)=index;mintemp2(2,2)=y2(index);
    mintab2 = [mintab2 ; mintemp2];
end
mintab2 = [mintab2 ; min2(end,:)];

packet_width2=[];
mean_packet_min2=[];
packet_integral2=zeros(1,length(maxtab2));
packet_boundary2=[];
for n=1:2:length(mintab2)-1
    packet_width2 = [packet_width2 mintab2(n+1,1)-mintab2(n,1)];
    mean_packet_min2 = [mean_packet_min2 mean([mintab2(n+1,2) mintab2(n,2)])];
end

integral_y2=[];
for n=1:length(maxtab2)
    packet_line_x2 = mintab2(n*2-1,1):mintab2(n*2,1);
    packet_slope2 = (mintab2(n*2,2)-mintab2(n*2-1,2))/(mintab2(n*2,1)-mintab2(n*2-1,1));
    packet_offset2 = mintab2(n*2-1,2) - packet_slope2 * mintab2(n*2-1,1);
    for index=1:length(packet_line_x2)
        packet_line_y2(index) = packet_line_x2(index)*packet_slope2 + packet_offset2;
    end
    packet_line2 = [packet_line_x2; packet_line_y2];
    packet_line2 = packet_line2';
    packet_boundary2 = [packet_boundary2;packet_line2];

    packet_amplitude2(n) = maxtab2(n,2) - mean_packet_min2(n);
    packet_amplitude_perdiff2(n) = packet_amplitude2(n)/mean_packet_min2(n);
for index = mintab2(n*2-1,1):mintab2(n*2,1)
        integral_y2=[integral_y2;y2(index)];
    end
    for index = 1:length(integral_y2)
    packet_integral2(n) = packet_integral2(n) + integral_y2(index)-packet_line_y2(index);
    end
    packet_integral_norm2(n) = packet_integral2(n)/mean_packet_min2(n);
    packet_line_y2=[];
    integral_y2=[];
end

fig2 = figure;
plot(y2, 'linewidth', 2)
hold on
plot(maxtab2(:,1),maxtab2(:,2),'r.',mintab2(:,1),mintab2(:,2),'g.', 'markersize', 15)
plot(packet_boundary2(:,1),packet_boundary2(:,2),'g.','markersize', 4)
xlabel('Frame');ylabel('Intensity');
%title([filename ' Plot 2']);


%print(fig2, '-dtiff', [filepath filename 'Plot2'])


%packets per min, assuming 10fps
packet_frequency2 = length(maxtab2)/(mintab2(end,1)-mintab2(2,1))*600;

avg_packet_width2 = mean(packet_width2);
avg_packet_amplitude2 = mean(packet_amplitude2);
avg_packet_amplitude_perdiff2 = mean(packet_amplitude_perdiff2);
avg_packet_integral2 = mean(packet_integral2);
avg_packet_integral_norm2 = mean(packet_integral_norm2);

%normalized per minute, assuming 10fps
packet_transport2 = sum(packet_integral2)/(mintab2(end,1)-mintab2(1,1))*600;
packet_transport_norm2 = sum(packet_integral_norm2)/(mintab2(end,1)-mintab2(1,1))*600;

output = [packet_frequency2; avg_packet_width2;avg_packet_amplitude2;avg_packet_integral2;packet_transport2];

%filename

packet_frequency2 

avg_packet_width2

avg_packet_amplitude2

avg_packet_amplitude_perdiff2

avg_packet_integral2

avg_packet_integral_norm2

packet_transport2

packet_transport_norm2

% dlmwrite([filepath filename 'Output' extension], output)