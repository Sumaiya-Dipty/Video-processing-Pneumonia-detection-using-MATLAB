clc; clear all; close all;

%% Input Video File
fprintf('Reading Video File...\n');
vid= VideoReader('E:\respiration files\pneumonia patients\patient-3\1.mp4');
fps= vid.Framerate;                 %Frames per second
Fs= 30;                    %sample per second
int= fps/Fs;              %sampling interval
nof= 600;  %total number of frames


%% Extracting RGB Values
i=1;                      %for storing RGB values
j=0;                      %for counting actual time in seconds
I = cell(1,3);
rects = zeros(3,4);
time=zeros(i,3);

r=zeros(i,3);
g=zeros(i,3);
b=zeros(i,3);

meanR = zeros(i,3);
meanG = zeros(i,3);
meanB = zeros(i,3);
meanRGB=zeros(i,3);

fprintf('\nSelect 1 ROI for calculating Respiration Rate. \n');
fprintf('and again select 2 ROI for calculating phase difference.\n');

for count = 1:int:nof
    img = read(vid,count);
    
    if count == 1
        for c = 1:3
            [I{c},rects(c,:)] = imcrop(img);
        end
    else
        for c = 1:3
            I{c} = imcrop(img,rects(c,:));
        end
    end
    
    for c = 1:3
        r= I{c}(:,:,1);
        g= I{c}(:,:,2);
        b= I{c}(:,:,3);
        
        meanR(i,c)= mean(mean(r));
        meanG(i,c)= mean(mean(g));
        meanB(i,c)= mean(mean(b));
        meanRGB(i,c)= (meanR(i,c)+meanG(i,c)+meanB(i,c))/3;
    end
    
    time(i) = j;
    j=j+(1/Fs);
    i=i+1;
end
fprintf('\nProcessing Data...\n');

%% Filtering RGB Values
Ww=meanRGB(:,1);                                  %RGB of Chest-abdominal region
Xx=meanRGB(:,2);                                  %chest
Yy=meanRGB(:,3);                                  %diaphragm
%Zz=meanRGB(:,4);                                  %abdomen

fl=0.833;                                           %lower cutoff 6 per minutes
fh=2.5;                                             %higher cutoff 120 per minutes
[b,a] = butter(3,[fl*2/Fs fh*2/Fs],'bandpass');   %butterworth bandpass filter

o = filtfilt(b,a,Ww);                              %filtered RGB (chest-abdominal) 
p = filtfilt(b,a,Xx);                              %filtered RGB (chest)
q = filtfilt(b,a,Yy);                              %filtered RGB (diaphragm)
%r = filtfilt(b,a,Zz); %filtered RGB (abdomen)

%% Respiration Rate calculation
L = length(o);
Y = fft(o,L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
[maxpeak,indexf]=max(P1);
maxvalf=f(indexf);
frequency = maxvalf;
rr=frequency*60;
fprintf('\nRespiration Rate:%.2f per minute.\n',rr);

%% Display Plots
figure(2)
subplot(121)
plot(time,Ww);  
title('RGB values of chest-abdominal region');
xlabel('time (seconds)');
ylabel('Amplitude');

subplot(122)
plot(time,o);  
title('Filtered RGB values of chest-abdominal region');
xlabel('time (seconds)');
ylabel('Amplitude');

figure(3)
subplot(2,3,1)
plot(time,p);
title('Filtered RGB values of chest region')
xlabel('time (seconds)');
ylabel('Amplitude');

subplot(2,3,2)
plot(time,q);
title('Filtered RGB values of diaphragm region')
xlabel('time (seconds)');
ylabel('Amplitude');

% subplot(2,3,3)
% plot(time,r);
% title('Filtered RGB values of abdomin region')
% xlabel('time (seconds)');
% ylabel('Amplitude');

subplot(2,3,[4,5])
plot(time,p,time,q);
xlabel('time (seconds)');
ylabel('Amplitude');


% phase difference measurement
PhDiff = phdiffmeasure(p, q);
PhDiff = rad2deg(PhDiff);
% display the phase difference
disp(['Phase difference Y->X = ' num2str(PhDiff) ' deg'])
commandwindow