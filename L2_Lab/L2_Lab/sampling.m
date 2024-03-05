clear; clc; close all;
% This is a class exercise of Modern Signal Analysis and Data Processing.
% Sampling Theorem
% Haoran Meng, Fengjiang Ju & Wenyue Xia, Feb 18, 2023

Fs0 = 32; % Reference sample frequency.
t0 = 0:1/Fs0:20;
y0 = cos(2*pi*t0); % The original function for sampling


%% Sampling example 1
figure
for i = 1:4
    subplot(4,1,i)
    Fs = (i+2)/3;% sample frequency
    t = 0:1/Fs:20;% the sample time 
    y = cos(2*pi*t);% the values after sampling
    plot(t0,y0,'-r','linewidth',2); % The original function
    hold on;
    scatter(t,y,50,'filled','MarkerFaceColor','k'); % The points after sampling
    plot(t,y,'-k','linewidth',2); % The function after sampling(continuous)
    ylim([-1 1]); xlim([0 10])
    ylabel('Amplitude');
    title(['cos(2\pif_0t), f_0 = 1 Hz, f_s = ' num2str(Fs) ' Hz'])
    set(gca,'fontsize',20)
end
xlabel('Time (sec)')

%% Sampling example 2: same to example 1
figure
scale = 3;% Scales for changeing sample frequency
for i = 1:4
    subplot(4,1,i)
    Fs = Fs0/scale^(i-1);
    t = 0:1/Fs:20;
    y = cos(2*pi*t);
    plot(t0,y0,'-r','linewidth',2); hold on;
    scatter(t,y,50,'filled','MarkerFaceColor','k'); 
    plot(t,y,'-k','linewidth',2); 
    ylim([-1 1]); xlim([0 10])
    ylabel('Amplitude');
    title(['cos(2\pif_0t), f_0 = 1 Hz, f_s = 32/' num2str(scale^(i-1)) ' Hz'])
    set(gca,'fontsize',20)
end
xlabel('Time (sec)')