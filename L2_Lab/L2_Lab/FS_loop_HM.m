clear; clc; clf; close all;
% This is a class exercise of Modern Signal Analysis and Data Processing.
% Fourier Series
% Haoran Meng, Fengjiang Ju & Wenyue Xia, Feb 20, 2023

%% Parameters
T = 2; % Period
f = 1/T; % Frequency
w0 = 2*pi*f; % Angular frequency
dt = 0.01; % Sample rate
t = -3:dt:3; % Length of signal for illustrating 
tao = -1:dt:1; % Length of the minimum period 

%% Functions to decompose.
ftao = zeros(size(tao)); flag = 0;

% delta function
ftao(floor(length(ftao)/2)+1) = 1;

% step function
% ftao(floor(length(ftao)/2)+1:end) = 1;
% ftao = ftao - mean(ftao); flag = 1;

% ramp function
%ftao(floor(length(ftao)/2)+1:end) = tao(floor(length(ftao)/2)+1:end);
%ftao = ftao - mean(ftao);

% sinsoid function 1
%ftao = sin(pi*tao); flag = 1;

% exponetial
%ftao = exp(-tao);

% expoential 2
%ftao(1:floor(length(ftao)/2)+1) = exp(tao(1:floor(length(ftao)/2)+1));  
%ftao(floor(length(ftao)/2)+1:end) = exp(-tao(floor(length(ftao)/2)+1:end));  
%ftao = ftao - mean(ftao); flag = 1;

% triangle
%ftao(1:floor(length(ftao)/2)+1) = tao(1:floor(length(ftao)/2)+1) + 0.5;
%ftao(floor(length(ftao)/2)+1:end) = -tao(floor(length(ftao)/2)+1:end) + 0.5;
%flag = 1;

% Parabola
%ftao = -tao.^2; ftao = ftao + 1/3; flag = 1;

%% Fourier Series
N = 10; % N-th Fourier Series

a0 = trapz(tao,ftao)/T;% trapz means trapezoidal integral
f = a0; f1 = a0;
an = zeros(1,N);
bn = zeros(1,N);
amp = zeros(1,N);% Amplitude
ph = zeros(1,N);% Phase
% Fourier Series
figure(1)
for n = 1:N
    subplot(3,4,[5 9])
    % Plot the cosine function term
    fcos = 1.*cos(n*w0*tao); 
    an(n)=trapz(tao,ftao.*fcos)*2/T;% Decomposition coefficient
    plot(tao,fcos/1.2-n*2,'-k','linewidth',1.5); hold on;
    xlabel('t (s)'); title('cos(n\omegat)'); ylim([-(n+1)*2 0])
    set(gca,'fontsize',15,'YTickLabel',[])

    subplot(3,4,[6 10])
    % Plot the coefficient of cosine function term
    scatter(an(n),-n*2,50,'filled','MarkerFaceColor','k'); hold on;
    xlim([-1 1]); xlabel('a_n'); ylim([-(n+1)*2 0]); box on;
    title('Coefficients')
    set(gca,'fontsize',15,'YTickLabel',[]);

    subplot(3,4,[7 11])
    % Plot the cosine function term
    fsin = 1.*sin(n*w0*tao) ; bn(n)=trapz(tao,ftao.*fsin)*2/T;
    plot(tao,fsin/1.2-n*2,'-k','linewidth',1.5); hold on;
    xlabel('t (s)'); title('sin(n\omegat)'); ylim([-(n+1)*2 0])
    set(gca,'fontsize',15,'YTickLabel',[])

    subplot(3,4,[8 12])
    % Plot the coefficient of sine function term
    scatter(bn(n),-n*2,50,'filled','MarkerFaceColor','k'); hold on;
    xlim([-1 1]); xlabel('b_n'); ylim([-(n+1)*2 0]); box on;
    title('Coefficients')
    set(gca,'fontsize',15,'YTickLabel',[])

    f = f + an(n)*cos(n*w0*t)+bn(n)*sin(n*w0*t); % The sum of N series. 

    % Amplitude and phase graph of Fourier series
    amp(n) = sqrt(an(n)^2+bn(n)^2);
    ph(n) =  atan2(-bn(n),an(n));
    if abs(bn(n))<1e-8 && abs(an(n))<1e-8
        ph(n) = 0;
    end

    if ph(n)<-pi
        ph(n)=ph(n)+2*pi;
    elseif ph(n)>pi
        ph(n)=ph(n)-2*pi;
    end

    f1 = f1 + amp(n)*cos(n*w0*t + ph(n));
end

subplot(3,4,[1 2 3 4])
% Plot the sum of N series.
plot([tao-T tao tao+T],[ftao ftao ftao],'r','linewidth',3); hold on;
plot(t,f,'k','linewidth',1.5); hold on; 
legend('Analytical','Numerical')
xlabel('t (s)'); ylabel('Amplitude'); ylim([min(ftao)-0.1 max(ftao)+0.1]); grid on;
title(['Fourier Series with N = ' num2str(N) ]);
set(gca,'fontsize',15)

%% Plot the Ampliude-Phase diagram
omega = (1:N)*w0;
figure(2)

subplot(311)
plot([tao-T tao tao+T],[ftao ftao ftao],'r','linewidth',3); hold on;
plot(t,f1,'k','linewidth',1.5); hold on; 
legend('Analytical','Numerical')
xlim([-3 3]); xlabel('t(s)'); ylim([min(ftao)-0.1 max(ftao)+0.1]);  ylabel('Ampitude'); box on;
title(['Fourier Series with N = ' num2str(N) ]);
set(gca,'fontsize',15)
subplot(312)
plot(omega,amp,'k','linewidth',1.5);hold on;
scatter(omega,amp,50,'filled','MarkerFaceColor','k');
xlabel('\omega(rad/s)');ylabel('Ampitude');
ylim([-0.1 max(amp)+0.1]);
title('A_n=sqrt(a_n^2+b_n^2)');
set(gca,'fontsize',15)
subplot(313)
plot(omega,ph,'k','linewidth',1.5);hold on;
scatter(omega,ph,50,'filled','MarkerFaceColor','k');
xlabel('\omega (rad/s)');ylabel('Phase (rad)');
title('\Phi=arctan(-b_n/a_n)');
ylim([-pi*1.1 pi*1.1])
yticks([-pi 0 pi]); yticklabels({'\pi','0','\pi'})
set(gca,'fontsize',15)


%% compute derivatives and intergals
if flag == 1
    %% Obtain the derivative of f(t)
    ftaodt = diff(ftao)/mean(diff(tao)); ftaodt(end+1) = ftaodt(end);
    fdt = trapz(tao,ftaodt)/T;% trapz means trapezoidal integral;
    fdt1 = trapz(tao,ftaodt)/T;
    for n=1:N

        fcos = 1.*cos(n*w0*tao);
        an_d(n)=trapz(tao,ftaodt.*fcos)*2/T;% Decomposition coefficient
        fsin = 1.*sin(n*w0*tao) ;
        bn_d(n)=trapz(tao,ftaodt.*fsin)*2/T;
        fdt = fdt + an_d(n)*cos(n*w0*t)+bn_d(n)*sin(n*w0*t); % The sum of N series.

        amp2_d(n) = sqrt(an_d(n)^2+bn_d(n)^2);
        ph2_d(n) =  atan2(-bn_d(n),an_d(n));


        amp_d(n) = amp(n)*n*w0;
        ph_d(n) = ph(n)+pi/2;

        if ph_d(n)<-pi
            ph_d(n)=ph_d(n)+2*pi;
        elseif ph_d(n)>pi
            ph_d(n)=ph_d(n)-2*pi;
        end

        if abs(bn(n))<1e-8 && abs(an(n))<1e-8
            ph2_d(n) = ph_d(n);
        end



        fdt1 = fdt1 + amp_d(n)*cos(n*w0*t + ph_d(n));

    end

    % Plot the Ampliude-Phase diagram of Derivative
    figure(3)
    subplot(311)
    plot([tao-T tao tao+T],[ftaodt ftaodt ftaodt],'r','linewidth',3); hold on;
    plot(t,fdt1,'k','linewidth',1.5); hold on;
    xlim([-3 3]); xlabel('t(s)'); ylabel('Ampitude'); box on;
    legend('Analytical','Numerical')
    title(['Derivative of Fourier Series with N = ' num2str(N) ]);
    set(gca,'fontsize',15)
    subplot(312)
    plot(omega,amp2_d,'r','linewidth',2);hold on;
    scatter(omega,amp_d,50,'filled','MarkerFaceColor','k');
    legend('Analytical','Numerical')
    xlabel('\omega(rad/s)');ylabel('Ampitude');
    title('\omegaA_n');
    set(gca,'fontsize',15)
    subplot(313)
    plot(omega,ph2_d,'r','linewidth',2);hold on;
    scatter(omega,ph_d,50,'filled','MarkerFaceColor','k');
    legend('Analytical','Numerical')
    xlabel('\omega (rad/s)');ylabel('Phase (rad)');
    title('\Phi+\pi/2');
    ylim([-pi*1.1 pi*1.1])
    yticks([-pi 0 pi]); yticklabels({'\pi','0','\pi'})

    %% Obtain the integral of f(t)
    ftao_i = cumsum(ftao)*mean(diff(tao));
    f = zeros(size(t));
    f_i = trapz(tao,ftao_i)/T;% trapz means trapezoidal integral;
    f_i1 = trapz(tao,ftao_i)/T;
    for n=1:N
        fcos = 1.*cos(n*w0*tao);
        an_i(n)=trapz(tao,ftao_i.*fcos)*2/T;% Decomposition coefficient
        fsin = 1.*sin(n*w0*tao) ;
        bn_i(n)=trapz(tao,ftao_i.*fsin)*2/T;
        f_i = f_i + an_i(n)*cos(n*w0*t)+bn_i(n)*sin(n*w0*t); % The sum of N series.

        amp2_i(n) = sqrt(an_i(n)^2+bn_i(n)^2);
        ph2_i(n) =  atan2(-bn_i(n),an_i(n));

        amp_i(n) = amp(n)/n/w0;
        ph_i(n) = -pi/2+ph(n);

        if ph_i(n)<-pi
            ph_i(n)=ph_i(n)+2*pi;
        elseif ph_i(n)>pi
            ph_i(n)=ph_i(n)-2*pi;
        end

        if abs(bn(n))<1e-8 && abs(an(n))<1e-8
            ph2_i(n) = ph_i(n);
        end


        f_i1 = f_i1 + amp_i(n)*cos(n*w0*t + ph_i(n));
    end

    % Plot the Ampliude-Phase diagram of integral
    figure(4)
    subplot(311)
    plot([tao-T tao tao+T],[ftao_i ftao_i ftao_i],'r','linewidth',3); hold on;
    plot(t,f_i1,'k','linewidth',1.5); hold on;
    xlim([-3 3]); xlabel('t(s)'); ylabel('Ampitude'); box on;
    legend('Analytical','Numerical')
    title(['Integral of Fourier Series with N = ' num2str(N) ]);
    set(gca,'fontsize',15)

    subplot(312)
    plot(omega,amp2_i,'r','linewidth',2);hold on;
    scatter(omega,amp_i,50,'filled','MarkerFaceColor','k');
    legend('Analytical','Numerical')
    xlabel('\omega (rad/s)');ylabel('Ampitude');
    ylim([-0.01 0.4]);
    title('A_n/\omega');
    set(gca,'fontsize',15)

    subplot(313)
    plot(omega,ph2_i,'r','linewidth',2);hold on;
    scatter(omega,ph_i,50,'filled','MarkerFaceColor','k');
    legend('Analytical','Numerical')
    xlabel('\omega (rad/s)');ylabel('Phase (rad)');
    title('\Phi-\pi/2');
    ylim([-pi*1.1 pi*1.1])
    yticks([-pi 0 pi]); yticklabels({'\pi','0','\pi'})
    set(gca,'fontsize',15)

end