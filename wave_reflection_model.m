clear

%% Parametros iniciais
hr = 25; % Altura da antena receptora em metros
e = 15; % Permissividade eletrica relativa do solo comum
e0 = 8.854*10^(-12); % Permissividade eletrica do vacuo em F/m
sigma = 5*10^(-3); % Condutividade da superficie refletora (solo comum) em Siemens/m
c = 299792458; % Velocidade da luz no vacuo em m/s
Freq = 2250000000; % Frequencia da onda em Hz
lambda = c/Freq; % Comprimento de onda em metros
omega = 2*pi*Freq; % Frequência angular
et = e-1j*sigma/(omega*e0); % Permissividade eletrica efetiva da superficie refletora

Pt = 20; % Potência do transmissor
VSWR = 1.5; % Voltage Standing Wave Radio
gama = 0.007; % Atenuacao especifica da atmosfera para freq = 2.25 GHz em dB/km
k_boltz = 1.38*10^(-23); % Constante de Boltzmann
BW = 2.32*10^6; % Largura de banda
Te = 211.4; % Temperatura de ruído na entrada do LNA

close all
clear alpha beta cos_theta_t delta_phi f f_com_abs f_ideal Lespaco Pr_dBm R1 R2 r_p SNR theta
% d : Distancia horizontal entre aeronave e antena receptora em metros
% ht: Altura da aeronave em metros
d = 100:100:750000;
ht = 40:20:8000;

%% Para plotagem
for i=1:length(d)
    for k=1:length(ht)
        R1(i,k) = sqrt(d(i)^2+(ht(k)-hr)^2); % Distancia da visada direta entre aeronave e antena receptora em metros
        R2(i,k) = sqrt(d(i)^2+(ht(k)+hr)^2); % Distancia percorrida pelo raio refletido em solo entre aeronave e antena receptora
        theta(i,k) = atan(d(i)/(hr+ht(k))); % Angulo de incidencia da reflexao com a vertical em radianos
        alpha(i,k) = atan((ht(k)-hr)/d(i)); % Angulo entre o raio da visada direta e a horizontal em radianos
        beta(i,k) = atan((ht(k)+hr)/d(i)); % Angulo entre o raio refletido e a horizontal em radianos
    end
end

delta_phi = 2*pi*(R2-R1)/lambda; % Diferença de fase
cos_theta_t = et-sin(theta).^2;
r_p = (et*cos(theta)-sqrt(cos_theta_t))./(et*cos(theta)+sqrt(cos_theta_t)); % coeficiente de reflexão relativo efetivo
f_com_abs_com_ref = 10*log10((lambda^2./(4*pi*R1).^2).*abs(1+r_p.^2.*R1.^2./R2.^2.*exp(1j*2*delta_phi).*10.^(-gama*(R2-R1)/10000)+2*R1./R2.*r_p.*exp(1j*delta_phi).*cos(alpha+beta).*10.^(-gama*(R2-R1)/20000)))-gama*R1/1000; % perda com absorção e reflexão
f_ideal = 10*log10(lambda^2./((4*pi)^2*R1.^2)); % Perda no espaço sem reflexão e sem absorção
f_sem_ref = 10*log10(lambda^2./((4*pi)^2*R1.^2)) - gama*R1/1000; % Perda no espaço sem reflexão

figure;
CO1(:,:,1) = zeros(length(d),length(ht)); % red
CO1(:,:,2) = zeros(length(d),length(ht)); % green
CO1(:,:,3) = ones(length(d),length(ht)); % blue
mesh(ht,d,f_com_abs_com_ref,CO1,'EdgeAlpha',0.15,'FaceAlpha',0.5) % Plota modelo mais completo
hold on
CO2(:,:,1) = 0.7*ones(length(d),length(ht)); % red
CO2(:,:,2) = zeros(length(d),length(ht)); % green
CO2(:,:,3) = zeros(length(d),length(ht)); % blue
mesh(ht,d,f_sem_ref,CO2) % Plota modelo do espaco livre com absorção
title('Modelos de perda devido ao caminho')
xlabel('h_{t} (m)')
ylabel('d (m)')
zlabel('-Perda (dB)')

%% Potencia recebida final
Pt_dBm = 10*log10(Pt)+30; % Potencia do transmissor
Gt = -13; % Ganho da antena transmissora
Lcon1TX = 0.1; % Perda de conector na transmissão
LcaboTX = 0.42235406; % Perda de cabo na transmissão
Lcon2TX = 0.1; % Perda de conector na transmissão
Lmismatch = -10*log10(1-((VSWR-1)/(VSWR+1))^2);
Lespaco = -f_com_abs_com_ref;
Lpol = 4; % Perda de diferença de polarição entre antenas
Ldl = 1; % Perda de desalinhamento entre antenas
% Labs esta contabilizado em f_com_abs_com_ref
Gr = 31.7; % Ganho da antena receptora
Gant = 22.96; % Ganho dos componentes da antena
L1 = 0.1; % Perda de conector na recepção
L3 = 6.39296045; % Perda de cabo na recepção
L2 = 0.4; % Perda de conectores na recepção

Pr_dBm = Pt_dBm+Gt-Lcon1TX-LcaboTX-Lcon2TX-Lmismatch-Lespaco-Lpol-Ldl+Gr+Gant-L1-L3-L2; % Potencia no receptor em dBm

%% Potencia do ruido na saida do receptor em dBm
Pn = 10*log10(k_boltz*BW*Te)+Gant-L1-L3-L2+30;
SNR = Pr_dBm - Pn;

% Plota contorno
p1=figure;
B = SNR >= 12;
contourf(d/1000,ht,B.');
hold on
Rterra = 6400000;
dterra = Rterra * (asin(sqrt(1-(Rterra./(Rterra + ht)).^2)) + asin(sqrt(1-(Rterra./(Rterra + hr)).^2)));
plot(dterra/1000,ht,'--r','LineWidth',3)
ylabel('h_{t} (m)')
xlabel('d (km)')
title('Envelope de operação para SNR \geq 12 dB com reflexão')
grid on
print('Envelope de Operação reflexao','-dpng')

% Plota contorno
p=figure;
B_sem_ref = SNR + Lespaco + f_sem_ref >= 12;
contourf(d/1000,ht,B_sem_ref.');
hold on
plot(dterra/1000,ht,'--r','LineWidth',3)
ylabel('h_{t} (m)')
xlabel('d (km)')
title('Envelope de operação para SNR \geq 12 dB sem reflexão')
grid on
print('Envelope de Operação ideal','-dpng')

