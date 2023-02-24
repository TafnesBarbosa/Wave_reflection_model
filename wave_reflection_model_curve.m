clear

%% Parametros iniciais
hrl = 25; % Altura da antena receptora em metros
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
% dl : Distancia horizontal entre aeronave e antena receptora em metros
% (paralela à superfície da Terra)
% htl: Altura da aeronave em metros
dl = 100:50:350000;
htl = 40:20:8000;
R = 6400000; % Raio da Terra em metros
at = R./(htl + R);
ar = R./(hrl + R);
b = dl / R;

for m=1:length(dl)
    for n=1:length(htl)
        if at(n) < cos(b(m)-acos(ar)) % Tem visada direta e reflete
            poly = [-4/3 + at(n)/6 + ar/6, 2*b(m) - b(m)*ar/2, 2 - b(m)^2 - at(n) - ar + ar*b(m)^2/2, -b(m) + ar*b(m) - ar*b(m)^3/6];
            polyroots = roots(poly);
            if polyroots(1) > 0 && b(m) - polyroots(1) > 0
                ang = polyroots(1);
            elseif polyroots(2) > 0 && b(m) - polyroots(2) > 0
                ang = polyroots(2);
            else
                ang = polyroots(3);
            end
            ht = (htl(n) + R)*cos(b(m) - ang) - R;
            hr = (hrl + R)*cos(ang) - R;
            d = (R + ht)*tan(b(m) - ang) + (R + hr)*tan(ang);
            
            R1(m,n) = sqrt(d^2+(ht-hr)^2); % Distancia da visada direta entre aeronave e antena receptora em metros
            R2(m,n) = sqrt(d^2+(ht+hr)^2); % Distancia percorrida pelo raio refletido em solo entre aeronave e antena receptora
            theta(m,n) = atan(d/(hr+ht)); % Angulo de incidencia da reflexao com a vertical em radianos
            alpha(m,n) = atan((ht-hr)/d); % Angulo entre o raio da visada direta e a horizontal em radianos
            beta(m,n) = atan((ht+hr)/d); % Angulo entre o raio refletido e a horizontal em radianos
        elseif at(n) == cos(b(m)-acos(ar)) % Tem visada direta e não reflete
            ht = 0;
            hr = 0;
            d = sqrt((R+htl(n))^2-R^2) + sqrt((R+hr)^2-R^2);
            
            R1(m,n) = sqrt(d^2+(ht-hr)^2); % Distancia da visada direta entre aeronave e antena receptora em metros
            R2(m,n) = Inf; % Distancia percorrida pelo raio refletido em solo entre aeronave e antena receptora
            theta(m,n) = pi/2; % Angulo de incidencia da reflexao com a vertical em radianos
            alpha(m,n) = 0; % Angulo entre o raio da visada direta e a horizontal em radianos
            beta(m,n) = 0; % Angulo entre o raio refletido e a horizontal em radianos
        else % Não tem visada direta
            R1(m,n) = Inf; % Distancia da visada direta entre aeronave e antena receptora em metros
            R2(m,n) = Inf; % Distancia percorrida pelo raio refletido em solo entre aeronave e antena receptora
            theta(m,n) = pi/2; % Angulo de incidencia da reflexao com a vertical em radianos
            alpha(m,n) = 0; % Angulo entre o raio da visada direta e a horizontal em radianos
            beta(m,n) = 0; % Angulo entre o raio refletido e a horizontal em radianos
        end
    end
    if mod(m,100) == 0
        m
    end
end

%% Para plotagem
delta_phi = 2*pi*(R2-R1)/lambda; % Diferença de fase
delta_phi(isnan(delta_phi)) = 0; % Remove NANs
cos_theta_t = et-sin(theta).^2;
r_p = (et*cos(theta)-sqrt(cos_theta_t))./(et*cos(theta)+sqrt(cos_theta_t)); % coeficiente de reflexão relativo efetivo
f_com_abs_com_ref = 10*log10((lambda^2./(4*pi*R1).^2).*abs(1+r_p.^2.*R1.^2./R2.^2.*exp(1j*2*delta_phi).*10.^(-gama*(R2-R1)/10000)+2*R1./R2.*r_p.*exp(1j*delta_phi).*cos(alpha+beta).*10.^(-gama*(R2-R1)/20000)))-gama*R1/1000; % perda com absorção e reflexão
f_com_abs_com_ref(isnan(f_com_abs_com_ref)) = -Inf; % Remove NANs
f_ideal = 10*log10(lambda^2./((4*pi)^2*R1.^2)); % Perda no espaço sem reflexão e sem absorção
f_sem_ref = 10*log10(lambda^2./((4*pi)^2*R1.^2)) - gama*R1/1000; % Perda no espaço sem reflexão

clear R1 R2 theta alpha beta

figure;
CO1(:,:,1) = zeros(length(dl),length(htl)); % red
CO1(:,:,2) = zeros(length(dl),length(htl)); % green
CO1(:,:,3) = ones(length(dl),length(htl)); % blue
mesh(htl,dl,f_com_abs_com_ref,CO1,'EdgeAlpha',0.0025,'FaceAlpha',0.05) % Plota modelo mais completo
hold on
CO2(:,:,1) = 0.7*ones(length(dl),length(htl)); % red
CO2(:,:,2) = zeros(length(dl),length(htl)); % green
CO2(:,:,3) = zeros(length(dl),length(htl)); % blue
mesh(htl,dl,f_sem_ref,CO2) % Plota modelo do espaco livre com absorção
title('Modelos de perda devido ao caminho')
xlabel('h_{t} (m)')
ylabel('d (m)')
zlabel('-Perda (dB)')
print('perda_reflexao','-dpng')
close

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
contourf(dl/1000,htl,B.');
hold on
Rterra = 6400000;
dterra = Rterra * (asin(sqrt(1-(Rterra./(Rterra + htl)).^2)) + asin(sqrt(1-(Rterra./(Rterra + hrl)).^2)));
plot(dterra/1000,htl,'--r','LineWidth',3)
ylabel('h_{t} (m)')
xlabel('d (km)')
title('Envelope de operação para SNR \geq 12 dB com reflexão')
grid on
print('snr_com_ref','-dpng')
close

% Plota contorno
p=figure;
B_sem_ref = SNR + Lespaco + f_sem_ref >= 12;
contourf(dl/1000,htl,B_sem_ref.');
hold on
plot(dterra/1000,htl,'--r','LineWidth',3)
ylabel('h_{t} (m)')
xlabel('d (km)')
title('Envelope de operação para SNR \geq 12 dB sem reflexão')
grid on
print('snr_sem_ref','-dpng')
close
