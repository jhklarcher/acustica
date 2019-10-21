%	
%	Script em MATLAB para calcular a Perda de Transmissão com um tubo de impedância
%	Quatro posições de microfone e suas funções de transferência
%	Utilizados arquivos .ita do ITA-Toolbox
%	Com base na norma ASTM E2611-09
%	
%	Autores:	Matheus Henrique Wachholz
%				José Henrique Kleinübing Larcher
%
%	Contato:	jhklarcher@gmail.com

clear all
clc

% Locais dos arquivos (ITA-Toolbox) para os mics de posição 1, 2, 3 e 4
% Terminação anecóica
func1a = '7dias - 6-3 - ap1.ita';
func2a = '7dias - 6-3 - ap2.ita';
func3a = '7dias - 6-3 - ap3.ita';
func4a = '7dias - 6-3 - ap4.ita';

% Terminação rígida / reflexiva
func1b = '7dias - 6-3 - rp1.ita';
func2b = '7dias - 6-3 - rp2.ita';
func3b = '7dias - 6-3 - rp3.ita';
func4b = '7dias - 6-3 - rp4.ita';

% Abre os arquivos como objeto do ITA-Toolbox
h1_ob_a = ita_read(func1a);
h2_ob_a = ita_read(func2a);
h3_ob_a = ita_read(func3a);
h4_ob_a = ita_read(func4a);

h1_ob_b = ita_read(func1b);
h2_ob_b = ita_read(func2b);
h3_ob_b = ita_read(func3b);
h4_ob_b = ita_read(func4b);


% Cortar o arquivo no domínio do tempo em (retirar o ruído do sinal)
h1_ob_a = ita_time_crop(h1_ob_a, [0 1.1]);
h2_ob_a = ita_time_crop(h2_ob_a, [0 1.1]);
h3_ob_a = ita_time_crop(h3_ob_a, [0 1.1]);
h4_ob_a = ita_time_crop(h4_ob_a, [0 1.1]);

h1_ob_b = ita_time_crop(h1_ob_b, [0 1.1]);
h2_ob_b = ita_time_crop(h2_ob_b, [0 1.1]);
h3_ob_b = ita_time_crop(h3_ob_b, [0 1.1]);
h4_ob_b = ita_time_crop(h4_ob_b, [0 1.1]);


% Associa o vetor no dominio da frequência para as variáveis h1, h2, h3, h4
h1_a = h1_ob_a.freq;
h2_a = h2_ob_a.freq;
h3_a = h3_ob_a.freq;
h4_a = h4_ob_a.freq;

h1_b = h1_ob_b.freq;
h2_b = h2_ob_b.freq;
h3_b = h3_ob_b.freq;
h4_b = h4_ob_b.freq;

% Dimensões
d = 2/100; % Largura da amostra [m]
l1 = 15/100; % Distância da amostra ao mic 1
s1 = 3/100; % Distância entre mic 1 e 2
l2 = (15+d)/100; % Distância da amostra ao mic 3
s2=s1; % Distância entre mic 3 e 4


%Temperatura do ar em ºC
T = 25;

%Pressão Atmosférica em Curitiba (considerado 914m de altitude)[kPa]
P = 90.8;

% Velocidade do som no ar [m/s]

c = 20.047*(273.15+T)^(1/2);
% Densidade do ar [kg/m³]
rho = 1.290*(P/101.325)*(273.15/(273.15+T));

% Faixa de frequência [Hz]
f = [20:(24000-20)/(length(h1_a)-1):24000];
f = transpose(f);

% Número de onda [rad/m]
k = (2*pi*f)/c;

% Amplitudes
A_a = 1i * ((h1_a.*exp(-1i*k*l1)-h2_a.*exp(-1i*k*(l1+s1)))./(2*sin(k*s1)));
B_a = 1i * ((h2_a.*exp(+1i*k*(l1+s1))-h1_a.*exp(+1i*k*l1))./(2*sin(k*s1)));
C_a = 1i * ((h3_a.*exp(+1i*k*(l2+s2))-h4_a.*exp(+1i*k*l2))./(2*sin(k*s2)));
D_a = 1i * ((h4_a.*exp(-1i*k*l2)-h3_a.*exp(-1i*k*(l2+s2)))./(2*sin(k*s2)));

A_b = 1i * ((h1_b.*exp(-1i*k*l1)-h2_b.*exp(-1i*k*(l1+s1)))./(2*sin(k*s1)));
B_b = 1i * ((h2_b.*exp(+1i*k*(l1+s1))-h1_b.*exp(+1i*k*l1))./(2*sin(k*s1)));
C_b = 1i * ((h3_b.*exp(+1i*k*(l2+s2))-h4_b.*exp(+1i*k*l2))./(2*sin(k*s2)));
D_b = 1i * ((h4_b.*exp(-1i*k*l2)-h3_b.*exp(-1i*k*(l2+s2)))./(2*sin(k*s2)));

% Pressão e velocidade normal incidente antes e depois da amostra
p0_a = A_a+B_a;
u0_a = (A_a-B_a)/(rho*c);
pd_a = C_a.*exp(-1i*k*d) + D_a.*exp(1i*k*d);
ud_a = (C_a.*exp(-1i*k*d) - D_a.*exp(1i*k*d))/(rho*c);

p0_b = A_b+B_b;
u0_b = (A_b-B_b)/(rho*c);
pd_b = C_b.*exp(-1i*k*d) + D_b.*exp(1i*k*d);
ud_b = (C_b.*exp(-1i*k*d) - D_b.*exp(1i*k*d))/(rho*c);

% Matriz de transferência
T_11 = (p0_a.*ud_b-p0_b.*ud_a) ./ (pd_a.*ud_b-pd_b.*ud_a);
T_12 = (p0_b.*pd_a-p0_a.*pd_b) ./ (pd_a.*ud_b-pd_b.*ud_a);
T_21 = (u0_a.*ud_b-u0_b.*ud_a) ./ (pd_a.*ud_b-pd_b.*ud_a);
T_22 = (pd_a.*u0_b-pd_b.*u0_a) ./ (pd_a.*ud_b-pd_b.*ud_a);

% Coeficiente de transmissão e Perda de Transmissibilidade
tau = 2*exp(1i*k*d)./(T_11+(T_12/(rho*c))+rho*c*T_21+T_22);
TL = -20*log10(tau);
TL = abs(TL);

figure
plot(f, TL);
title('Transmission Loss $$[dB]$$','interpreter','latex')