% parametros S 

%S11
mod=0.37;
ang=270;
s11=mod*cosd(ang)+i*sind(ang)*mod;

%S22
mod=0.76;
ang=343;
s22=mod*cosd(ang)+i*sind(ang)*mod;

%S12
mod=0.046;
ang=66;
s12=mod*cosd(ang)+i*sind(ang)*mod;

%S21
mod=4.5;
ang=96;
s21=mod*cosd(ang)+i*sind(ang)*mod;

%---------------------------------%

%determinante
det = s11*s22-s12*s21;

%Valor de K rollet
%   K < 1 Pot inestable 
K = (1-abs(s11)^2-abs(s22)^2+abs(det)^2)/(2*abs(s12*s21))

%----------------------------%

%Circulo estabilidad Entrada
c1=s11-det*conj(s22)
D=(((abs(s11))^2)-((abs(det))^2)) 

radioentrada = abs((s12*s21)/((abs(s11)^2)-(abs(det)^2)))
centroentrada = (conj(c1))/D

%Circulo estabilidad Salida
radiosalida = abs(s12*s21/(((abs(s22)^2))-(abs(det)^2)))
centrosalida = conj((s22-det*conj(s11)))/((abs(s22)^2)-(abs(det)^2))


figure();
smithchart();
hold on;

%entrada --> |s11| < 1 -> por fuera estable
viscircles([real(real(centroentrada)) imag(imag(centroentrada)*i)], [radioentrada], 'color', 'b')

%salida --> |s22| < 1 -> por fuera estable
viscircles([real(real(centrosalida)) imag(imag(centrosalida)*i)], [radiosalida], 'color', 'g')


%CIRCULOS DE GANANCIA CONSTANTE

%SI ROLLET ES < 1 NO SE PUEDE CALCULAR LA MAG
%SE USA MSG = |S21|/|S12|
MSG = abs(s21)/abs(s12);
MSGdb = 10*log10(MSG)
%para darle margen 
MSGdb = MSGdb - 2
MSG = 10^(MSGdb/10)
%ganancia deseada en veces
Gain = MSG; 
D2 = (abs(s22))^2-(abs(det))^2;
C2 = s22-(det*conj(s11));
G = Gain / (abs(s21))^2;
%centro
c = (G*conj(C2))/(1+D2*G);
%radio
r = (sqrt(1-(2*K*abs(s12*s21)*G)+((abs(s12*s21))^2)*G^2))/(1+D2*G);
viscircles([real(real(c)) imag(imag(c)*i)], [r], 'color', 'r')

%ELIJO UN PUNTO SOBRE EL CIRCULO DE G-CTE

%despues de buscarlo en el circulo de ganancia

rl = 0.6096 + 0.3385i;
%En polar  --> rl = 0.6972|29.042° 
rs = conj(s11+((s12*s21*rl)/(1-(rl*s22))))
%En polar  --> rs = 0.5567|118.63°

Gt = ((abs(s21)^2)*(1-abs(rs)^2)*(1-abs(rl)^2))/(abs((1-s11*rs)*(1-s22*rl)-(s12*s21*rl*rs))^2)
Gtdb = 10*log10(Gt)