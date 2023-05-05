%datos del problema
f=60e6;
B=1e6;
Vceq=10;
Icq=6e-3;

b11=6j*1e-3;
b12= - 2*pi*60e6*0.25e-12j;
b21=-35j*1e-3;
b22=450j*1e-6;
g11=5e-3;
g12=0;              %falta colocar valor
g22=70e-6;
g21=175e-3;

y11=b11+g11;
y12=b12+g12;
y21=b21+g21;
y22=b22+g22;

y11
y22
y12
y21

%CALCULO DE ESTABILIDAD 
%Y21 -> Yf .... Y12 -> Yr

C = (abs(y21*y12))/((2*g11*g22)-(real(y21*y12)))

%ELECCION DE K
K = 3

%CALCULO DE YL YS

%YS YL
gs = (sqrt((K*(abs(y12*y21)+real(y12*y21)))/(2))*sqrt(g11/g22))-g11
gl = (sqrt((K*(abs(y12*y21)+real(y12*y21)))/(2))*sqrt(g22/g11))-g22






bl=-i*imag(y22);

for a=5:-1:0
yl=gl+bl;
y1=y11-((y12*y21)/(y22+yl));
bs=-i*imag(y1);
ys=gs+bs;
y2=y22-((y12*y21)/(y11+ys));
bl=-i*imag(y2);
end


ys
yl
Zs=1/ys
Zl=1/yl

%MAG
MAG = ((abs(y21)^2))/(4*real(y11)*real(y22))
MAGdb = 10*log10(MAG)

%GT
GT = (4*real(ys)*real(yl)*(abs(y21))^2)/((abs((y11+ys)*(y22+yl)-(y12*y21)))^2)
GTdb = 10*log10(GT)


% Calculo de red de polarizacion
% Considero que Vb=2V
Re=(Vb-0.7)/Icq
Rb=(hfe*Re)/10
R1=(Rb*Vcc)/Vb
R2=(Rb*Vcc)/(Vcc-Vb)
Rc=(Vcc-Vceq-Icq*Re)/Icq