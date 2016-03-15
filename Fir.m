function params1= Fir_pub
%%
Ca_out = 5;%%mM
Cl_ex=114.5; %mM
Cl_o=1; %nM
Cp=6 ;
E_ex=0;
EK_o=-97;
EL_o=-29;
ENa_o=45;
Far =96485; 
g_inhib=12.3;
gKo=9.7;
gNa =36;
k = 1.4*10^(-23) ;
K_R=34;
KCa= 22;
KCl1=4;
KCl2=1;
Kex1=574050000;
Kex2=1;
Kgk=10;
KKCa=0.16;
Kout = 1; %mM
Naout = 145; %mM
nca=2.2;
nCl=-0.2;
nex1=2.5;
nex2=-1;
nKCa=-1;
PCa=0.05;
PCl= 0.3;
PK_o=1.1;
PNa=0.036;
q = 1.6*10^(-19) ;% C
R=8.314 ;%J.K-1.mol-1
T = 37 + 273.15; % Kelvin
T_room= 22+273.15;
V_R=0.41;
V_theta=20;
vCa=12.3;
vCl1=15.5;
vCl2=19;
Vex1= 101.;
Vex2=3.5;
vgk=10;
vKCa=3;
vPK=1.9;
KPK=1;
nPK=-2;
K_BZ=16;
K_G=88000;
L0=40000;
d=0.1*10^(-1);
c=0.002;

params1=[Ca_out Cl_ex Cl_o Cp E_ex EK_o EL_o ENa_o Far g_inhib gKo gNa k K_R KCa KCl1 KCl2...
    Kex1 Kex2 Kgk KKCa Kout Naout...
    nca nCl nex1 nex2 nKCa PCa PCl PK_o PNa q R T T_room V_R V_theta...
    vCa vCl1 vCl2 Vex1 Vex2 vgk vKCa vPK KPK nPK K_BZ K_G L0 d c];
 