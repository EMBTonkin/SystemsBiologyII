function [v ]= FiringRates(Ca_in,GABA,F,CC,BC,MP,Cl_o) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
Ca_out = F(1);%15;%mM
Cl_ex=F(2); % AC: extracellular Cl^- concentration
% Cl_o=F(3); %AC: basal intracellular Cl^- concentration
Cp=F(4);
E_ex=F(5);
EK_o=F(6);
EL_o=F(7);
ENa_o=F(8);%
g_inhib= F(10); %AC: possibly g_gaba: ISPC conductance
gKo=F(11);
gNa =F(12);
k = F(13);
K_R=F(14);
KCa= F(15);
KCl1=F(16); % AC: saturation constant of PER controlled CL^- release into cytosol
KCl2=F(17); % AC: saturation constant of GABA controlled CL^- release into cytosol
Kex1=F(18);
Kex2=F(19);
Kgk=F(20);
KKCa=F(21);
nca=F(24);
nCl=F(25); % AC: cooperativity coefficient
nex1=F(26); 
nex2=F(27);
nKCa=F(28);
q = F(33);
T = F(35); % AC: Body temperature
T_room= F(36); %AC: Room temperature
V_R=F(37);
V_theta=F(38);
vCa=F(39);
vCl1=F(40); % AC: max rate of PER controlled CL^- release into cytosol
vCl2=F(41); % AC: max rate of GABA controlled CL^- release into cytosol
Vex1= F(42);
Vex2=F(43);
vgk=F(44);
vKCa=F(45);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REVERSAL POTENTIALS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ENa = ENa_o * T/(T_room);% AC: Sodium Reversal Potential
EK = (EK_o)*T/(T_room); % AC: Potassium Reversal Potential
EL =EL_o* T/(T_room); % AC: Leakage Current Reversal Potential
ECa =k*T/(2*q)*log(Ca_out./Ca_in)*1000; % AC: Calcium Reversal Potential

Cl_in=Cl_o+(MP./(KCl1+MP)*vCl1)+(GABA.^nCl./(KCl2+GABA.^nCl))*vCl2; % AC: Eq (12)
%Cl_in=Cl_o+(MP./(KCl1+MP)*vCl1)+(GABA.^nCl./(0.2+GABA.^nCl))*11.;
E_inhib = -k*T/(q)*log(Cl_ex./Cl_in)*1000'; % AC: inhibitory reversal potential (E_gaba in paper)
%   min(Cl_in)
%   max(Cl_in)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEMBRANE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vrest =properties (Ca_in,ENa,EK,q,k,T,Ca_out,Cl_ex,Cl_in,F,BC); % 
theta =Vrest + V_theta ; 
Vreset= Vrest+4;
R=V_R*(Vrest)./(K_R + Vrest);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONDUCTANCES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INa = gNa.*(Vrest-ENa);%// 36nS Jackson  %AC: Eq (5)
gK=(gKo+MP./(Kgk+MP)*vgk);
g_ex=(Vex1*abs(INa).^nex1./(Kex1+abs(INa).^nex1)+ ((Ca_in).^nex2)./(Kex2+(Ca_in).^nex2).*Vex2);
%g_ex=(Vex1*abs(INa).^nex1./(Kex1+abs(INa).^nex1)+ VEX_Ca);
gL = 1./R ;%
gCa=vCa*(MP.^nca./(KCa+MP.^nca));%
gKCa =vKCa*(CC.^nKCa./(KKCa+CC.^nKCa));%


 I_inhib = g_inhib.*(Vrest- E_inhib);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FINAL CALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AC: Vasalou Henson 2011 Supporting materials 
I_star = (-g_inhib.*E_inhib -g_ex.*E_ex+gNa.*ENa + gCa.*ECa + gK.*EK + gL.*EL+ gKCa.*EK) ;% pA (10^-9 A) %AC: Eq (7)
R_star = 1./(gNa+ gK + gL + gCa + gKCa- g_inhib - g_ex) ;% GOhm (10^9 Ohm) %AC: Eq (8)
tau_m =Cp.*R_star;% milliseconds %AC: Eq (9)

v= -(tau_m.*log( (theta-(R_star.*I_star))./(Vreset-R_star.*I_star))).^-1 ;%msec-1 %AC: Eq (6)

