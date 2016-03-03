function dydt= ODEs(t,y,p,xa,ba,A,sumal,A1,sumalGABA,ncell,ns,vsP0,vsB,vmB,Cl_o,vvip,t1)
dydt=zeros(ns*ncell,1);

i=1:ncell;
Ca=y((i-1)*ns+1);
Ca_store =  y((i-1)*ns+2);
MP = y((i-1)*ns+3);
MC = y((i-1)*ns+4);
MB = y((i-1)*ns+5);
PC = y((i-1)*ns+6);
CC = y((i-1)*ns+7);
PCP = y((i-1)*ns+8);
CCP = y((i-1)*ns+9);
PCC = y((i-1)*ns+10);
PCN = y((i-1)*ns+11);
PCCP = y((i-1)*ns+12);
PCNP = y((i-1)*ns+13);
BC = y((i-1)*ns+14);
BCP = y((i-1)*ns+15);
BN = y((i-1)*ns+16);
BNP = y((i-1)*ns+17);
IN = y((i-1)*ns+18);
CB = y((i-1)*ns+19);
vVIP = y((i-1)*ns+20);
gGABA=y((i-1)*ns+21);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%
k1 = p(1);
k2 = p(2);
k3 = p(3);
k4 = p(4);
k5 = p(5);
k6 = p(6);
k7 = p(7);
k8 = p(8);
KAP = p(9);
KAC = p(10);
KIB = p(11);
kdmb = p(12);
kdmc = p(13);
kdmp = p(14);
kdnc = p(15);
kdn = p(16);
Kd = p(17);
Kdp = p(18);
Kp = p(19);
KmB = p(20);
KmC = p(21);
KmP = p(22);
ksB = p(23);
ksC = p(24);
ksP = p(25);
n = p(26);
m = p(27);
Vphos = p(28);
V1P = Vphos;
V1PC = Vphos;
V3PC = Vphos;
V1B = p(29);
V1C = p(30);
V2B = p(33);
V2C = p(34);
V2P = p(35);
V2PC = p(36);
V3B = p(37);
V4B = p(39);
V4PC = p(40);
vdBC = p(41);
vdBN = p(42);
vdCC = p(43);
vdIN = p(44);
vdPC = p(45);
vdPCC = p(46);
vdPCN = p(47);
%vmB= p(48);%
vmC = p(49);
vmP = p(50);
%vsB = p(51);%1;
vsC = p(52);
KD = p(55);%2;%*scale;
vP = p(59);
VMK=p(60);
K_1 = p(62);
K_2 = p(63);
WT = p(64);
CT = p(65);
KC = p(66);
%vsP0 = p(67);%1.1;
kf=p(69);
IP3=p(70);
VM3= p(71);
M3= p(72);
KR= p(73);
KA = p(74);
pA = p(75);
VM2=p(76);
K2=p(77);
M2= p(78);
kMK=p(79);
V_b=p(80);
k_b=p(81);

%%%%%%%%%%%% Neurotransmitters %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k = 1.4*10^(-23) ;
% T = 37 + 273.15; % Kelvin
% q = 1.6*10^(-19) ;% C
% Ca_out = 5;%;%mM
% Kex2=1;
% nex2=-1;
% Vex2=3.5;%
% nex1=2.5;
% Vex1= 101.;%
% Kex1=574050000;



if t<0
  switz=0;
  beta1=zeros(ncell,1);
  S_GABA=0.1*ones(ncell,1);
else 
    switz=5;
    %%%%% VIP calcs %%%%%%%%%%%%%
    VIP=vVIP';
    VIP=VIP(ones(1,ncell),:);
    S_VIP=sum(VIP.*A,2)'.*sumal+vvip;
    clear VIP
    beta1 = (0+S_VIP)'./(KD+(0+S_VIP))';

%%% GABA calcs%%%%%%%
    GABA=0.1+gGABA';
    GABA=GABA(ones(1,ncell),:);
    S_GABA=ba+sum(GABA.*A1,2).*sumalGABA'; % GABA binding to cell surface

end

%%%% Phase Response 
% if t>t1 && t<t1+2
% S_GABA=ba+sum(GABA.*A1,2).*sumalGABA'+10^5;
%  %   beta1 = (100+S_VIP)'./(KD+(100+S_VIP))';
% end
% %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if t>t1 && t<t1+2
%     wow =0.1;
% else
     wow=0;
% end
[vv]= FiringRates((xa*Ca+wow),S_GABA, Fir,CC,BC,MP,Cl_o);
%mean(vv)


vv3 =1*(VM3.*(Ca_store.^M3)./(KR^M3+Ca_store.^M3)).*(Ca.^pA)./(KA^pA+Ca.^pA);%
vv2 =VM2 *(Ca.^M2)./(K2^M2+Ca.^M2);
kk1=3.3*CC.^0.1./(.02+CC.^0.1);
vo=(0.09*BC.^4.5./(4.5+BC.^4.5));

vK = switz*(VMK.*(xa*Ca+wow)./(kMK+(xa*Ca+wow))+V_b*beta1./(k_b+beta1));%
vsPc =vsP0+ CT.*CB./(KC+CB);


for i =1:ncell
    dydt((i-1)*ns+1,1) = vo(i)+ 0.0003*IP3-vv2(i)+vv3(i)+kf*Ca_store(i) - kk1(i)*Ca(i).^(2); %19*Ca.^(4.6); Ca_intracellular+0.01*Ca_store - 30*Ca^(4.6)); %
    dydt((i-1)*ns+2,1) = (vv2(i)-vv3(i)-kf*Ca_store(i)); % Ca in non_IP3 store-0.01*Ca_store
    dydt((i-1)*ns+3,1) = vsPc(i)*BN(i)^n/(KAP^n+BN(i)^n)-vmP*MP(i)/(KmP+MP(i))-kdmp*MP(i);
    dydt((i-1)*ns+4,1) = vsC*BN(i)^n/(KAC^n+BN(i)^n)-vmC*MC(i)/(KmC+MC(i))-kdmc*MC(i);
    dydt((i-1)*ns+5,1) = vsB(i)*KIB^m/(KIB^m+BN(i)^m)-vmB(i)*MB(i)/(KmB+MB(i))-kdmb*MB(i);
    dydt((i-1)*ns+6,1) = ksP*MP(i)-V1P*PC(i)/(Kp+PC(i))+V2P*PCP(i)/(Kdp+PCP(i))+k4*PCC(i)-k3*PC(i)*CC(i)-kdn*PC(i);
    dydt((i-1)*ns+7,1) = ksC*MC(i)-V1C*CC(i)/(Kp+CC(i))+V2C*CCP(i)/(Kdp+CCP(i))+k4*PCC(i)-k3*PC(i)*CC(i)-kdnc*CC(i);
    dydt((i-1)*ns+8,1) = V1P*PC(i)/(Kp+PC(i))-V2P*PCP(i)/(Kdp+PCP(i))-vdPC*PCP(i)/(Kd+PCP(i))-kdn*PCP(i);
    dydt((i-1)*ns+9,1) = V1C*CC(i)/(Kp+CC(i))-V2C*CCP(i)/(Kdp+CCP(i))-vdCC*CCP(i)/(Kd+CCP(i))-kdn*CCP(i);

    dydt((i-1)*ns+10,1) = -V1PC*PCC(i)/(Kp+PCC(i))+V2PC*PCCP(i)/(Kdp+PCCP(i))-k4*PCC(i)+k3*PC(i)*CC(i)+k2*PCN(i)-k1*PCC(i)-kdn*PCC(i);
    dydt((i-1)*ns+11,1) = -V3PC*PCN(i)/(Kp+PCN(i))+V4PC*PCNP(i)/(Kdp+PCNP(i))-k2*PCN(i)+k1*PCC(i)-k7*BN(i)*PCN(i)+k8*IN(i)-kdn*PCN(i);
    dydt((i-1)*ns+12,1) = V1PC*PCC(i)/(Kp+PCC(i))-V2PC*PCCP(i)/(Kdp+PCCP(i))-vdPCC*PCCP(i)/(Kd+PCCP(i))-kdn*PCCP(i);
    dydt((i-1)*ns+13,1) = V3PC*PCN(i)/(Kp+PCN(i))-V4PC*PCNP(i)/(Kdp+PCNP(i))-vdPCN*PCNP(i)/(Kd+PCNP(i))-kdn*PCNP(i);

    dydt((i-1)*ns+14,1) = ksB*MB(i)-V1B*BC(i)/(Kp+BC(i))+V2B*BCP(i)/(Kdp+BCP(i))-k5*BC(i)+k6*BN(i)-kdn*BC(i);
    dydt((i-1)*ns+15,1) = V1B*BC(i)/(Kp+BC(i))-V2B*BCP(i)/(Kdp+BCP(i))-vdBC*BCP(i)/(Kd+BCP(i))-kdn*BCP(i);
    dydt((i-1)*ns+16,1) = -V3B*BN(i)/(Kp+BN(i))+V4B*BNP(i)/(Kdp+BNP(i))+k5*BC(i)-k6*BN(i)-k7*BN(i)*PCN(i)+k8*IN(i)-kdn*BN(i);
    dydt((i-1)*ns+17,1) = V3B*BN(i)/(Kp+BN(i))-V4B*BNP(i)/(Kdp+BNP(i))-vdBN*BNP(i)/(Kd+BNP(i))-kdn*BNP(i);
    dydt((i-1)*ns+18,1) = -k8*IN(i)+k7*BN(i)*PCN(i)-vdIN*IN(i)/(Kd+IN(i))-kdn*IN(i);
    dydt((i-1)*ns+19,1) = (vP/WT)*(vK(i)/vP*(1-CB(i))/(K_1+1-CB(i))-CB(i)/(K_2+CB(i)));
    dydt((i-1)*ns+20,1) = 0.5*vv(i).^1.9/(20+vv(i).^1.9)-0.5*vVIP(i).^0.2;%
    dydt((i-1)*ns+21,1) = 0.5*vv(i).^1.9/(20+vv(i).^1.9)-0.5*gGABA(i).^0.2;%
%    dydt((i-1)*ns+21,1) = .75*vv(i).^1.9/(15+vv(i).^1.9)-0.5*gGABA(i).^0.2;%

end




