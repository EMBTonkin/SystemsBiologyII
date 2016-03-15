function dydt= ODEs_pub(t,y,p,xa,ba,GABA_background,A,sumal,A1,sumalGABA,...
    ncell,ns,vsP0,vsB,vmB,Cl_distribution,vvip,light_hours_1, light_hours_2,...
    Coupling_start,Light_start, Shift_start, Phase_shift, End_of_simulation,l_rand,v_VIP,v_GABA,Lp_start,Lp_end)
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
k1=p(1);
k2 = p(2);
k3 = p(3);
k4 = p(4);
k5 = p(5);
k6 = p(6);
k7 = p(7);
k8 = p(8);
KAP = p(9); %%%made variable with IN to simulate competitive inhibition
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
%vmB= p(48);
vmC = p(49);
vmP = p(50);
%vsB = p(51);
vsC = p(52);
KD = p(55);
vP = p(59);
VMK=p(60);
K_1 = p(62);
K_2 = p(63);
WT = p(64);
CT = p(65);
KC = p(66);
%vsP0 = p(67);
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
K_b=p(81);
%v_GABA=p(82)*GABA_switch;
KGABA=p(83);
nGABA=p(84);
kdGABA=p(85);
ndGABA=p(86);
%v_VIP=p(87);
KVIP=p(88);
nVIP=p(89);
kdVIP=p(90);
ndVIP=p(91);
vGlu=p(92);
KGlu=p(93);
KGluR=p(94);
vGluR=p(95);
vv1=p(96);
vkk=p(97);
nkk=p(98);
Kkk=p(99);
v_Ca=p(100);
vvo=p(101);
nvo=p(102);
Kvo=p(103);


%%%%%% Light switch  %%%%%%%%%%%%%%%%%%%%

%%%dark if 0, light if 1
%%%Light switch starts with the full light period
%%%If using a value for Shift_start evenly divisible by the light period (24h),
%%%Negative phase shifts interrupt night period, and positive phase shifts
%%%lengthen light periods.  A Shift_start offset by 12h would inversely
%%%have negative phase shifts interrupt the first light period and
%%% positive phase shifts would lengthen the night period.  In general,
%%% positive phase shifts generate phase delays and negative phase shifts
%%% create a phase advance.
Light_period=24;
if  (t>=Light_start) && (t<=Shift_start+Phase_shift)
    if (t-Light_start)/Light_period-fix((t-Light_start)/Light_period)-light_hours_1/Light_period>=0
        daylight=0;
    else
        daylight=1;
    end
elseif (t>Shift_start+Phase_shift)
     if (t-Light_start-Phase_shift)/Light_period-fix((t-Light_start-Phase_shift)/Light_period)-light_hours_2/Light_period>=0
         daylight=0;
     else
         daylight=1;
     end
else
daylight=0;
end

if (t>=Lp_start) && (t<Lp_end)  %%%Dark pulse
    daylight=0; %%%If using light pulse instead of dark pulse, set to 1
end

if t<Coupling_start
    switz=0;
    beta1=zeros(ncell,1);
    S_GABA=GABA_background*ones(ncell,1);
else
    switz=1;
    %%%%% VIP calculations %%%%%%%%%%%%%
    VIP=vVIP';
    VIP=VIP(ones(1,ncell),:);
    S_VIP=sum(VIP.*A,2)'.*sumal+vvip;
    if daylight==0
        beta1 = (0+S_VIP)'./(KD+(0+S_VIP))';
    elseif daylight==1
        beta1 = 1*l_rand+(0+S_VIP)'./(KD+(0+S_VIP))'.*(1-l_rand);
    end
    %%% GABA calculations%%%%%%% 
     gGABA(gGABA<0)=0;
     GABA=GABA_background+gGABA';
     GABA=GABA(ones(1,ncell),:);
     S_GABA=ba+sum(GABA.*A1,2).*sumalGABA'; % GABA binding to cell surface
end
S_GABA(S_GABA==0)=GABA_background;  %%%Keeps calculations finite when no GABA connections are present

wow=0;

if Light_start>End_of_simulation
    vo=vvo*BC.^nvo./(Kvo+BC.^nvo); %%%calcium influx term
else
    %%%Glutamate modeling
    Glu=vGlu*MP./(KGlu+MP);
    
    if daylight==0
        bGluR=Glu./(KGluR+Glu);
    elseif daylight==1  
        bGluR=1*l_rand+Glu./(KGluR+Glu).*(1-l_rand);
    end
    vo=vGluR*bGluR;  %%%calcium influx term
end

vK = switz*(VMK.*(xa*Ca+wow)./(kMK+(xa*Ca+wow))+V_b*beta1./(K_b+beta1));

vsPc=vsP0+CT.*CB./(KC+CB);

vv =FiringRates_pub((1*Ca+wow),S_GABA, Fir_pub,CC,BC,MP,Cl_distribution);

vv2 =VM2 *(Ca.^M2)./(K2^M2+Ca.^M2); %calcium from cytosol to ryanodine stores
vv3 =1*(VM3.*(Ca_store.^M3)./(KR^M3+Ca_store.^M3)).*(Ca.^pA)./(KA^pA+Ca.^pA);%NJK: calcium from ryanodine stores to cytosol (source: CV's Plos Eq (17))

kk1=vkk*CC.^nkk./(Kkk+CC.^nkk); %calcium efflux

i = 1:ncell;

dydt((i-1)*ns+1,1) = vo + vv1*IP3-vv2+vv3 + kf*Ca_store ...
    - kk1.*Ca.^(v_Ca); 

dydt((i-1)*ns+2,1) = vv2-vv3-kf*Ca_store; % Ca in non_IP3 store-0.01*Ca_store

dydt((i-1)*ns+3,1) = (vsPc.*BN.^n)./(KAP.^n+BN.^n)-vmP.*MP./(KmP+MP)-kdmp*MP;

dydt((i-1)*ns+4,1) = (vsC.*BN.^n)./(KAC^n+BN.^n)-vmC.*MC./(KmC+MC)-kdmc*MC;

dydt((i-1)*ns+5,1) = KIB^m*vsB./(KIB^m+BN.^m)-(vmB.*MB)./(KmB+MB)-kdmb*MB;

dydt((i-1)*ns+6,1) = ksP*MP - V1P*PC./(Kp+PC)+ ...
    V2P*PCP./(Kdp+PCP) + k4*PCC -k3*PC.*CC - kdn*PC;

dydt((i-1)*ns+7,1) = ksC*MC - V1C*CC./(Kp+CC)+ ...
    V2C*CCP./(Kdp+CCP) + k4*PCC - k3*PC.*CC - kdnc*CC;

dydt((i-1)*ns+8,1) = V1P*PC./(Kp+PC) - V2P*PCP./(Kdp+PCP)- ...
    vdPC*PCP./(Kd+PCP) - kdn*PCP;

dydt((i-1)*ns+9,1) = V1C*CC./(Kp+CC) - V2C*CCP./(Kdp+CCP) - ...
    vdCC*CCP./(Kd+CCP) - kdn*CCP;

dydt((i-1)*ns+10,1) = -V1PC*PCC./(Kp+PCC)+ ...
    V2PC*PCCP./(Kdp+PCCP) - k4*PCC + k3*PC.*CC + ...
    k2*PCN - k1*PCC - kdn*PCC;

dydt((i-1)*ns+11,1) = -V3PC*PCN./(Kp+PCN)+ ...
    V4PC*PCNP./(Kdp+PCNP) - k2*PCN + k1*PCC- ...
    k7*BN.*PCN + k8*IN - kdn*PCN;

dydt((i-1)*ns+12,1) = V1PC*PCC./(Kp+PCC)- ...
    V2PC*PCCP./(Kdp+PCCP)-vdPCC*PCCP./(Kd+PCCP)-kdn*PCCP;

dydt((i-1)*ns+13,1) = V3PC*PCN./(Kp+PCN)- ...
    V4PC*PCNP./(Kdp+PCNP)-vdPCN*PCNP./(Kd+PCNP)-kdn*PCNP;

dydt((i-1)*ns+14,1) = ksB*MB - V1B*BC./(Kp+BC)+ ...
    V2B*BCP./(Kdp+BCP) - k5*BC + k6*BN - kdn*BC;

dydt((i-1)*ns+15,1) = V1B*BC./(Kp+BC)- ...
    V2B*BCP./(Kdp+BCP) - vdBC*BCP./(Kd+BCP) - kdn*BCP;

dydt((i-1)*ns+16,1) = -V3B*BN./(Kp+BN)+ ...
    V4B*BNP./(Kdp+BNP) + k5*BC - k6*BN - k7*BN.*PCN + k8*IN - kdn*BN;

dydt((i-1)*ns+17,1) = V3B*BN./(Kp+BN)- ...
    V4B*BNP./(Kdp+BNP) - vdBN*BNP./(Kd+BNP) - kdn*BNP;

dydt((i-1)*ns+18,1) = -k8*IN + k7*BN.*PCN - vdIN*IN./(Kd+IN) - kdn*IN;

dydt((i-1)*ns+19,1) = vP/WT*(vK/vP.*(1-CB)./(K_1+1-CB) - CB./(K_2+CB));

dydt((i-1)*ns+20,1) = (v_VIP*vv.^nVIP)./(KVIP+vv.^nVIP) - kdVIP*vVIP.^ndVIP;

dydt((i-1)*ns+21,1) = (v_GABA.*vv.^nGABA)./(KGABA+vv.^nGABA) - kdGABA*gGABA.^ndGABA;