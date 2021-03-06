function params1= Fir

% AC: info from Vasalou 2010 paper, params summarized in table 2

Ca_out = 5;%;%mM % AC: Ca_ex in paper, Ca^2+ concentration in extracellular space
Cl_ex=114.5; % AC: extracellular Cl^- concentration, from previous studies
Cl_o=1; % AC: basal in intracellular CL^- concentration
Cp=6 ;%5.0%;--> changes from previous
E_ex=0; % AC: reversal potential of the excitatory synaptic current, assumed constant
EK_o=-97; % AC: reversal potential of potassium
EL_o=-29; % AC: Just EL in paper, resting potential of leakage current, obtained from another paper so no need to change
ENa_o=45;% % AC: Just ENa in paper, reversal potential of sodium conductance
Far =96485;% AC: faraday constant
g_inhib=12.3; % AC: g_gaba in paper, value of ISPC conductance, obtained from another paper
gKo=9.7; % AC: basal value of potassium conductance
gNa =36; % AC: values of sodium conductance
k = 1.4*10^(-23) ;
K_R=34; % AC: saturation constant of the membrane resistance oscillations
KCa= 22; % AC: saturation constant of calcium chanel dynamics
KCl1=4; % AC: saturation constant of PER controlled CL^- release into the cytosol
KCl2=1; % AC: saturation constant of GABA induced CL^- release into the cytosol
Kex1=574050000; % AC: saturation constant of AMPA- induced EPSCs
Kex2=1; % AC: saturation constant of NMDA-induced EPSCs
Kgk=10; % AC: saturation constant of potassium channel dynamics
KKCa=0.16; % AC: saturation constant of Ca^2+ activated K^+ channel dynamics
Kout = 1; % AC: Kex in paper, concentration of K^+ in extracellular space
Naout = 145; % AC: Na_ex in paper, concentration of Na^+ in extracellular space
nca=2.2; % AC: cooperativity coefficient of calcium
nCl=-0.2; % AC: cooperativity coefficient of GABA induced CL^- release into the cytosol
nex1=2.5; % AC: cooperativity coefficient of AMPA- induced EPSCs
nex2=-1; % AC: cooperativity coefficient NMDA-induced EPSCs
nKCa=-1; % AC: cooperativity coefficient of Ca^2+ activated K^+
PCa=0.05;% AC: membrane permeability of Ca^2+
PCl= 0.3;% AC: membrane permeability of Cl^-
PK_o=1.1;
PNa=0.036; % AC: membrane permeability of Na^+
q = 1.6*10^(-19) ;% C
R=8.314 ;%J.K-1.mol-1
T = 39 + 273.15; % Kelvin % AC: Body temperature
T_room= 22+273.15; % Room Temperature
V_R=0.41; % AC: maximum value of the membrane resistance oscillations
V_theta=20;
vCa=12.3; % AC: maximum rate of calcium conductance
vCl1=15.5; % AC: maximum rate of PER controlled CL^- release into the cytosol
vCl2=19; % AC: maximum rate of GABA induced CL^- release into the cytosol
Vex1= 101.;%;--> changes from previous (105 prin), AC: maximum rate of AMPA- induced EPSCs
Vex2=3.5;% AC: maximum rate of NMDA-induced EPSCs, 4.4ns in paper
vgk=10; % AC: maximum rate of potassium conductance
vKCa=3; % AC: maximum rate of Ca^2+ activated potassium conductance
vPK=1.9; % AC: maximum value of the P_K (membrane permeability of K^+) oscillations, modeled with Kuhlman et.al.
KPK=1; % AC: saturation constant of the P_K (membrane permeability of K^+) oscillations, modeled with Kuhlman et.al.
nPK=-2; % AC: cooperativity coefficient of the P_K (membrane permeability of K^+) oscillations, modeled with Kuhlman et.al.

params1=[Ca_out Cl_ex Cl_o Cp E_ex EK_o EL_o ENa_o Far g_inhib gKo gNa k K_R KCa KCl1 KCl2...
    Kex1 Kex2 Kgk KKCa Kout Naout...
    nca nCl nex1 nex2 nKCa PCa PCl PK_o PNa q R T T_room V_R V_theta...
    vCa vCl1 vCl2 Vex1 Vex2 vgk vKCa vPK KPK nPK];
 