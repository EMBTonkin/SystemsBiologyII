function vv = properties (Ca_in,ENa,EK,q,k,T,Ca_out,Cl_out,Cl_in,F,BC)
% OC = original comment
% equations refer to Vasalou 2010, page 13, "Membrane Properties" section

Far =F(9); % Faraday Constant
Kout = F(22);% K^+ concentration in the extracellular space, Kex in paper
Naout = F(23);% Na^+ concentration in the extracellular space, Naex in paper
PCa=F(29); % Membrane permeability of Ca^2+
PCl= F(30); % Membrane permeability of Cl^-
PNa=F(32); % Membrane permeability of Na^+
R=F(34); % gas constant
vPK=F(46); % maximum PK
KPK=F(47); % PK saturation constant
nPK=F(48); % cooperativity coefficient of the PK oscillations

% OC: Permeability 
PK=(vPK*BC.^nPK./(KPK+BC.^nPK));  %OC: --> y(15,:)%PK=(.5+1*MPP.^-1./(5+MPP.^-1)); %==> Function of MP, AC: Eq (27), PK = Membrane permeability of K^+
thetaNa=exp(ENa*q/k/T/1000);
thetaK = exp(EK*q/k/T/1000);
Kin = Kout/thetaK; % K^+ concentration in the cytosol, "computed by inversion of the Nerst eq."
Nain = Naout/thetaNa; % Na^+ concentration in the cytosol, "computed by inversion of the Nerst eq."
Cl_in=Cl_in'; 

%OC: Spangler's equations
alpha = 4*PCa.*(Ca_in)*10^-3 + PK.*(Kin) + PNa.*(Nain)+PCl.*Cl_out; %Eq (24)
bitaa = PK*(Kin)-PK*Kout + PNa*(Nain) - PNa*Naout +PCl*Cl_out-PCl*(Cl_in)'; %Eq (25)
c=-(4*PCa.*(Ca_out)*10^-3 + PK*(Kout) + PNa*(Naout)+PCl*(Cl_in)'); %Eq (26)
 psi =  (-bitaa + sqrt(bitaa.^2 - 4*alpha.*c))./(2.*alpha); %Eq (23) part


 %OC: Membrane Potential
 vv= R*T/Far*2.303*log10(psi)*1000; %OC: mVolt  %AC: Eq (23), similar to Vrest