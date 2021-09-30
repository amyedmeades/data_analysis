%% 0%PEGDA/100%EMIMTF2N from 05-2019 me
%all of the kinetics stuff is technically *not* necessary. It's just here
%because CO2 data fits best with it. I tried to see if I could make this
%code work without it so lets see what that looks like.
%%
startup
%% load dataMatrix
cd('C:\Users\yjd21365\OneDrive - Science and Technology Facilities Council\Documents\2DIR\SGR\') %set this to data path
load('0P100I_dry_cjk2019.mat') %will load in as data, crop_data
%% NEW NEW NEW
%% start CO2
clear dyn
dyn = aDArrayBnd;%first time call empty, later can call with indices

%Orientation relaxation
label = 'Orientation_rel';
oPara = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o2)).*(1+4/5.*exp(-t2./p.tau_o2)).*exp(-T3./(3*p.tau_o2));
oPerp = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o2)).*(1-2/5.*exp(-t2./p.tau_o2)).*exp(-T3./(3*p.tau_o2));
oCrossed = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o2)).*exp(-t2./p.tau_o2).*exp(-T3./(3*p.tau_o2));
oMagic = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o2)).*exp(-T3./(3*p.tau_o2));%need to checkprefactors
fun_array = {oPara};
ind_array = {[1:32]};
dynParams2 = fitParamBnd('tau_o2',23,1,150,'');
dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams2,'free');%change this to dyn(2) when you put kinetics back in

%% SO if you're not going to fit for kinetics its best to turn *OFF* the diagrams relating to anything
%other than the main bands.

%cheat to turn some diagrams off and on

 label = 'cheat';
 on = @(T1,t2,T3,p) 1;
 off = @(T1,t2,T3,p) 0;
 fun_array = {off,on};
 ind_array = {[20 23 25 26 27 28 29 30 31 32 8 11 7 10  9 12 14 17 13 16 15 18 19 22 21 24  ...
             ],[ 1 4 2 5 3 6 ]};
 dynParams3 = fitParamBnd('onoff',0,1,0,'fixed');
 dyn(2) = aDArrayBnd(label,fun_array,ind_array,dynParams3,'');
%%
clear lsfParams
%singleexp
lsfParams = fitParamBnd('Delta1_cm',1.62,0.1,5,''); 
lsfParams(2) = fitParamBnd('tau1',15.8,5,250,'');
lsfParams(3) = fitParamBnd('T2',3.13,0.1,10,'');
lsf = lsf1exp1fastBnd(lsfParams,'free');


%
scans = [1:length(crop_data)];
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2342.9,2338,2348,'fixed');
aRFoptions.dw_sb_cm = fitParamBnd('dw_sb_cm',12.2,10,15,'fixed');
%aRFoptions.mu01sq = fitParamBnd('mu01sq',6.21,0,25,'fixed');x
aRFoptions.mu01sq = fitParamBnd('mu01sq',0.071,0.0000001,1.5,'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.82,0,2.5,'fixed');
aRFoptions.anh_cm = fitParamBnd('anh_cm',23.5,0,30,'fixed');
aRFoptions.phase_deg = fitParamBnd('phase_deg',-1.1,-90,90,'fixed');
aRFoptions.damping = lsf;
aRFoptions.dyn = dyn;
aRFoptions.dE_cm = fitParamBnd('dE_cm',667,650,700,'fixed');
aRFoptions.temperature = fitParamBnd('temperature',295,290,505,'fixed');
aRFoptions.t2_array = [crop_data(scans).t2]./1000;
aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
aRFoptions.w3_in = (crop_data(1).w3);%2296:2:2350;
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;



s = aRFCO2Bnd(aRFoptions)

%
s = s.calcSpectrum(s.p0)

s.dataMatrix = prepareGlobalFitData(crop_data(scans))/(1e4); 

%
nscans = zeros(1,length(scans));
for ii = 1:length(scans)
    nscans(ii) = crop_data(scans(ii)).PARAMS.nScans;
end

%
wt = nScansToWeights(s(1).dataMatrix,nscans);
s.weightMatrix = wt;

 %% chi by eye

% compare across all of them
close all
matrixPlot2DIR(s.dataMatrix, s.w1_in,s.w3_in,s.t2_array,[3 9],'zlimit',1) %data
matrixPlot2DIR(s.simMatrix,s.w1_in,s.w3_in,s.t2_array,[3 9],'zlimit',1) %sim
matrixPlot2DIR(s.simMatrix-s.dataMatrix,s.w1_in,s.w3_in,s.t2_array,[3 9]) %residuals
%% lets see if the fitting algorithm can get us closer
s = s.globalFitW;
%%
s = s.calcSpectrum(s.pfit); %now re-simulate the data with the 'new' numbers

%% So lets turn those diagrams back on so we can actually ... do things.

%% case 3 - mixed filling of |10>
%four level system
%ku = 00-10 and 01-11
%kd = 10-00 and 11-01
%k_hgs10 = 01-10
%k_hgs00 = 01-00
%k_hgs2 = 11-10
% 
syms k_u k_d k_hgs10 k_hgs00 t


K = [-k_u,      k_hgs00,  k_d,           0;
        0, -k_hgs00-k_hgs10-k_u,    0,         k_d;
      k_u,           k_hgs10, -k_d,      k_hgs00;
        0,         k_u,    0, -k_hgs00-k_d];

[V3, D3] = eig(K);
ev3 = diag(D3);


V3_0 = [1; 0; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

A =  matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
B = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
C = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
D = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);

V3_0 = [0; 1; 0; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

E =  matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
F = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
G = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
H = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);

V3_0 = [0; 0; 1; 0];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

I = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
J = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
K = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
L = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);

V3_0 = [0; 0; 0; 1];

p3 = V3*expm(D3*t)*V3^(-1);
V3_t = p3*V3_0;

M = matlabFunction(V3_t(1),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
N = matlabFunction(V3_t(2),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
O = matlabFunction(V3_t(3),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);
P = matlabFunction(V3_t(4),'Vars',[k_u,k_d,k_hgs10,k_hgs00,t]);

%%
AA = @(T1,t2,T3,p) A(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
BB = @(T1,t2,T3,p) B(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
CC = @(T1,t2,T3,p) C(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
DD = @(T1,t2,T3,p) D(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
EE = @(T1,t2,T3,p) E(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
FF = @(T1,t2,T3,p) F(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
GG = @(T1,t2,T3,p) G(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
HH = @(T1,t2,T3,p) H(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
II = @(T1,t2,T3,p) I(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
JJ = @(T1,t2,T3,p) J(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
KK = @(T1,t2,T3,p) K(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
LL = @(T1,t2,T3,p) L(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
MM = @(T1,t2,T3,p) M(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
NN = @(T1,t2,T3,p) N(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
OO = @(T1,t2,T3,p) O(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);
PP = @(T1,t2,T3,p) P(p.k_u,p.k_u/(2*exp(-p.dE_cm/(0.69503476*p.temperature))),p.k_hgs10,p.k_hgs00,t2);

%% load dataMatrix
cd('C:\Users\yjd21365\OneDrive - Science and Technology Facilities Council\Documents\2DIR\SGR\') %MAC (for CJ's data)
load('0P100I_dry_cjk2019.mat') %will load in as data, crop_data
%% NEW NEW NEW
%% start CO2
clear dyn
dyn = aDArrayBnd;%first time call empty, later can call with indices
label = 'kinetics';
fun_array = {AA BB CC DD ... 
             EE FF GG HH ...
             II JJ KK LL...
             MM NN OO PP};
         
ind_array = {[1 4], [], [13 16], [],...
            [27 28], [2 5 3 6], [25 26], [14 17 15 18],...
            [19 22], [], [7 10], [],...
            [31 32], [20 23 21 24], [29 30], [8 11 9 12]};
dynParams1    = fitParamBnd('k_u',0.0105,0.00001,0.1,'');
dynParams1(2) = fitParamBnd('k_hgs10',0.0222,1e-10,0.1,'');
dynParams1(3) = fitParamBnd('k_hgs00',0.00008,1e-10,0.5,'');
dyn(1) = aDArrayBnd(label,fun_array,ind_array,dynParams1,'fixed');

%Orientation relaxation
label = 'Orientation_rel';
oPara = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o2)).*(1+4/5.*exp(-t2./p.tau_o2)).*exp(-T3./(3*p.tau_o2));
oPerp = @(T1,t2,T3,p) (1/9).*exp(-T1./(3.*p.tau_o2)).*(1-2/5.*exp(-t2./p.tau_o2)).*exp(-T3./(3*p.tau_o2));
oCrossed = @(T1,t2,T3,p) 1/15.*exp(-T1./(3.*p.tau_o2)).*exp(-t2./p.tau_o2).*exp(-T3./(3*p.tau_o2));
oMagic = @(T1,t2,T3,p) 1/9.*exp(-T1./(3.*p.tau_o2)).*exp(-T3./(3*p.tau_o2));%need to checkprefactors
fun_array = {oPara};
ind_array = {[1:32]};
dynParams2 = fitParamBnd('tau_o2',38,1,150,'');
dyn(2) = aDArrayBnd(label,fun_array,ind_array,dynParams2,'free');%change this to dyn(2) when you put kinetics back in

%%
clear lsfParams
% %singleexp
lsfParams = fitParamBnd('Delta1_cm',1.62,0.1,5,''); %%%check how calcd!!!
lsfParams(2) = fitParamBnd('tau1',15.8,5,250,'');
%lsfParams (3) = fitParamBnd('Delta2_cm',1.94,0.1,5,''); %%%check how calcd!!!
%lsfParams(4) = fitParamBnd('tau2',179,10,250,'');
lsfParams(3) = fitParamBnd('T2',3.13,0.1,10,'');
lsf = lsf1exp1fastBnd(lsfParams,'free');

%
scans = [1:length(crop_data)];
clear aRFoptions s;
aRFoptions.dt = 0.2;
aRFoptions.n_t = 64;
aRFoptions.n_zp = 2*aRFoptions.n_t;
aRFoptions.flag_rotating_frame = true;
aRFoptions.w_01_cm = fitParamBnd('w_01_cm',2342.9,2338,2348,'fixed');
aRFoptions.dw_sb_cm = fitParamBnd('dw_sb_cm',12.2,10,15,'fixed');
%aRFoptions.mu01sq = fitParamBnd('mu01sq',6.21,0,25,'fixed');x
aRFoptions.mu01sq = fitParamBnd('mu01sq',0.071,0.0000001,1.5,'fixed');
aRFoptions.mu12sqRatio = fitParamBnd('mu12sqRatio',1.82,0,2.5,'fixed');
aRFoptions.anh_cm = fitParamBnd('anh_cm',23.5,0,30,'fixed');
aRFoptions.phase_deg = fitParamBnd('phase_deg',-1.1,-90,90,'fixed');
aRFoptions.damping = lsf;
aRFoptions.dyn = dyn;
aRFoptions.dE_cm = fitParamBnd('dE_cm',667,650,700,'fixed');
aRFoptions.temperature = fitParamBnd('temperature',295,290,505,'fixed');
aRFoptions.t2_array = [crop_data(scans).t2]./1000;
aRFoptions.w1_in = crop_data(1).w1; %2320:0.5:2350;
aRFoptions.w3_in = (crop_data(1).w3);%2296:2:2350;
aRFoptions.useParallel=false;  
aRFoptions.nboot = 100;



s = aRFCO2Bnd(aRFoptions)

%
s = s.calcSpectrum(s.p0)


%s.dataMatrix = prepareGlobalFitData(PEGDA_2d_crop); 
s.dataMatrix = prepareGlobalFitData(crop_data(scans))/(1e4); 

%
nscans = zeros(1,length(scans));
for ii = 1:length(scans)
    nscans(ii) = crop_data(scans(ii)).PARAMS.nScans;
end

%
wt = nScansToWeights(s(1).dataMatrix,nscans);
s.weightMatrix = wt;

 %% chi by eye

% compare across all of them
close all
matrixPlot2DIR(s.dataMatrix, s.w1_in,s.w3_in,s.t2_array,[3 9],'zlimit',1) %data
matrixPlot2DIR(s.simMatrix,s.w1_in,s.w3_in,s.t2_array,[3 9],'zlimit',1) %sim
matrixPlot2DIR(s.simMatrix-s.dataMatrix,s.w1_in,s.w3_in,s.t2_array,[3 9]) %residuals
%% lets see if the fitting algorithm can get us closer
s = s.globalFitW;
%%
s = s.calcSpectrum(s.pfit); %now re-simulate the data with the 'new' numbers


