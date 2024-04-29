%% Alcott, Mills and Poulton 2019, Science
% Model based on Slomp and VanCappellen, 2007, Biogeosciences; Tsandev et al., 2008, GBC; Tsandev and Slomp, 2009, EPSL.
% Model equations (do not run this script)

function dy = Alcott_et_al_2019_Science(t,y)
% Set up dy array
dy = zeros(21,1) ;

%%% Set up Global parameters
global stepnumber
global pars
global workingstate
global starting
global per
global present
global forcings

%% Reservoirs
Water_P = y(1) ;
Water_D = y(2) ;
Water_S = y(3) ;
Water_DP = y(4) ;
POC_P = y(5) ;
POC_D = y(6) ;
POC_S = y(7) ;
POC_DP = y(8) ;
O2_DP = y(12) ;
SRP_P = y(13) ;
OP_P = y(14) ;
SRP_D = y(15) ;
OP_D = y(16) ;
SRP_S = y(17) ;
OP_S = y(18) ;
SRP_DP = y(19) ;
OP_DP = y(20) ;
O2_A = y(21) ;

%% Forcings

%%%% P weathering multiplier
Pforce = per.P ;

% %%% extra P input from weathering
% if isempty(forcings.Weather_P) == 0
%     P_WEATHER = interp1q( forcings.Weather_P(:,11).*1e6 , forcings.Weather_P(:,forcings.pchoice) , t ) ;
% else
%     P_WEATHER = 0 ;
% end

%%% extra P input from weathering
if isempty(forcings.Weather_P) == 0
    P_WEATHER = interp1( forcings.Time_P , forcings.Weather_P(forcings.pchoice,:) , t/1e6 ) ;
else
    P_WEATHER = 0 ;
end

%%% P input from ridges
if isempty(forcings.Ridge_P) == 0
    P_RIDGE = interp1( forcings.Time_P , forcings.Ridge_P(forcings.pchoice,:) , t/1e6 ) ;
else
    P_RIDGE = 0 ;
end
    
generic_reductant_flux_i = 0 ;
% if t > 0
%     generic_reductant_flux_i = 0 ;
% else
%     generic_reductant_flux_i = interp1([-4e9 0], [45e12 0],t)  ;
% end

%% Concentrations 
present.Conc_O2_deep = present.O2_DP / starting.Water_DP ;
O2_DPconc    = y(12)/y(4);                        
OP_Dconc  = y(16)/y(2);                           
OP_Pconc  = y(14)/y(1);                           
SRP_DPconc = y(19)/y(4);                           
SRP_Dconc = y(15)/y(2);                            
SRP_Pconc = y(13)/y(1);                            
SRP_Sconc = y(17)/y(3);                            
    

%% Oceanic Water Cycle
%%% As in Slomp and Van Cappellen, 2007

%River flow entering Water_P
River_Water = 37e12 ; %Berner and Berner1996

%flow from Water_P to Water_D
Water_P_D = River_Water ; 

%Coastal Upwelling
% Water_DP_D = pars.kWF6 * Water_DP ; %12Sv -  Brink et al1995
Water_DP_D = pars.Water_DP_D_1 ;
%flow from Water_D to Water_S
% Water_D_S = Water_P_D + Water_DP_D ;
Water_D_S = pars.Water_D_S_2 ;
%Oceanic Upwelling
% Water_DP_S = Water_DP* pars.kWF5 * pars.vmix ; %120Sv - Brink et al1995
Water_DP_S = pars.Deep_Water_S_0 * pars.vmix ;

%downwelling (Water_S to Water_DP)
% Water_S_DP = pars.kWF4 * Water_S * pars.vmix ;
% Water_S_DP = ( pars.kWF4 * Water_S * pars.vmix ) + pars.Water_S_DP_diff ;
Water_S_DP = pars.Water_S_DP_0 * pars.vmix ;
%Evaporation from Water_S
Evaporation_Water = River_Water ;






%% Hydrological Cycle

%Proximal Zone
dy(1) = River_Water - Water_P_D ;

%Distal Zone
dy(2) = Water_P_D + Water_DP_D - Water_D_S ;

%Surface Ocean
dy(3) = Water_DP_S - Water_S_DP + Water_D_S - Evaporation_Water ;
% dy(3) = Water_DP_S  + Water_D_S - Evaporation_Water ;

%Deep Ocean
dy(4) = -Water_DP_S + Water_S_DP - Water_DP_D  ;
% dy(4) = Water_DP_S - Water_DP_D  ;


%% Normalised to present day reservoirs
Norm_SRP_P = SRP_P / present.SRP_P ;
Norm_SRP_D = SRP_D / present.SRP_D ;
Norm_SRP_S = SRP_S / present.SRP_S ;
Norm_OP_P = OP_P / present.OP_P ;
Norm_OP_D = OP_D / present.OP_D ;
Norm_OP_S = OP_S / present.OP_S ;
Norm_POC_P = POC_P / present.POC_P ;
Norm_POC_D = POC_D / present.POC_D ;
Norm_POC_S = POC_S / present.POC_S ;
Norm_POC_DP = POC_DP / present.POC_DP ;
Norm_O2_A = O2_A / present.O2_A ;


%% Marine Carbon Cycle

% Primary production in Proximal
PP_P = pars.kPhotoprox * Norm_SRP_P * pars.Redfield_CP ; 

% POC mineralisation in Proximal
POC_Min_P = pars.kminprox * Norm_POC_P ;

% POC export from Proximal to Distal
OP_P_D = Water_P_D * OP_Pconc ;
XP_P_D = OP_P_D * pars.Redfield_CP ; 

% Proximal sediment POC burial
POC_P_Burial = pars.Prox_C_Bur * PP_P ; 

% Primary Production in Distal
PP_D = pars.kPhotodist * Norm_SRP_D * pars.Redfield_CP ; 

% POC mineralisation in Distal
POC_Min_D = pars.kmindist * Norm_POC_D ; 

% POC Export from Distal to Surface
OP_D_S = Water_D_S * OP_Dconc ;
XP_D_S = OP_D_S * pars.Redfield_CP ; 

% Distal sediment POC burial
POC_D_Burial = pars.Dist_C_Bur * ( XP_P_D + PP_D ) ; 

% Primary Production in Surface
PP_S = pars.kPhotosurf * Norm_SRP_S * pars.Redfield_CP ; 

% POC mineralisation in Surface
POC_Min_S = pars.kminsurf * Norm_POC_S ; 

% POC export from Surface to Deep
XP_S_DP = pars.Surf_Deep_XP * ( XP_D_S + PP_S ) ; 

% POC respiration in Water_DP
POC_DP_Resp = pars.kCF12 * Norm_POC_DP ; 


%% Deep Carbon Burial

% C:P ratios for oxic and anoxic conditions
pars.CPoxic = 250 ;
pars.CPanoxic_prox = 1100 ;
pars.CPanoxic_dist = 1100 ;
pars.CPanoxic_deep = 1100  ;

if O2_DPconc < present.Conc_O2_deep
    
    % Deep sediment POP burial
    % Redox dependancy is set up according to Slomp and VC, 2007 and Tsandev et al., 2009. 
    OP_DP_Burial = ( pars.kPOP_Bur_Deep * XP_S_DP / pars.CPoxic ) * ( (1-per.POP_deep_feedback) + (per.POP_deep_feedback * O2_DPconc / present.Conc_O2_deep ));
    
    % Deep sediment FeP burial
    P_FeP_DP = pars.kFeP_Deep * ( O2_DPconc / present.Conc_O2_deep );
    
    % Deep sediment POC burial
    POC_DP_Burial = OP_DP_Burial * ( ( pars.CPanoxic_deep * pars.CPoxic ) / ( ( ( O2_DPconc / present.Conc_O2_deep ) * pars.CPanoxic_deep ) + ( ( 1 - ( O2_DPconc / present.Conc_O2_deep ) ) * pars.CPoxic ) ) ) ;
        
else
     OP_DP_Burial= pars.kPOP_Bur_Deep * XP_S_DP / 250;  
     P_FeP_DP = pars.kFeP_Deep ;
     POC_DP_Burial = 250 * OP_DP_Burial ;    
end

%Carbon Proximal Zone
dy(5) = PP_P - POC_Min_P - POC_P_Burial  - XP_P_D ;

%Carbon Distal Zone
dy(6) = XP_P_D + PP_D - POC_D_Burial - POC_Min_D - XP_D_S ; 

%Carbon Surface Ocean
dy(7) = XP_D_S + PP_S - POC_Min_S - XP_S_DP ;

%Carbon Deep Ocean
dy(8) = XP_S_DP - POC_DP_Resp - POC_DP_Burial ;


%% Oxygen Cycle

% Linear relationship for oxygen content in boxes in contact with atmosphere 
O2_P = ( present.O2_P * ( O2_A / present.O2_A ) ) ;
O2_S = ( present.O2_S * ( O2_A / present.O2_A ) ) ;
O2_D = ( present.O2_D * ( O2_A / present.O2_A ) );
Norm_O2_D = O2_D / present.O2_D ;

% fanoxic parameters (From Watson et al., 2017)
kanox = 10 ; 
O2O20 = Norm_O2_A ;
kU = 0.4 ;

% fanoxic calculation from Watson et al., 2017
fanoxicdist = 1 / ( 1 + exp(-kanox * ( kU * Norm_SRP_D - O2O20 ) ) ) ; 
fanoxicprox = 1 / ( 1 + exp(-kanox * ( kU * Norm_SRP_P - O2O20 ) ) ) ; 

present.fanoxicprox = 0.0025 ;
present.fanoxicdist = 0.0025 ;

%Concentration of oxygen 
O2_Sconc = O2_S/Water_S ;

%O2 Downwelling
O2_S_DP = Water_S_DP * O2_Sconc ;  

%O2 coastal upwelling
O2_DP_D =  Water_DP_D * O2_DPconc; 

%O2 oceanic upwelling
O2_DP_S = Water_DP_S * O2_DPconc; 

%Aerobic O2 respiration
Respiration_O21 = pars.kCF12 * Norm_POC_DP / pars.Redfield_CO2; % Preliminary Respiration quanitity before Monod inclusion
KmO2 = 0.0001 ; % monod constant for oxic respiration in mol/m3 
Mon_O2_deep = O2_DPconc / ( KmO2 + O2_DPconc ) ;
Respiration_O2 = Respiration_O21 * Mon_O2_deep ;

%% Oxygen dys

%Oxygen Deep Ocean 
dy(12) = - Respiration_O2 - O2_DP_S - O2_DP_D + O2_S_DP ;

Atmos_Weather = pars.O2_A_Weathering * (sqrt(O2_A/ present.O2_A)) ;

generic_reductant_flux = generic_reductant_flux_i * sigmf((log10(Norm_O2_A)),[3,-5]) ; 

Total_POC_Burial = POC_P_Burial + POC_D_Burial + POC_DP_Burial ; 

locb = starting.locb ;


%% Oxygen atmosphere
dy(21) = Total_POC_Burial - Atmos_Weather  - generic_reductant_flux + locb ;
% dy(21) = 0 ;

%% P Cycle

% SRP upwelling Water_DP to Water_S
SRP_DP_S = SRP_DPconc * Water_DP_S ; 

% SRP upwelling Water_DP to Water_D
SRP_DP_D = Water_DP_D * SRP_DPconc ; 


%% Proximal Coastal
%O2 limit for scavenging
scavlim = 1e-3 ;

% Forcing to vary Riverine P input.
River_SRP = pars.River_SRP_0 * Pforce ;

% Primary Production proximal
P_PP_P =  PP_P/pars.Redfield_CP ; 

% POP mineralisation
OP_P_Min = Norm_OP_P * pars.kPrel_prox ; 

% SRP transport from Proximal to Distal
SRP_P_D = SRP_Pconc * Water_P_D ;

% Proximal sediment POP burial 
OP_P_Burial = pars.Prox_C_Bur * PP_P * ( ( ( 1-fanoxicprox ) / pars.CPoxic ) + ( fanoxicprox / pars.CPanoxic_prox ) ) ;

% Proximal FeP burial
P_FeP_P = pars.kFePprox * SRP_P ;

% Proximal CaP burial
P_AuthP_P = pars.kPrel_prox * OP_P * pars.kCaP_prox ; 


%% Distal Coastal

% Primary Production in Water_D
P_PP_D = PP_D / pars.Redfield_CP ; 

% POP mineralisation in Water_D
OP_D_Min = Norm_OP_D * pars.kPrel_dist ; 

% SRP transport Water_D to Water_S
SRP_D_S = SRP_Dconc * Water_D_S ;

% Distal sediment POP burial
OP_D_Burial = pars.kPOPDOADist * (PP_D + XP_P_D) * ( ( ( 1-fanoxicdist ) / pars.CPoxic ) + ( fanoxicdist / pars.CPanoxic_dist ) ) ;

% Distal sediment FeP burial
P_FeP_D = pars.kFePDOADist * SRP_D * ( 1- fanoxicdist ) ; 

% Distal Sediment CaP burial
P_AuthP_D = pars.kCaPDOADist * Norm_O2_D * ( 1 - fanoxicdist ) ;

%% Surface Ocean

% Primary Production 
P_PP_S = PP_S / pars.Redfield_CP; % Links P to C cycle

% POP mineralisation 
OP_S_Min = Norm_OP_S * pars.kPrel_surf ; 

% SRP downwelling from Water_S to Water_DP
SRP_S_DP = SRP_Sconc * Water_S_DP ;

% POP export from Water_S to Water_DP
OP_S_DP = pars.kCF11 * ( PP_S + XP_D_S ) / pars.Redfield_CP ; % Links P to C cycle

% Scavenging flux of FeP from surface ocean.
eSCAV_SURF = 1 / ( 1 + exp(5000*(O2_DPconc - scavlim ))) ; % Simple logistic curve, similar to Reinhard et al., 2017
Fe_SCAV_SURF = min( eSCAV_SURF * SRP_DP_S * per.sig_SCAV , 4.5e10 ); % upper limit variable

%% Deep Ocean

% POP mineralisation in Water_DP
OP_DP_Min = pars.kPrel_deep * OP_DP ;

% Deep sediment CaP burial
% Redox dependancy is set up according to Slomp and VC, 2007 and Tsandev et al., 2009. 
P_AuthP_DP =  pars.fPF34 * OP_DP_Min * ( (1-per.CaP_deep_feedback) + ( per.CaP_deep_feedback * ( O2_DPconc/present.Conc_O2_deep ) ) ) ; 



%% Phosphorus Differentials

    %SRP Proximal
    dy(13) = River_SRP - P_PP_P + OP_P_Min - P_FeP_P - P_AuthP_P - SRP_P_D + P_WEATHER ; 

    %POP Proximal
    dy(14) = P_PP_P - OP_P_Min - OP_P_Burial - OP_P_D  ;

    %SRP Distal
    dy(15) = SRP_P_D - P_PP_D + OP_D_Min - P_FeP_D - P_AuthP_D - SRP_D_S + SRP_DP_D ;

    %POP Distal
    dy(16) = OP_P_D + P_PP_D - OP_D_Min - OP_D_Burial - OP_D_S ;

    %SRP Surface Ocean
    dy(17) = SRP_D_S - P_PP_S + OP_S_Min - SRP_S_DP + SRP_DP_S ;

    %POP Surface Ocean
    dy(18) = OP_D_S + P_PP_S - OP_S_Min - OP_S_DP;

    %SRP Deep Ocean       
    dy(19) = SRP_S_DP + OP_DP_Min - P_FeP_DP - P_AuthP_DP - SRP_DP_S - SRP_DP_D - Fe_SCAV_SURF + P_RIDGE ;

    %POP Deep Ocean
    dy(20) = OP_S_DP - OP_DP_Min - OP_DP_Burial ;


%% Saving data
workingstate.generic_reductant_flux(stepnumber,1) = generic_reductant_flux ;
workingstate.Water_P(stepnumber,1) = Water_P ;
workingstate.Water_D(stepnumber,1) = Water_D ;
workingstate.Water_S(stepnumber,1) = Water_S ;
workingstate.Water_DP(stepnumber,1) = Water_DP;
workingstate.River_Water(stepnumber,1) = River_Water ;
workingstate.Prox_Water_D(stepnumber,1) = Water_P_D ;
workingstate.Dist_Water_S(stepnumber,1) = Water_D_S ;
workingstate.Surf_Water_DP(stepnumber,1) = Water_S_DP ;
workingstate.Water_DP_S(stepnumber,1) = Water_DP_S ;
workingstate.Water_DP_D(stepnumber,1) = Water_DP_D ;
workingstate.Evaporation_Water(stepnumber,1) = Evaporation_Water ;
workingstate.POC_P(stepnumber,1) = POC_P ;
workingstate.POC_D(stepnumber,1) = POC_D ;
workingstate.POC_S(stepnumber,1) = POC_S ;
workingstate.POC_DP(stepnumber,1) = POC_DP ;
workingstate.PP_P(stepnumber,1) = PP_P ;
workingstate.POC_Min_P(stepnumber,1) = POC_Min_P ;
workingstate.POC_P_Burial(stepnumber,1) = POC_P_Burial ;
workingstate.XP_P_D(stepnumber,1) = XP_P_D ;
workingstate.PP_D(stepnumber,1) = PP_D ;
workingstate.POC_Min_D(stepnumber,1) = POC_Min_D ;
workingstate.POC_D_Burial(stepnumber,1) = POC_D_Burial ;
workingstate.XP_D_S(stepnumber,1) = XP_D_S ;
workingstate.PP_S(stepnumber,1) = PP_S ;
workingstate.POC_Min_S(stepnumber,1) = POC_Min_S ;
workingstate.XP_S_DP(stepnumber,1) = XP_S_DP ;
workingstate.POC_DP_Resp(stepnumber,1) = POC_DP_Resp ;
workingstate.POC_DP_Burial(stepnumber,1) = POC_DP_Burial ;
workingstate.O2_P(stepnumber,1) = O2_P ;
workingstate.O2_D(stepnumber,1) = O2_D ;
workingstate.O2_S(stepnumber,1) = O2_S ;
workingstate.O2_DP(stepnumber,1) = O2_DP ;
workingstate.O2_A(stepnumber,1) = O2_A ;
workingstate.O2_S_DP(stepnumber,1) = O2_S_DP ;
workingstate.Respiration_O2(stepnumber,1) = Respiration_O2 ;
workingstate.O2_DP_D(stepnumber,1) = O2_DP_D ;
workingstate.O2_DP_S(stepnumber,1) = O2_DP_S ;
workingstate.SRP_P(stepnumber,1) = SRP_P ;
workingstate.OP_P(stepnumber,1) = OP_P ;
workingstate.SRP_D(stepnumber,1) = SRP_D ;
workingstate.OP_D(stepnumber,1) = OP_Dconc ;
workingstate.SRP_S(stepnumber,1) = SRP_S ;
workingstate.OP_S(stepnumber,1) = OP_S ;
workingstate.SRP_DP(stepnumber,1) = SRP_DP ;
workingstate.OP_DP(stepnumber,1) = OP_DP ;
workingstate.River_SRP(stepnumber,1) = River_SRP ;
workingstate.P_PP_P(stepnumber,1) = P_PP_P ;
workingstate.OP_P_Min(stepnumber,1) = OP_P_Min ;
workingstate.OP_P_Burial(stepnumber,1) = OP_P_Burial ;
workingstate.P_FeP_P(stepnumber,1) = P_FeP_P ;
workingstate.P_AuthP_P(stepnumber,1) = P_AuthP_P ;
workingstate.SRP_P_D(stepnumber,1) = SRP_P_D ;
workingstate.OP_P_D(stepnumber,1) = OP_P_D ;
workingstate.P_PP_D(stepnumber,1) = P_PP_D ;
workingstate.OP_D_Min(stepnumber,1) = OP_D_Min ;
workingstate.OP_D_Burial(stepnumber,1) = OP_D_Burial ;
workingstate.P_FeP_D(stepnumber,1) = P_FeP_D ;
workingstate.P_AuthP_D(stepnumber,1) = P_AuthP_D ;
workingstate.SRP_D_S(stepnumber,1) = SRP_D_S ;
workingstate.OP_D_S(stepnumber,1) = OP_D_S ;
workingstate.P_PP_S(stepnumber,1) = P_PP_S ;
workingstate.OP_S_Min(stepnumber,1) = OP_S_Min ;
workingstate.SRP_S_DP(stepnumber,1) = SRP_S_DP ;
workingstate.OP_S_DP(stepnumber,1) = OP_S_DP ;
workingstate.OP_DP_Min(stepnumber,1) = OP_DP_Min ;
workingstate.OP_DP_Burial(stepnumber,1) = OP_DP_Burial ;
workingstate.P_FeP_DP(stepnumber,1) = P_FeP_DP ;
workingstate.P_AuthP_DP(stepnumber,1) = P_AuthP_DP ;
workingstate.SRP_DP_S(stepnumber,1) = SRP_DP_S ;
workingstate.SRP_DP_D(stepnumber,1) = SRP_DP_D ;
workingstate.CP_Dist(stepnumber,1) = ( ( ( 1-fanoxicdist ) * pars.CPoxic ) + ( fanoxicdist * pars.CPanoxic_dist ) ) ;
workingstate.CP_Deep(stepnumber,1) = POC_DP_Burial / OP_DP_Burial ;
workingstate.CP_Prox(stepnumber,1) = ( ( ( 1-fanoxicprox ) * pars.CPoxic ) + ( fanoxicprox * pars.CPanoxic_prox ) ) ;
workingstate.Deep_Preac_Burial(stepnumber,1) = P_AuthP_DP + P_FeP_DP + OP_DP_Burial ;
workingstate.Dist_Preac_Burial(stepnumber,1) = P_AuthP_D + P_FeP_D + OP_D_Burial ;
workingstate.fanoxicdist(stepnumber,1) = fanoxicdist ;
workingstate.fanoxicprox(stepnumber,1) = fanoxicprox ;
workingstate.Atmos_Weather(stepnumber,1) = Atmos_Weather ;
workingstate.generic_reductant_flux(stepnumber,1) = generic_reductant_flux ;
workingstate.Fe_SCAV_SURF(stepnumber,1)= Fe_SCAV_SURF ;
workingstate.P_RIDGE(stepnumber,1)= P_RIDGE ;
workingstate.P_WEATHER(stepnumber,1)= P_WEATHER ;

workingstate.O2_DPconc(stepnumber,1) = O2_DPconc ;

%%%%%%%% record time
workingstate.time(stepnumber,1) = t ;

%%%% final action: record current model step
stepnumber = stepnumber + 1 ;



end
