function MWS = setup_ShaftGB_PwrSys(MWS)

% setup variables for power system and shaft/gearbox. Creates substructures
% PowerSystem and Shaft.

% Inputs DEM and VEATE define parameters for the dedicated EM and VEATE
% approaches respectively

% Pick of and re-assign key parameters for convenience -------------------%

HybConfig = MWS.In.Options.HybridConfig;
GBConfig = MWS.In.Options.PGBConfig;

% Load in needed data ----------------------------------------------------%

% speed an power data to help size the motors and gearbox
%   NOTE: chose data at 5000ft and Mach 0.2 because it same to more or less
%   result in the most constraining weights.
ShaftData = load('Data_5kft02MN.mat');
if HybConfig == 1
    Nlps = ShaftData.STD.Nlps;
    Nhps = ShaftData.STD.Nhps;
    Nlps_EPT = ShaftData.STD.Nlps_EPT;
    Nhps_EPT = ShaftData.STD.Nhps_EPT;
    Plps_nom = ShaftData.STD.Plps_nom;
    Phps_nom = ShaftData.STD.Phps_nom;
    Plps_nom_EPT = ShaftData.STD.Plps_nom_EPT;
    Phps_nom_EPT = ShaftData.STD.Phps_nom_EPT;
elseif HybConfig == 2
    Nlps = ShaftData.Boost.Nlps;
    Nhps = ShaftData.Boost.Nhps;
    Nlps_EPT = ShaftData.Boost.Nlps_EPT;
    Nhps_EPT = ShaftData.Boost.Nhps_EPT;
    Plps_nom = ShaftData.Boost.Plps_nom;
    Phps_nom = ShaftData.Boost.Phps_nom;
    Plps_nom_EPT = ShaftData.Boost.Plps_nom_EPT;
    Phps_nom_EPT = ShaftData.Boost.Phps_nom_EPT;
else
    Nlps = ShaftData.PEx.Nlps;
    Nhps = ShaftData.PEx.Nhps;
    Nlps_EPT = ShaftData.PEx.Nlps_EPT;
    Nhps_EPT = ShaftData.PEx.Nhps_EPT;
    Plps_nom = ShaftData.PEx.Plps_nom;
    Phps_nom = ShaftData.PEx.Phps_nom;
    Plps_nom_EPT = ShaftData.PEx.Plps_nom_EPT;
    Phps_nom_EPT = ShaftData.PEx.Phps_nom_EPT;
end

% General Parameters -----------------------------------------------------%

% Engine Shafts
% -- Low pressure shaft inertia, slugs*ft2
% NOTE: LPS inertia defined in PGBoptFun.m --> J_LPS = ((450102/(3.1*3.1)) + 17796 + 16237) /32.2/144.; 
MWS.Shaft.LPS_Eff = 0.99;
% -- High pressure shaft inertia, slugs*ft2
% NOTE: HPS inertia defined in PGBoptFun.m --> J_HPS = 8627/32.2/144.;

% Design parameters - applies to DEM and VEATE
Nspool_max = [6831 23549]; %max speed of LPS and HPS, rpm
PEx_ACsys = 175;
PX_EPT = 250; %EPT power transfer (LPS --> HPS), hp
Pin_TEEMa = 500; %TEEM power input to HPS during accels, hp
PExLMax_TEEMa = 200; %Maximum power extraction from LPS during TEEM accel (prevent LPS speed dip), hp
PX_TEEMd = 200; %TEEM power transfer (LPS --> HPS) during decels, hp
PmaxPerEM = 5000; %Maximum power allowed to a single EM
energyDenESD = 0.7297; %energy density of energy storage, for super-cap a good 
    % range is 0.7297 - 1.459 kW-hr/slug (typical - sporty) with 2.652 kW-hr/slug
    % being a stretch
PwrElecPwrDen = 146.78; % power electronics (inverter/rectifier/converter) power 
    % density. Based on "Projecting Power Converter Specific Power Through
    % 2050 for Aerospace Applications", by Hall, Pastra, et. al. the SOA is
    % 7.5kW/kg (146.78hp/slug) and 2050 projection could range from 12 -
    % 52.9kW/kg (234.85 - 1035.29 hp/slug)
transientDur = 5; % transient duration, sec
Psf = 1.2; % power safety factor
Trqsf = 1.2; % torque safety factor
Esf = 1.4; % energy storage safety factor
GBsf = 1.5; % gearbox safety factor
DiIS = 0.25; % inner diameter of inner shaft, in
BearingGap = 0.25; % gap between tower shafts for bearings, in
UseToggle = [1 1 1 1]; %[Nominal EPT TEEMaccel TEEMdecel]

% Design parameters - DEM only
Nem_max_DEM = [30000 30000]; %max speed of electric machines [sun ring carrier], rpm
if HybConfig == 1
    NtsqNhps_DEM = 1.0; %HPS tower shafts speed / HPS speed
    NtsqNlps_DEM = 2.0365; %LPS tower shafts speed / LPS speed
elseif HybConfig == 2
    NtsqNhps_DEM = 1.0; %HPS tower shafts speed / HPS speed
    NtsqNlps_DEM = 1.3621; %LPS tower shafts speed / LPS speed
else
    NtsqNhps_DEM = 1.0; %HPS tower shafts speed / HPS speed
    NtsqNlps_DEM = 1.7976; %LPS tower shafts speed / LPS speed
end

% Design parameters - Planetary Gearbox Only (VEATE only)
rR = 0.375; %ring gear radius, ft
Nem_max_PGB = [30000 30000 30000]; %max speed of electric machines [sun ring carrier], rpm
NPGB_max = [50000 50000 50000 100000]; %maximum speed of sun, ring, carrier, and planets
PemCoup_maxTEEM = 200; %Maximum off-nom power with coupling EM for TEEM, hp
if HybConfig == 3
    PemCoup_maxNom = 0; %Maximum nominal power with coupling EM, hp
else
    PemCoup_maxNom = 400; %Maximum nominal power with coupling EM, hp
end
if HybConfig == 1
    if GBConfig == 1
        rSqrR = 0.5398; %0.443; %radius of sun gear / radius of ring gear
        NLPGBqNem = -1.4088; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.1653; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 2.1009; %LPS tower shafts speed / LPS speed
    elseif GBConfig == 2
        rSqrR = 0.6413; %0.4152; %radius of sun gear / radius of ring gear
        NLPGBqNem = 0.9621; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.3525; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1.9308; %LPS tower shafts speed / LPS speed
    elseif GBConfig == 4
        rSqrR = 0.6098; %0.75; %radius of sun gear / radius of ring gear
        NLPGBqNem = 1.2360; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.1742; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 2.0775; %LPS tower shafts speed / LPS speed
    else
        rSqrR = 0.5; %radius of sun gear / radius of ring gear
        NLPGBqNem = 1; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1; %LPS tower shafts speed / LPS speed
        disp('Warning: No valid solutions have been found for PGB configurations 3, 5, and 6')
    end
elseif HybConfig == 2
    if GBConfig == 1
        rSqrR = 0.5054; %0.384; %radius of sun gear / radius of ring gear
        NLPGBqNem = -2.6855; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.0; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1.6547; %LPS tower shafts speed / LPS speed
    elseif GBConfig == 2
        rSqrR = 0.6483; %0.6319; %radius of sun gear / radius of ring gear
        NLPGBqNem = 2.1176; %1.0145; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.0; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1.6708; %LPS tower shafts speed / LPS speed
    elseif GBConfig == 4
        rSqrR = 0.6089; %radius of sun gear / radius of ring gear
        NLPGBqNem = 2.4875; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.0; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1.7151; %LPS tower shafts speed / LPS speed
    else
        rSqrR = 0.5; %radius of sun gear / radius of ring gear
        NLPGBqNem = 1; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1; %LPS tower shafts speed / LPS speed
        disp('Warning: No valid solutions have been found for PGB configurations 3, 5, and 6')
    end 
else
    if GBConfig == 1
        rSqrR = 0.6142; %0.36; %radius of sun gear / radius of ring gear
        NLPGBqNem = -2.3688; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.0; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1.8302; %LPS tower shafts speed / LPS speed
    elseif GBConfig == 2
        rSqrR = 0.6715; %radius of sun gear / radius of ring gear
        NLPGBqNem = 1.6362; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.0; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1.8580; %LPS tower shafts speed / LPS speed
    elseif GBConfig == 4
        rSqrR = 0.4623; %0.75; %radius of sun gear / radius of ring gear
        NLPGBqNem = 1.7017; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1.0; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1.8530; %LPS tower shafts speed / LPS speed
    else
        rSqrR = 0.5; %radius of sun gear / radius of ring gear
        NLPGBqNem = 1; %speed of PGB gearbox component connected to LPS / speed of EM connected to the LPS
        NtsqNhps_PGB = 1; %HPS tower shafts speed / HPS speed
        NtsqNlps_PGB = 1; %LPS tower shafts speed / LPS speed
        disp('Warning: No valid solutions have been found for PGB configurations 3, 5, and 6')
    end
end
NHPGBqNem = 1; %speed of PGB gearbox component connected to HPS / speed of EM connected to the HPS

% Sizing --------------------------------------------------------------%

% DEM sizing
[Mass_DEM,info_DEM] = DEMSizerFun(HybConfig,Nem_max_DEM,Nspool_max, ...
    Nlps,Nhps,Nhps_EPT,Nlps_EPT,Plps_nom,Phps_nom,Plps_nom_EPT,Phps_nom_EPT, ...
    PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PmaxPerEM, ...
    energyDenESD,PwrElecPwrDen,transientDur,Psf,Trqsf,Esf,GBsf,DiIS,BearingGap, ...
    NtsqNhps_DEM,NtsqNlps_DEM,UseToggle);
% Gear ratios
MWS.Shaft.DEM.LPSGB.GR = info_DEM.NemLqNlps; %GR = Nem/Nspool
MWS.Shaft.DEM.HPSGB.GR = info_DEM.NemHqNhps; %GR = Nem/Nspool
% Power system
% -- LPS EM
MWS.PowerSystem.EML.NMax = info_DEM.NemLrange(end);
MWS.PowerSystem.EML.PMaxC = info_DEM.PemL;
MWS.PowerSystem.EML.PMax = info_DEM.PemL*Psf;
MWS.PowerSystem.EML.TrqMaxC = info_DEM.TrqemL;
MWS.PowerSystem.EML.TrqMax = info_DEM.TrqemL*Trqsf;
MWS.PowerSystem.EML.Mass = info_DEM.MemL;
MWS.PowerSystem.EML.J = info_DEM.JemL;
% -- HPS EM
MWS.PowerSystem.EMH.NMax = info_DEM.NemHrange(end);
MWS.PowerSystem.EMH.PMaxC = info_DEM.PemH;
MWS.PowerSystem.EMH.PMax = info_DEM.PemH*Psf;
MWS.PowerSystem.EMH.TrqMaxC = info_DEM.TrqemH;
MWS.PowerSystem.EMH.TrqMax = info_DEM.TrqemH*Trqsf;
MWS.PowerSystem.EMH.Mass = info_DEM.MemH;
MWS.PowerSystem.EMH.J = info_DEM.JemH;
MWS.PowerSystem.EMH.Mass = info_DEM.MemH;
% -- Power Electronics
MWS.PowerSystem.EMLInverter.Mass = info_DEM.MpeL;
MWS.PowerSystem.EMHInverter.Mass = info_DEM.MpeH;
% -- Effective Shaft Inertias
MWS.Shaft.DEM.LPS_Inertia = info_DEM.JlpsEff;
MWS.Shaft.DEM.HPS_Inertia = info_DEM.JhpsEff;

% PGB (VEATE) Sizing
[Mass_PGB,info_PGB] = PGBoptFun([rSqrR NLPGBqNem NtsqNhps_PGB NtsqNlps_PGB],GBConfig,HybConfig,rR,NHPGBqNem,Nem_max_PGB,Nspool_max, ...
    NPGB_max,Nlps,Nhps,Nhps_EPT,Nlps_EPT,Plps_nom,Phps_nom,Plps_nom_EPT,Phps_nom_EPT, ...
    PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,PmaxPerEM, ...
    energyDenESD,PwrElecPwrDen,transientDur,Psf,Trqsf,Esf,GBsf,DiIS,BearingGap,UseToggle);
% Gearbox Parameters
% -- Sun Gear
MWS.Shaft.PGB.Sun.GR_gb = info_PGB.NSqNemS; %GR = NS/Nint
if GBConfig == 1 || GBConfig == 2 % sun gear connected to HPS
    MWS.Shaft.PGB.Sun.GR_int = info_PGB.NemqNhps; %GR = Nint/Ncomp
elseif GBConfig == 3 || GBConfig == 5 % sun gear connected to LPS
    MWS.Shaft.PGB.Sun.GR_int = info_PGB.NemqNlps; %GR = Nint/Ncomp
else
    MWS.Shaft.PGB.Sun.GR_int = 1;
end
MWS.Shaft.PGB.Sun.r = rSqrR*rR; %radius of sun gear, ft
MWS.Shaft.PGB.Sun.J_go = info_PGB.JSo; %inertia of gear only, slug-ft2
% -- Ring Gear
MWS.Shaft.PGB.Ring.GR_gb = info_PGB.NRqNemR; %GR = NS/Nint
if GBConfig == 3 || GBConfig == 4 % sun gear connected to HPS
    MWS.Shaft.PGB.Ring.GR_int = info_PGB.NemqNhps; %GR = Nint/Ncomp
elseif GBConfig == 1 || GBConfig == 6 % sun gear connected to LPS
    MWS.Shaft.PGB.Ring.GR_int = info_PGB.NemqNlps; %GR = Nint/Ncomp
else
    MWS.Shaft.PGB.Ring.GR_int = 1;
end
MWS.Shaft.PGB.Ring.r = rR; %radius of sun gear, ft
MWS.Shaft.PGB.Ring.J_go = info_PGB.JRo; %inertia of gear only, slug-ft2
% -- Carrier
MWS.Shaft.PGB.Carrier.GR_gb = info_PGB.NCqNemC; %GR = NS/Nint
if GBConfig == 5 || GBConfig == 6 % sun gear connected to HPS
    MWS.Shaft.PGB.Carrier.GR_int = info_PGB.NemqNhps; %GR = Nint/Ncomp
elseif GBConfig == 2 || GBConfig == 4 % sun gear connected to LPS
    MWS.Shaft.PGB.Carrier.GR_int = info_PGB.NemqNlps; %GR = Nint/Ncomp
else
    MWS.Shaft.PGB.Carrier.GR_int = 1;
end
MWS.Shaft.PGB.Carrier.J_go = info_PGB.JCo; %inertia of gear only, slug-ft2
% -- Planets
MWS.Shaft.PGB.Planet.GR_gb = 1; %GR = NP/Nem
MWS.Shaft.PGB.Planet.nP = info_PGB.nP; %number of planets
MWS.Shaft.PGB.Planet.mP = info_PGB.mP; %mass of planets, slugs
MWS.Shaft.PGB.Planet.J_go = info_PGB.JP; %inertia of gear only, slug-ft2
% Power System
% -- Sun Gear EM
MWS.PowerSystem.EMS.NMax = max(abs(info_PGB.NemSrange));
MWS.PowerSystem.EMS.PMaxC = info_PGB.PemS;
MWS.PowerSystem.EMS.PMax = info_PGB.PemS*Psf;
MWS.PowerSystem.EMS.TrqMaxC = info_PGB.TrqemS;
MWS.PowerSystem.EMS.TrqMax = info_PGB.TrqemS*Trqsf;
MWS.PowerSystem.EMS.Mass = info_PGB.MemS;
MWS.PowerSystem.EMS.J = info_PGB.JemS;
% -- Ring Gear EM
MWS.PowerSystem.EMR.NMax = max(abs(info_PGB.NemRrange));
MWS.PowerSystem.EMR.PMaxC = info_PGB.PemR;
MWS.PowerSystem.EMR.PMax = info_PGB.PemR*Psf;
MWS.PowerSystem.EMR.TrqMaxC = info_PGB.TrqemR;
MWS.PowerSystem.EMR.TrqMax = info_PGB.TrqemR*Trqsf;
MWS.PowerSystem.EMR.Mass = info_PGB.MemR;
MWS.PowerSystem.EMR.J = info_PGB.JemR;
% -- Carrier EM
MWS.PowerSystem.EMC.NMax = max(abs(info_PGB.NemCrange));
MWS.PowerSystem.EMC.PMaxC = info_PGB.PemC;
MWS.PowerSystem.EMC.PMax = info_PGB.PemC*Psf;
MWS.PowerSystem.EMC.TrqMaxC = info_PGB.TrqemC;
MWS.PowerSystem.EMC.TrqMax = info_PGB.TrqemC*Trqsf;
MWS.PowerSystem.EMC.Mass = info_PGB.MemC;
MWS.PowerSystem.EMC.J = info_PGB.JemC;
% -- Planet EM
MWS.PowerSystem.EMP.NMax = 30000;
MWS.PowerSystem.EMP.PMaxC = 0;
MWS.PowerSystem.EMP.PMax = 0;
MWS.PowerSystem.EMP.TrqMaxC = 0;
MWS.PowerSystem.EMP.TrqMax = 0;
MWS.PowerSystem.EMP.Mass = 0;
MWS.PowerSystem.EMP.J = 0;
% -- Torque Coefficients
MWS.Shaft.PGB.TrqCoeff.CTrqSR = info_PGB.TrqCoeff.CTrqSR;
MWS.Shaft.PGB.TrqCoeff.CTrqSC = info_PGB.TrqCoeff.CTrqSC;
MWS.Shaft.PGB.TrqCoeff.CTrqSP = info_PGB.TrqCoeff.CTrqSP;
MWS.Shaft.PGB.TrqCoeff.CTrqRS = info_PGB.TrqCoeff.CTrqRS;
MWS.Shaft.PGB.TrqCoeff.CTrqRC = info_PGB.TrqCoeff.CTrqRC;
MWS.Shaft.PGB.TrqCoeff.CTrqRP = info_PGB.TrqCoeff.CTrqRP;
MWS.Shaft.PGB.TrqCoeff.CTrqCS = info_PGB.TrqCoeff.CTrqCS;
MWS.Shaft.PGB.TrqCoeff.CTrqCR = info_PGB.TrqCoeff.CTrqCR;
MWS.Shaft.PGB.TrqCoeff.CTrqCP = info_PGB.TrqCoeff.CTrqCP;
MWS.Shaft.PGB.Sun.J = info_PGB.JS;
MWS.Shaft.PGB.Ring.J = info_PGB.JR;
MWS.Shaft.PGB.Carrier.J = info_PGB.JC;
MWS.Shaft.PGB.Planet.J = info_PGB.JP;
MWS.Shaft.PGB.Sun.JEff = info_PGB.JSEff;
MWS.Shaft.PGB.Ring.JEff = info_PGB.JREff;
MWS.Shaft.PGB.Carrier.JEff = info_PGB.JCEff;
% TEEM Power Constraints
MWS.PowerSystem.PwrMaxCoup_TEEMa = info_PGB.PemCoup_TEEMonlyA;
MWS.PowerSystem.PwrMaxCoup_TEEMd = info_PGB.PemCoup_TEEMonlyD;
MWS.PowerSystem.PwrMaxH_TEEM = info_PGB.PemH_TEEMonly;

% Effective Spool Inertias
if MWS.In.Options.EngineEMInt == 1
    MWS.Shaft.JlpsEff = MWS.Shaft.DEM.LPS_Inertia;
    MWS.Shaft.JhpsEff = MWS.Shaft.DEM.HPS_Inertia;
else
    MWS.Shaft.JlpsEff = info_PGB.JlpsEff;
    MWS.Shaft.JhpsEff = info_PGB.JhpsEff;
end

if HybConfig == 2
    MWS.PowerSystem.ESD.EnergyCapacity = 150;
else
    if MWS.In.Options.EngineEMInt == 1
        MWS.PowerSystem.ESD.EnergyCapacity = info_DEM.EnergyStorage;
    else
        MWS.PowerSystem.ESD.EnergyCapacity = info_PGB.EnergyStorage;
    end
end

% Final Mass Estimates ------------------------------------------------%

MWS.PowerSystem.PwrSysMass_DEM = info_DEM.Mps;
MWS.PowerSystem.PwrSysMass_PGB = info_PGB.Mps;
MWS.PowerSystem.TotalMass_DEM = info_DEM.Mtot;
MWS.PowerSystem.TotalMass_PGB = info_PGB.Mtot;

% END --------------------------------------------------------------------%