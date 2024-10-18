function [Mass,info] = DEMSizerFun(HybConfig,Nem_max,Nspool_max, ...
    Nlps,Nhps,Nhps_EPT,Nlps_EPT,Plps_nom,Phps_nom,Plps_nom_EPT,Phps_nom_EPT, ...
    PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PmaxPerEM, ...
    energyDenESD,PwrElecPwrDen,transientDur,Psf,Trqsf,Esf,GBsf,DiIS,BearingGap, ...
    NtsqNhps,NtsqNlps,UseToggle)

[N,M] = size(Nlps);
if N > M
    Nlps = Nlps';
end
[N,M] = size(Nhps);
if N > M
    Nhps = Nhps';
end
[N,M] = size(Nlps_EPT);
if N > M
    Nlps_EPT = Nlps_EPT';
end
[N,M] = size(Nhps_EPT);
if N > M
    Nhps_EPT = Nhps_EPT';
end
[N,M] = size(Plps_nom);
if N > M
    Plps_nom = Plps_nom';
end
[N,M] = size(Phps_nom);
if N > M
    Phps_nom = Phps_nom';
end
[N,M] = size(Plps_nom_EPT);
if N > M
    Plps_nom_EPT = Plps_nom_EPT';
end
[N,M] = size(Phps_nom_EPT);
if N > M
    Phps_nom_EPT = Phps_nom_EPT';
end

tol = 250;
NemH_max = 0;
PemH_max = 100;
NemL_max = 0;
PemL_max = 100;
iter_max = 10;
iter = 0;
k = 1;

while k == 1 || abs(NemH_max-min([5e5/sqrt(PemH_max*1.34102) Nem_max(1)])) > tol ...
    || abs(NemL_max-min([5e5/sqrt(PemL_max*1.34102) Nem_max(2)])) > tol

    if abs(NemH_max-min([5e5/sqrt(PemH_max*1.34102) Nem_max(1)])) > tol
        NemH_max = (1-k)*NemH_max + k*min([5e5/sqrt(PemH_max*1.34102) Nem_max(1)]);
    end

    if abs(NemL_max-min([5e5/sqrt(PemL_max*1.34102) Nem_max(2)])) > tol
        NemL_max = (1-k)*NemL_max + k*min([5e5/sqrt(PemL_max*1.34102) Nem_max(2)]);
    end

    % Maximum shaft speeds
    Nlps_max = Nspool_max(1); % maximum LPS speed, rpm
    Nhps_max = Nspool_max(2); % maximum HPS speed, rpm
    
    % shaft inertias
    Jlps = ((450102/(3.1*3.1)) + 17796 + 16237) /32.2/144.; %LPS inertia
    Jhps = 8627/32.2/144.; %HPS inertia
    
    % EM inertia approximation
    Vtip = 100/0.3048; %maximum EM rotor tip velocity
        
    % determine speeds of components
    NemHqNhps = NemH_max/Nhps_max;
    NemLqNlps = NemL_max/Nlps_max;
    NemH = Nhps*NemHqNhps;
    NemL = Nlps*NemLqNlps;
    NemH_EPT = Nhps_EPT*NemHqNhps;
    NemL_EPT = Nlps_EPT*NemLqNlps;
    
    % Powers, Torques, Energy (TEEM)
    [PemH_nom_MA, PemH_EPT_MA, PemH_TEEMa_MA, PemH_TEEMd_MA, ...
     TrqemH_nom_MA, TrqemH_EPT_MA, TrqemH_TEEMa_MA, TrqemH_TEEMd_MA, ...
     PemL_nom_MA, PemL_EPT_MA, PemL_TEEMa_MA, PemL_TEEMd_MA, ...
     TrqemL_nom_MA, TrqemL_EPT_MA, TrqemL_TEEMa_MA, TrqemL_TEEMd_MA, ...
     EnergyUse_TEEMa] ...
     = DEM_PwrTrq(HybConfig,NemH,NemL,NemH_EPT,NemL_EPT,...
                  Plps_nom,Phps_nom,Plps_nom_EPT,Phps_nom_EPT,PEx_ACsys,PX_EPT, ...
                  Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,transientDur,Psf,Trqsf);
    
    % Electric Power System Mass
    PemH_max = max(UseToggle.*[PemH_nom_MA PemH_EPT_MA PemH_TEEMa_MA PemH_TEEMd_MA]);%Psf
    TrqemH_max = max(UseToggle.*[TrqemH_nom_MA TrqemH_EPT_MA TrqemH_TEEMa_MA TrqemH_TEEMd_MA]);%*Trqsf;
    PemL_max = max(UseToggle.*[PemL_nom_MA PemL_EPT_MA PemL_TEEMa_MA PemL_TEEMd_MA]);%Psf
    TrqemL_max = max(UseToggle.*[TrqemL_nom_MA TrqemL_EPT_MA TrqemL_TEEMa_MA TrqemL_TEEMd_MA]);%*Trqsf;
    numEMH = ceil(PemH_max/PmaxPerEM);
    MemH = numEMH*1.2866*0.0685218*0.3462*(TrqemH_max/max([numEMH, 1]))^0.7486;
    numEML = ceil(PemL_max/PmaxPerEM);
    MemL = numEML*1.2866*0.0685218*0.3462*(TrqemL_max/max([numEML, 1]))^0.7486;
    if HybConfig == 2
        Mesd = 0; %boost leverages large existing battery (not included)
    else
        Mesd = max([EnergyUse_TEEMa 0])/energyDenESD * Esf;
    end
    MpeH = PemH_max/PwrElecPwrDen;
    MpeL = PemL_max/PwrElecPwrDen;
    Pesd = (EnergyUse_TEEMa/(transientDur*745.7*2.77778e-7));
    if HybConfig == 2
        MpeConv = max(Plps_nom)/PwrElecPwrDen; %boost leverages converter sized for boost (assumed to handle TEEM for short durations)
    else
        MpeConv = Pesd/PwrElecPwrDen;
    end
    MassE = MemH + MemL + Mesd + MpeH + MpeL + MpeConv; %mass of electrical system, slug

    % Properties
    sigy = 290*145.038; %yield stress, MPa->psi
    rho = 7900*1.12287e-6; %density, kg/m3->slug/in3
    L = 60; %length of shaft, in
    minSWT = 0.125; %minimum shaft wall thickness, in

    % Spool powers and torques
    [Phps_nom_MA, Phps_EPT_MA, Phps_TEEMa_MA, Phps_TEEMd_MA, ...
          Trqhps_nom_MA, Trqhps_EPT_MA, Trqhps_TEEMa_MA, Trqhps_TEEMd_MA, ...
          Plps_nom_MA, Plps_EPT_MA, Plps_TEEMa_MA, Plps_TEEMd_MA, ...
          Trqlps_nom_MA, Trqlps_EPT_MA, Trqlps_TEEMa_MA, Trqlps_TEEMd_MA] ...
          = PwrTrqGB(HybConfig,Nhps,Nlps,Nhps_EPT,Nlps_EPT, ...
            Plps_nom,Phps_nom,Plps_nom_EPT,Phps_nom_EPT,PEx_ACsys,PX_EPT, ...
            Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,Psf,Trqsf);
    Trqhps_max = max(UseToggle.*[Trqhps_nom_MA, Trqhps_EPT_MA, Trqhps_TEEMa_MA, Trqhps_TEEMd_MA]);
    Trqlps_max = max(UseToggle.*[Trqlps_nom_MA, Trqlps_EPT_MA, Trqlps_TEEMa_MA, Trqlps_TEEMd_MA]);

    % Gear/Speed ratios
    NemLqNts = abs(NemLqNlps/NtsqNlps);
    NemHqNts = abs(NemHqNhps/NtsqNhps);
    
    % Tower Shafts
    % HP Shaft (inner shaft)
    Do = DiIS*linspace(1.1,5,100);
    p1 = polyfit(log(Do./DiIS),log((Do.^2-DiIS^2)),1); 
    mtsh1 = p1(1);
    btsh1 = exp(p1(2));
    Trqmax_sh = Trqhps_max/NtsqNhps;
    DoqDitshI = ((Trqmax_sh*12)*(16/pi)*(GBsf/sigy)*(1/btsh1))^(1/mtsh1);
    if (DoqDitshI-1)*DiIS < 2*minSWT %enforce min shaft wall thickness
        DoqDitshI = (DiIS+2*minSWT)/DiIS;
    end
    Mtsh = rho*L*(pi/4)*((DoqDitshI*DiIS)^2-DiIS^2); %mass of HP tower shaft, slug
    % LP Shaft (outer shaft)
    DiOS = DoqDitshI*DiIS + 2*BearingGap;
    Do = DiOS*linspace(1.1,5,100);
    p2 = polyfit(log(Do./DiOS),log((Do.^2-DiOS^2)),1); 
    mtsh2 = p2(1);
    btsh2 = exp(p2(2));
    Trqmax_sl = Trqlps_max/NtsqNlps;
    DoqDitshO = ((Trqmax_sl*12)*(16/pi)*(GBsf/sigy)*(1/btsh2))^(1/mtsh2);
    if (DoqDitshO-1)*DiOS < 2*minSWT %enforce min shaft wall thickness
        DoqDitshO = (DiOS+2*minSWT)/DiOS;
    end
    Mtsl = rho*L*(pi/4)*((DoqDitshO*DiOS)^2-DiOS^2); %mass of LP tower shaft, slug
    % Total
    Mts = Mtsh + Mtsl;

    % Bevel Gears
    % --HP shaft - spool side
    nbg_h = min(NtsqNhps,1/NtsqNhps);
    Trqmax = max(abs([Trqhps_max/NtsqNhps, Trqhps_max]));
    Q_h = (Trqmax/5252.113)*(nbg_h+1)^3/nbg_h^2;
    Mbgh = (55/32.174)*Q_h^0.91;
    % --HP shaft - EM side
    nbg_emh = min(NemHqNts,1/NemHqNts);
    Trqmax = max(abs([Trqhps_max/NtsqNhps, Trqhps_max/NemHqNhps]));
    Q_emh = (Trqmax/5252.113)*(nbg_emh+1)^3/nbg_emh^2;
    Mbgemh = (55/32.174)*Q_emh^0.91;
    % --LP shaft - spool side
    nbg_l = min(NtsqNlps,1/NtsqNlps);
    Trqmax = max(abs([Trqlps_max/NtsqNlps, Trqlps_max]));
    Q_l = (Trqmax/5252.113)*(nbg_l+1)^3/nbg_l^2;
    Mbgl = (55/32.174)*Q_l^0.91;
    % --LP shaft - EM side
    nbg_eml = min(NemLqNts,1/NemLqNts);
    Trqmax = max(abs([Trqlps_max/NtsqNlps, Trqlps_max/NemLqNlps]));
    Q_eml = (Trqmax/5252.113)*(nbg_eml+1)^3/nbg_eml^2;
    Mbgeml = (55/32.174)*Q_eml^0.91;
    % --Total
    Mbg = Mbgemh + Mbgh + Mbgeml + Mbgl; %slug

    % Total Mech System Mass
    Mmech = Mts + Mbg; %slug

    % Overall System Mass
    Mass = MassE + Mmech; %slug

    %iteration check
    iter = iter + 1;
    k = 0.95*k;
    if iter >= iter_max
        disp('Unable to converge in Nem while-loop')
        break;
    end

end

% Inertia
Drotor_EMH = 60*Vtip/(pi*NemH_max);
JemH = (1/2)*MemH*(Drotor_EMH/2)^2;
Drotor_EML = 60*Vtip/(pi*NemL_max);
JemL = (1/2)*MemL*(Drotor_EML/2)^2;
JlpsEff = Jlps + numEML*JemL*NemLqNlps^2;
JhpsEff = Jhps + numEMH*JemH*NemHqNhps^2;

% pick off mass drivers for the motors and power electronics
%   1 - nominial use 
%   2 - EPT
%   3 - TEEM-accel
%   4 - TEEM-decel
[dumbie,MassDriver_emH] = max(UseToggle.*[TrqemH_nom_MA TrqemH_EPT_MA TrqemH_TEEMa_MA TrqemH_TEEMd_MA]);
[dumbie,MassDriver_emL] = max(UseToggle.*[TrqemL_nom_MA TrqemL_EPT_MA TrqemL_TEEMa_MA TrqemL_TEEMd_MA]);
[dumbie,MassDriver_peH] = max(UseToggle.*[PemH_nom_MA PemH_EPT_MA PemH_TEEMa_MA PemH_TEEMd_MA]);
[dumbie,MassDriver_peL] = max(UseToggle.*[PemL_nom_MA PemL_EPT_MA PemL_TEEMa_MA PemL_TEEMd_MA]);

info.NemHmax = NemH_max;
info.NemLmax = NemL_max;
info.NemHqNhps = NemHqNhps;
info.NemLqNlps = NemLqNlps;
info.NemHrange = [min(NemH) max(NemH)];
info.NemLrange = [min(NemL) max(NemL)];
info.PemH = PemH_max;
info.PemL = PemL_max;
info.Pesd = Pesd;
info.TrqemH = TrqemH_max;
info.TrqemL = TrqemL_max;
info.MemH = MemH;
info.MemL = MemL;
info.Mesd = Mesd;
info.MpeH = MpeH;
info.MpeL = MpeL;
info.MpeConv = MpeConv;
info.Mps = MassE;
info.JemH = JemH;
info.JemL = JemL;
info.JlpsEff = JlpsEff;
info.JhpsEff = JhpsEff;
info.MassDriver_emH = MassDriver_emH;
info.MassDriver_emL = MassDriver_emL;
info.MassDriver_peH = MassDriver_peH;
info.MassDriver_peL = MassDriver_peL;
info.EnergyStorage = EnergyUse_TEEMa*Esf;
info.Mts = Mts;
info.Mbg = Mbg;
info.Mmech = Mts + Mbg;
info.Mtot = Mass;
if iter >= iter_max
    info.NemMaxLoopPass = 0;
else
    info.NemMaxLoopPass = 1;
end

end