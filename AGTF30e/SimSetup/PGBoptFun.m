function [Mass,info] = PGBoptFun(X,GBConfig,HybConfig,rR,NHPGBqNem,Nem_max,Nspool_max, ...
    NPGB_max,Nlps,Nhps,Nhps_EPT,Nlps_EPT,Plps_nom,Phps_nom,Plps_nom_EPT,Phps_nom_EPT, ...
    PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,PmaxPerEM, ...
    energyDenESD,PwrElecPwrDen,transientDur,Psf,Trqsf,Esf,GBsf,DiIS,BearingGap,UseToggle)

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

% gearbox design parameters
rSqrR = X(1); % radius of sun / radius of ring
NLPGBqNem = X(2); % speed of PGB component connected to LPS / speed of EM connected to LPS
NtsqNhps = X(3); % speed of HPS tower shaft / speed of HPS
NtsqNlps = X(4); % speed of LPS tower shaft / speed of LPS
% NHPGBqNem = X(3); % speed of PGB component connected to HPS / speed of EM connected to HPS
% PemCoup_maxTEEM = X(4); % maximum power of the coupling EM with TEEM, hp
% if HybConfig > 1
%     PemCoup_maxNom = X(5); % maximum power of the coupling EM for nominal operation, hp
% else
%     PemCoup_maxNom = 0; % this is not an independent in the optimization for the standard engine
% end

% dumbie initialization
tol = 250;
NemS_max = Nem_max(1);
NemR_max = Nem_max(2);
NemC_max = Nem_max(3);
iter_max = 40;
iter = 0;
PXeffH = -1; 
if HybConfig == 1
    PemS_max = 100;
    PemR_max = 100;
    PemC_max = 100;
elseif HybConfig == 2
    if GBConfig == 1
        PemS_max = 400;
        PemR_max = 1200;
        PemC_max = PemCoup_maxNom;
    elseif GBConfig == 2
        PemS_max = 400;
        PemR_max = PemCoup_maxNom;
        PemC_max = 1200;
    elseif GBConfig == 3
        PemS_max = 1200;
        PemR_max = 400;
        PemC_max = PemCoup_maxNom;
    elseif GBConfig == 4
        PemS_max = PemCoup_maxNom;
        PemR_max = 400;
        PemC_max = 1200;
    elseif GBConfig == 5
        PemS_max = 1200;
        PemR_max = PemCoup_maxNom;
        PemC_max = 400;
    else
        PemS_max = PemCoup_maxNom;
        PemR_max = 1200;
        PemC_max = 400;
    end
else
    if GBConfig == 1
        PemS_max = 1200;
        PemR_max = 400;
        PemC_max = PemCoup_maxNom;
    elseif GBConfig == 2
        PemS_max = 1200;
        PemR_max = PemCoup_maxNom;
        PemC_max = 400;
    elseif GBConfig == 3
        PemS_max = 400;
        PemR_max = 1200;
        PemC_max = PemCoup_maxNom;
    elseif GBConfig == 4
        PemS_max = PemCoup_maxNom;
        PemR_max = 1200;
        PemC_max = 400;
    elseif GBConfig == 5
        PemS_max = 400;
        PemR_max = PemCoup_maxNom;
        PemC_max = 1200;
    else
        PemS_max = PemCoup_maxNom;
        PemR_max = 400;
        PemC_max = 1200;
    end
end
k = 1;

while k == 1 || abs(NemS_max-min([5e5/sqrt(PemS_max*1.34102) Nem_max(1)])) > tol ...
    || abs(NemR_max-min([5e5/sqrt(PemR_max*1.34102) Nem_max(2)])) > tol ...
    || abs(NemC_max-min([5e5/sqrt(PemC_max*1.34102) Nem_max(3)])) > tol

    if max(PXeffH) >= 0 % no point optimizing a design that won't work
        break;
    end

    if abs(NemS_max-min([5e5/sqrt(PemS_max*1.34102) Nem_max(1)])) > tol
        NemS_max = (1-k)*NemS_max + k*min([5e5/sqrt(PemS_max*1.34102) Nem_max(1)]);
    end

    if abs(NemR_max-min([5e5/sqrt(PemR_max*1.34102) Nem_max(2)])) > tol
        NemR_max = (1-k)*NemR_max + k*min([5e5/sqrt(PemR_max*1.34102) Nem_max(2)]);
    end

    if abs(NemC_max-min([5e5/sqrt(PemC_max*1.34102) Nem_max(3)])) > tol
        NemC_max = (1-k)*NemC_max + k*min([5e5/sqrt(PemC_max*1.34102) Nem_max(3)]);
    end

    % % EM maximum speed parameters
    % NemS_max = Nem_max(1); % maximum speed of the sun gear EM, rpm
    % NemR_max = Nem_max(2); % maximum speed of the ring gear EM, rpm
    % NemC_max = Nem_max(3); % maximum speed of the carrier, rpm
    
    % Maximum shaft speeds
    Nlps_max = Nspool_max(1); % maximum LPS speed, rpm
    Nhps_max = Nspool_max(2); % maximum HPS speed, rpm
    
    % Maximum gearbox speeds
    NS_max = NPGB_max(1);
    NR_max = NPGB_max(2);
    NC_max = NPGB_max(3);
    NP_max = NPGB_max(4);
    
    % shaft inertias
    Jlps = ((450102/(3.1*3.1)) + 17796 + 16237) /32.2/144.; %LPS inertia
    Jhps = 8627/32.2/144.; %HPS inertia
    
    % EM inertia approximation
    Vtip = 100/0.3048; %maximum EM rotor tip velocity
    
    % important radii
    rS = rSqrR*rR; % radius of sun
    rP = (1/2)*(rR-rS);
    
    % number of planets
    rSqrR_vec = [0 0.3 0.33 0.5 0.625 0.9 1];
    nP_vec = [3 3 4 4 5 8 8];
    nP = round(interp1(rSqrR_vec,nP_vec,rSqrR));
    
    % parameters for approximating compononent weights and inertias
    rho = 490/32.174; % density of steel, slug/ft^3
    tR = 0.5/12; % thickness of the ring gear, ft
    tC = 0.5/12; %thickness of the carrier plate, ft
    thicknessPGB = 1.5/12; %thickness of ring gear, planets, and sun gear
    
    % masses and isolated component inertias
    mP = rho*(pi*rP^2*thicknessPGB); % mass of planet, slug
    JP = (1/2)*mP*rP^2;
    mS = rho*(pi*rS^2*thicknessPGB); % mass of sun gear, slug
    JSo = (1/2)*mS*rS^2;
    mR1 = rho*(pi*((rR+tR)^2-rR^2)*thicknessPGB); % mass of ring gear
    mR2 = rho*(pi*(rR^2)*tR); % mass of ring gear housing
    JRo = (1/2)*mR1*(rR^2+(rR+tR)^2) + (1/2)*mR2*rR^2;
    if nP == 3 % carrier is a triangle
        height = (3/2)*(rS+rP);
        base = 2*height*tand(30);
        mC = rho*(1/2)*base*height*tC; 
        JCo = (1/18)*mC*((3/2)*(rS+rP))^2;
    elseif nP == 4 % carreir is a square
        sides = (rS+rP)*tand(45)*2;
        mC = rho*sides^2*tC;
        JCo = (1/6)*mC*sides^2;
    elseif nP == 5 % carrier is a pentagon
        sides = (rS+rP)*sind(36)*2;
        mC = rho*(1/4)*sqrt(5*(5+2*sqrt(5)))*sides^2*tC;
        JCo = (2+3/sqrt(5))*(1/6)*mC*sides^2;
    elseif nP == 6 % carrier is a hexagon
        sides = (rS+rP)*sind(30)*2;
        mC = rho*3*sqrt(3)/2*sides^2*tC;
        JCo = (5/6)*mC*sides^2;
    else % carrier is a circle
        mC = rho*(pi*(rS+rP)^2*tC);
        JCo = (1/2)*mC*(rS+rP)^2;
    end
    
    % TEEM simplification 
    %   splits the transient into 5 section and weights the the power transfer
    %   effect for each section.
    dNlps = (Nlps(end)-Nlps(1))/5;
    Nlps_bins = Nlps(1) + [0 dNlps 2*dNlps 3*dNlps 4*dNlps 5*dNlps];
    idx1 = 1;
    for i = 1:length(Nlps)
        if Nlps(i) <= Nlps_bins(2)
            idx2 = i;
        elseif Nlps(i) <= Nlps_bins(3)
            idx3 = i;
        elseif Nlps(i) <= Nlps_bins(4)
            idx4 = i;
        elseif Nlps(i) <= Nlps_bins(5)
            idx5 = i;
        else
            break
        end
    end
    idx6 = length(Nlps);
    TEEMaWeights = [0.4 0.3 0.2 0.1 0]; % acceleraiton weightings
    TEEMdWeights = [0.2 0.2 0.2 0.2 0.2]; % deceleraiton weightings

    if GBConfig == 1 % Sun - HP, Ring - LP, Carrier - Coupling
        
        % determine speeds of components
        % -- can infer the speed ratio of the components connected to the LPS
        %    and HPS relative to the LPS and HPS
        NemSqNhps = NemS_max/Nhps_max;
        NemRqNlps = NemR_max/Nlps_max;
        NemqNhps = NemSqNhps;
        NemqNlps = NemRqNlps;
        NSqNemS = NHPGBqNem;
        NRqNemR = NLPGBqNem;
        % -- planetary gearbox speeds
        NemS = Nhps*NemSqNhps;
        NS = NemS*NSqNemS;
        NemR = Nlps*NemRqNlps;
        NR = NemR*NRqNemR;
        NC = (rSqrR*NS+NR)./(rSqrR+1);
        NCqNemC = max(abs(NC))/NemC_max;
        NemC = NC./NCqNemC;
        NemS_EPT = Nhps_EPT*NemSqNhps;
        NS_EPT = NemS_EPT*NSqNemS;
        NemR_EPT = Nlps_EPT*NemRqNlps;
        NR_EPT = NemR_EPT*NRqNemR;
        NC_EPT = (rSqrR*NS_EPT+NR_EPT)./(rSqrR+1);
        NemC_EPT = NC_EPT./NCqNemC;
        % Cast EM speeds for PGB_PwrTrq.m
        NemH = NemS;
        NemL = NemR;
        NemCoup = NemC;
        NemH_EPT = NemS_EPT;
        NemCoup_EPT = NemC_EPT;
        % Check PGB speeds vs. max
        % -- planetary gearbox inertias
        JS = JSo + Jhps/abs(NemSqNhps*NSqNemS)^2; % does not include EM(s) yet
        JC = JCo; % does not include EM(s) yet
        JR = JRo + Jlps/abs(NemRqNlps*NRqNemR)^2; % does not include EM(s) yet
        JErr = inf;
        iterCnt = 0;
        while JErr > 0.02 && iterCnt <= 7
            % Gearbox effect
            [TrqCoeff,JEff] = PGBimpact(rS,rR,nP,mP,JS,JR,JC,JP);
            PXeffH = TrqCoeff.CTrqSC*NS./NC;
            PXeffL = TrqCoeff.CTrqRC*NR./NC;
            PXeffH_EPT = TrqCoeff.CTrqSC*NS_EPT./NC_EPT;
            PXeffL_EPT = TrqCoeff.CTrqRC*NR_EPT./NC_EPT;
            PXeffH_bins = [trapz(Nlps(idx1:idx2),PXeffH(idx1:idx2))/dNlps, trapz(Nlps(idx2:idx3),PXeffH(idx2:idx3))/dNlps, ...
                trapz(Nlps(idx3:idx4),PXeffH(idx3:idx4))/dNlps, trapz(Nlps(idx4:idx5),PXeffH(idx4:idx5))/dNlps, ...
                trapz(Nlps(idx5:idx6),PXeffH(idx5:idx6))/dNlps];
            % -- EM power
            if max(PXeffH) >= 0 
                MassE = 1000 + (mS + mR1 + mR2 + mC + mP); %slug
                disp('WARNING: The gearbox design does not favor the LPS across the entire engine power range and will not be consistent with the control algorithm')
            else
                % Powers, Torques, Energy (TEEM)
                [PemS_nom_MA, PemS_EPT_MA, PemS_TEEMa_MA, PemS_TEEMd_MA, ...
                  TrqemS_nom_MA, TrqemS_EPT_MA, TrqemS_TEEMa_MA, TrqemS_TEEMd_MA, ...
                  PemR_nom_MA, PemR_EPT_MA, PemR_TEEMa_MA, PemR_TEEMd_MA, ...
                  TrqemR_nom_MA, TrqemR_EPT_MA, TrqemR_TEEMa_MA, TrqemR_TEEMd_MA, ...
                  PemC_nom_MA, PemC_EPT_MA, PemC_TEEMa_MA, PemC_TEEMd_MA, ...
                  TrqemC_nom_MA, TrqemC_EPT_MA, TrqemC_TEEMa_MA, TrqemC_TEEMd_MA, ...
                  PemCoup_TEEMonlyA_MA, PemCoup_TEEMonlyD_MA, PemH_TEEMonly, EnergyUse_TEEMa] ...
                  = PGB_PwrTrq(HybConfig,GBConfig,NemH,NemL,NemCoup,NemH_EPT,NemCoup_EPT, ...
                    PXeffH,PXeffL,PXeffH_EPT,PXeffL_EPT,Plps_nom,Phps_nom,Plps_nom_EPT, ...
                    Phps_nom_EPT,PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,...
                    TEEMaWeights,TEEMdWeights,PXeffH_bins,transientDur,Psf,Trqsf);
                % Mass
                PemS_max = max(UseToggle.*[PemS_nom_MA PemS_EPT_MA PemS_TEEMa_MA PemS_TEEMd_MA]);%Psf;
                TrqemS_max = max(UseToggle.*[TrqemS_nom_MA TrqemS_EPT_MA TrqemS_TEEMa_MA TrqemS_TEEMd_MA]);%*Trqsf;
                PemR_max = max(UseToggle.*[PemR_nom_MA PemR_EPT_MA PemR_TEEMa_MA PemR_TEEMd_MA]);%*Psf;
                TrqemR_max = max(UseToggle.*[TrqemR_nom_MA TrqemR_EPT_MA TrqemR_TEEMa_MA TrqemR_TEEMd_MA]);%*Trqsf;
                PemC_max = max(UseToggle.*[PemC_nom_MA PemC_EPT_MA PemC_TEEMa_MA PemC_TEEMd_MA]);%*Psf;
                TrqemC_max = max(UseToggle.*[TrqemC_nom_MA TrqemC_EPT_MA TrqemC_TEEMa_MA TrqemC_TEEMd_MA]);%*Trqsf;
                numEMS = ceil(PemS_max/PmaxPerEM);
                MemS = numEMS*1.2866*0.0685218*0.3462*(TrqemS_max/max([numEMS, 1]))^0.7486;
                numEMR = ceil(PemR_max/PmaxPerEM);
                MemR = numEMR*1.2866*0.0685218*0.3462*(TrqemR_max/max([numEMR, 1]))^0.7486;
                numEMC = ceil(PemC_max/PmaxPerEM);
                MemC = numEMC*1.2866*0.0685218*0.3462*(TrqemC_max/max([numEMC, 1]))^0.7486;
                if HybConfig == 2
                    Mesd = 0; %boost leverages large existing battery (not included)
                else
                    Mesd = max([EnergyUse_TEEMa 0])/energyDenESD * Esf;
                end
                MpeS = PemS_max/PwrElecPwrDen;
                MpeR = PemR_max/PwrElecPwrDen;
                MpeC = PemC_max/PwrElecPwrDen;
                if UseToggle(3) ~= 1
                    Pesd = 0;
                else
                    Pesd = (EnergyUse_TEEMa/(transientDur*745.7*2.77778e-7));
                end
                if HybConfig == 2
                    MpeConv = max(Plps_nom)/PwrElecPwrDen; %boost leverages converter sized for boost (assumed to handle TEEM for short durations)
                else
                    MpeConv = Pesd/PwrElecPwrDen;
                end
                MassE = MemS + MemR + MemC + Mesd + MpeS + MpeR + MpeC + MpeConv; %mass of electrical system, slug
                % Inertia update
                JSold = JS;
                JRold = JR;
                JCold = JC;
                Drotor_EMS = 60*Vtip/(pi*NemS_max);
                JemS = (1/2)*MemS*(Drotor_EMS/2)^2;
                Drotor_EMR = 60*Vtip/(pi*NemR_max);
                JemR = (1/2)*MemR*(Drotor_EMR/2)^2;
                Drotor_EMC = 60*Vtip/(pi*NemC_max);
                JemC = (1/2)*MemC*(Drotor_EMC/2)^2;
                JS = JSo + Jhps/abs(NemSqNhps*NSqNemS)^2 + numEMS*JemS/abs(NSqNemS)^2; 
                JC = JCo + numEMC*JemC/abs(NCqNemC)^2;
                JR = JRo + Jlps/abs(NemRqNlps*NRqNemR)^2 + numEMR*JemR/abs(NRqNemR)^2; 
                JErr = max(abs(([JS JR JC]-[JSold JRold JCold])./[JSold JRold JCold]));
            end
            iterCnt = iterCnt + 1;
        end 
        JlpsEff = JEff.Ring*abs(NemRqNlps*NRqNemR)^2;
        JhpsEff = JEff.Sun*abs(NemSqNhps*NSqNemS)^2;
    
    elseif GBConfig == 2 % Sun - HP, Ring - Coupling, Carrier - LP
    
        % determine speeds of components
        % -- can infer the speed ratio of the components connected to the LPS
        %    and HPS relative to the LPS and HPS
        NemSqNhps = NemS_max/Nhps_max;
        NemCqNlps = NemC_max/Nlps_max;
        NemqNhps = NemSqNhps;
        NemqNlps = NemCqNlps;
        NSqNemS = NHPGBqNem;
        NCqNemC = NLPGBqNem;
        % -- planetary gearbox speeds
        NemS = Nhps*NemSqNhps;
        NS = NemS*NSqNemS;
        NemC = Nlps*NemCqNlps;
        NC = NemC*NCqNemC;
        NR = (rSqrR+1)*NC-rSqrR*NS;
        NRqNemR = max(abs(NR))/NemR_max;
        NemR = NR./NRqNemR;
        NemS_EPT = Nhps_EPT*NemSqNhps;
        NS_EPT = NemS_EPT*NSqNemS;
        NemC_EPT = Nlps_EPT*NemCqNlps;
        NC_EPT = NemC_EPT*NCqNemC;
        NR_EPT = (rSqrR+1)*NC_EPT-rSqrR*NS_EPT;
        NemR_EPT = NR_EPT./NRqNemR;
        % Cast EM speeds for PGB_PwrTrq.m
        NemH = NemS;
        NemL = NemC;
        NemCoup = NemR;
        NemH_EPT = NemS_EPT;
        NemCoup_EPT = NemR_EPT;
        % -- planetary gearbox inertias
        JS = JSo + Jhps/abs(NemSqNhps*NSqNemS)^2; % does not include EM(s) yet
        JC = JCo + Jlps/abs(NemCqNlps*NCqNemC)^2; % does not include EM(s) yet
        JR = JRo; % does not include EM(s) yet
        JErr = inf;
        iterCnt = 0;
        while JErr > 0.02 && iterCnt <= 7
            % Gearbox effect
            [TrqCoeff,JEff] = PGBimpact(rS,rR,nP,mP,JS,JR,JC,JP);
            PXeffH = TrqCoeff.CTrqSR*NS./NR;
            PXeffL = TrqCoeff.CTrqCR*NC./NR;
            PXeffH_EPT = TrqCoeff.CTrqSR*NS_EPT./NR_EPT;
            PXeffL_EPT = TrqCoeff.CTrqCR*NC_EPT./NR_EPT;
            PXeffH_bins = [trapz(Nlps(idx1:idx2),PXeffH(idx1:idx2))/dNlps, trapz(Nlps(idx2:idx3),PXeffH(idx2:idx3))/dNlps, ...
                trapz(Nlps(idx3:idx4),PXeffH(idx3:idx4))/dNlps, trapz(Nlps(idx4:idx5),PXeffH(idx4:idx5))/dNlps, ...
                trapz(Nlps(idx5:idx6),PXeffH(idx5:idx6))/dNlps];
            % -- EM power
            if max(PXeffH) >= 0 
                MassE = 1000 + (mS + mR1 + mR2 + mC + mP); %slug
                disp('WARNING: The gearbox design does not favor the LPS across the entire engine power range and will not be consistent with the control algorithm')
            else
                % Powers, Torques, Energy (TEEM)
                [PemS_nom_MA, PemS_EPT_MA, PemS_TEEMa_MA, PemS_TEEMd_MA, ...
                  TrqemS_nom_MA, TrqemS_EPT_MA, TrqemS_TEEMa_MA, TrqemS_TEEMd_MA, ...
                  PemR_nom_MA, PemR_EPT_MA, PemR_TEEMa_MA, PemR_TEEMd_MA, ...
                  TrqemR_nom_MA, TrqemR_EPT_MA, TrqemR_TEEMa_MA, TrqemR_TEEMd_MA, ...
                  PemC_nom_MA, PemC_EPT_MA, PemC_TEEMa_MA, PemC_TEEMd_MA, ...
                  TrqemC_nom_MA, TrqemC_EPT_MA, TrqemC_TEEMa_MA, TrqemC_TEEMd_MA, ...
                  PemCoup_TEEMonlyA_MA, PemCoup_TEEMonlyD_MA, PemH_TEEMonly, EnergyUse_TEEMa] ...
                  = PGB_PwrTrq(HybConfig,GBConfig,NemH,NemL,NemCoup,NemH_EPT,NemCoup_EPT, ...
                    PXeffH,PXeffL,PXeffH_EPT,PXeffL_EPT,Plps_nom,Phps_nom,Plps_nom_EPT, ...
                    Phps_nom_EPT,PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,...
                    TEEMaWeights,TEEMdWeights,PXeffH_bins,transientDur,Psf,Trqsf);
                % Mass
                PemS_max = max(UseToggle.*[PemS_nom_MA PemS_EPT_MA PemS_TEEMa_MA PemS_TEEMd_MA]);%Psf;
                TrqemS_max = max(UseToggle.*[TrqemS_nom_MA TrqemS_EPT_MA TrqemS_TEEMa_MA TrqemS_TEEMd_MA]);%*Trqsf;
                PemR_max = max(UseToggle.*[PemR_nom_MA PemR_EPT_MA PemR_TEEMa_MA PemR_TEEMd_MA]);%*Psf;
                TrqemR_max = max(UseToggle.*[TrqemR_nom_MA TrqemR_EPT_MA TrqemR_TEEMa_MA TrqemR_TEEMd_MA]);%*Trqsf;
                PemC_max = max(UseToggle.*[PemC_nom_MA PemC_EPT_MA PemC_TEEMa_MA PemC_TEEMd_MA]);%*Psf;
                TrqemC_max = max(UseToggle.*[TrqemC_nom_MA TrqemC_EPT_MA TrqemC_TEEMa_MA TrqemC_TEEMd_MA]);%*Trqsf;
                numEMS = ceil(PemS_max/PmaxPerEM);
                MemS = numEMS*1.2866*0.0685218*0.3462*(TrqemS_max/max([numEMS, 1]))^0.7486;
                numEMR = ceil(PemR_max/PmaxPerEM);
                MemR = numEMR*1.2866*0.0685218*0.3462*(TrqemR_max/max([numEMR, 1]))^0.7486;
                numEMC = ceil(PemC_max/PmaxPerEM);
                MemC = numEMC*1.2866*0.0685218*0.3462*(TrqemC_max/max([numEMC, 1]))^0.7486;
                if HybConfig == 2
                    Mesd = 0; %boost leverages large existing battery (not included)
                else
                    Mesd = max([EnergyUse_TEEMa 0])/energyDenESD * Esf;
                end
                MpeS = PemS_max/PwrElecPwrDen;
                MpeR = PemR_max/PwrElecPwrDen;
                MpeC = PemC_max/PwrElecPwrDen;
                if UseToggle(3) ~= 1
                    Pesd = 0;
                else
                    Pesd = (EnergyUse_TEEMa/(transientDur*745.7*2.77778e-7));
                end
                if HybConfig == 2
                    MpeConv = max(Plps_nom)/PwrElecPwrDen; %boost leverages converter sized for boost (assumed to handle TEEM for short durations)
                else
                    MpeConv = Pesd/PwrElecPwrDen;
                end
                MassE = MemS + MemR + MemC + Mesd + MpeS + MpeR + MpeC + MpeConv; %mass of electrical system, slug
                % Inertia update
                JSold = JS;
                JRold = JR;
                JCold = JC;
                Drotor_EMS = 60*Vtip/(pi*NemS_max);
                JemS = (1/2)*MemS*(Drotor_EMS/2)^2;
                Drotor_EMR = 60*Vtip/(pi*NemR_max);
                JemR = (1/2)*MemR*(Drotor_EMR/2)^2;
                Drotor_EMC = 60*Vtip/(pi*NemC_max);
                JemC = (1/2)*MemC*(Drotor_EMC/2)^2;
                JS = JSo + Jhps/abs(NemSqNhps*NSqNemS)^2 + numEMS*JemS/abs(NSqNemS)^2; 
                JC = JCo + Jlps/abs(NemCqNlps*NCqNemC)^2 + numEMC*JemC/abs(NCqNemC)^2;
                JR = JRo + numEMR*JemR/abs(NRqNemR)^2; 
                JErr = max(abs(([JS JR JC]-[JSold JRold JCold])./[JSold JRold JCold]));
            end
            iterCnt = iterCnt + 1;
        end 
        JlpsEff = JEff.Carrier*abs(NemCqNlps*NCqNemC)^2;
        JhpsEff = JEff.Sun*abs(NemSqNhps*NSqNemS)^2;
    
    elseif GBConfig == 3 % Sun - LP, Ring - HP, Carrier - Coupling
    
        % determine speeds of components
        % -- can infer the speed ratio of the components connected to the LPS
        %    and HPS relative to the LPS and HPS
        NemRqNhps = NemR_max/Nhps_max;
        NemSqNlps = NemS_max/Nlps_max;
        NemqNhps = NemRqNhps;
        NemqNlps = NemSqNlps;
        NRqNemR = NHPGBqNem;
        NSqNemS = NLPGBqNem;
        % -- planetary gearbox speeds
        NemR = Nhps*NemRqNhps;
        NR = NemR*NRqNemR;
        NemS = Nlps*NemSqNlps;
        NS = NemS*NSqNemS;
        NC = (rSqrR*NS+NR)./(rSqrR+1);
        NCqNemC = max(abs(NC))/NemC_max;
        NemC = NC./NCqNemC;
        NemR_EPT = Nhps_EPT*NemRqNhps;
        NR_EPT = NemR_EPT*NRqNemR;
        NemS_EPT = Nlps_EPT*NemSqNlps;
        NS_EPT = NemS_EPT*NSqNemS;
        NC_EPT = (rSqrR*NS_EPT+NR_EPT)./(rSqrR+1);
        NemC_EPT = NC_EPT./NCqNemC;
        % Cast EM speeds for PGB_PwrTrq.m
        NemH = NemR;
        NemL = NemS;
        NemCoup = NemC;
        NemH_EPT = NemR_EPT;
        NemCoup_EPT = NemC_EPT;
        % -- planetary gearbox inertias
        JS = JSo + Jlps/abs(NemSqNlps*NSqNemS)^2; % does not include EM(s) yet
        JC = JCo; % does not include EM(s) yet
        JR = JRo + Jhps/abs(NemRqNhps*NRqNemR)^2; % does not include EM(s) yet
        JErr = inf;
        iterCnt = 0;
        while JErr > 0.02 && iterCnt <= 7
            % Gearbox effect
            [TrqCoeff,JEff] = PGBimpact(rS,rR,nP,mP,JS,JR,JC,JP);
            PXeffH = TrqCoeff.CTrqRC*NR./NC;
            PXeffL = TrqCoeff.CTrqSC*NS./NC;
            PXeffH_EPT = TrqCoeff.CTrqRC*NR_EPT./NC_EPT;
            PXeffL_EPT = TrqCoeff.CTrqSC*NS_EPT./NC_EPT;
            PXeffH_bins = [trapz(Nlps(idx1:idx2),PXeffH(idx1:idx2))/dNlps, trapz(Nlps(idx2:idx3),PXeffH(idx2:idx3))/dNlps, ...
                trapz(Nlps(idx3:idx4),PXeffH(idx3:idx4))/dNlps, trapz(Nlps(idx4:idx5),PXeffH(idx4:idx5))/dNlps, ...
                trapz(Nlps(idx5:idx6),PXeffH(idx5:idx6))/dNlps];
            % -- EM power
            if max(PXeffH) >= 0 
                MassE = 1000 + (mS + mR1 + mR2 + mC + mP); %slug
                disp('WARNING: The gearbox design does not favor the LPS across the entire engine power range and will not be consistent with the control algorithm')
            else
                % Powers, Torques, Energy (TEEM)
                [PemS_nom_MA, PemS_EPT_MA, PemS_TEEMa_MA, PemS_TEEMd_MA, ...
                  TrqemS_nom_MA, TrqemS_EPT_MA, TrqemS_TEEMa_MA, TrqemS_TEEMd_MA, ...
                  PemR_nom_MA, PemR_EPT_MA, PemR_TEEMa_MA, PemR_TEEMd_MA, ...
                  TrqemR_nom_MA, TrqemR_EPT_MA, TrqemR_TEEMa_MA, TrqemR_TEEMd_MA, ...
                  PemC_nom_MA, PemC_EPT_MA, PemC_TEEMa_MA, PemC_TEEMd_MA, ...
                  TrqemC_nom_MA, TrqemC_EPT_MA, TrqemC_TEEMa_MA, TrqemC_TEEMd_MA, ...
                  PemCoup_TEEMonlyA_MA, PemCoup_TEEMonlyD_MA, PemH_TEEMonly, EnergyUse_TEEMa] ...
                  = PGB_PwrTrq(HybConfig,GBConfig,NemH,NemL,NemCoup,NemH_EPT,NemCoup_EPT, ...
                    PXeffH,PXeffL,PXeffH_EPT,PXeffL_EPT,Plps_nom,Phps_nom,Plps_nom_EPT, ...
                    Phps_nom_EPT,PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,...
                    TEEMaWeights,TEEMdWeights,PXeffH_bins,transientDur,Psf,Trqsf);
                % Mass
                PemS_max = max(UseToggle.*[PemS_nom_MA PemS_EPT_MA PemS_TEEMa_MA PemS_TEEMd_MA]);%Psf;
                TrqemS_max = max(UseToggle.*[TrqemS_nom_MA TrqemS_EPT_MA TrqemS_TEEMa_MA TrqemS_TEEMd_MA]);%*Trqsf;
                PemR_max = max(UseToggle.*[PemR_nom_MA PemR_EPT_MA PemR_TEEMa_MA PemR_TEEMd_MA]);%*Psf;
                TrqemR_max = max(UseToggle.*[TrqemR_nom_MA TrqemR_EPT_MA TrqemR_TEEMa_MA TrqemR_TEEMd_MA]);%*Trqsf;
                PemC_max = max(UseToggle.*[PemC_nom_MA PemC_EPT_MA PemC_TEEMa_MA PemC_TEEMd_MA]);%*Psf;
                TrqemC_max = max(UseToggle.*[TrqemC_nom_MA TrqemC_EPT_MA TrqemC_TEEMa_MA TrqemC_TEEMd_MA]);%*Trqsf;
                numEMS = ceil(PemS_max/PmaxPerEM);
                MemS = numEMS*1.2866*0.0685218*0.3462*(TrqemS_max/max([numEMS, 1]))^0.7486;
                numEMR = ceil(PemR_max/PmaxPerEM);
                MemR = numEMR*1.2866*0.0685218*0.3462*(TrqemR_max/max([numEMR, 1]))^0.7486;
                numEMC = ceil(PemC_max/PmaxPerEM);
                MemC = numEMC*1.2866*0.0685218*0.3462*(TrqemC_max/max([numEMC, 1]))^0.7486;
                if HybConfig == 2
                    Mesd = 0; %boost leverages large existing battery (not included)
                else
                    Mesd = max([EnergyUse_TEEMa 0])/energyDenESD * Esf;
                end
                MpeS = PemS_max/PwrElecPwrDen;
                MpeR = PemR_max/PwrElecPwrDen;
                MpeC = PemC_max/PwrElecPwrDen;
                if UseToggle(3) ~= 1
                    Pesd = 0;
                else
                    Pesd = (EnergyUse_TEEMa/(transientDur*745.7*2.77778e-7));
                end
                if HybConfig == 2
                    MpeConv = max(Plps_nom)/PwrElecPwrDen; %boost leverages converter sized for boost (assumed to handle TEEM for short durations)
                else
                    MpeConv = Pesd/PwrElecPwrDen;
                end
                MassE = MemS + MemR + MemC + Mesd + MpeS + MpeR + MpeC + MpeConv; %mass of electrical system, slug
                % Inertia update
                JSold = JS;
                JRold = JR;
                JCold = JC;
                Drotor_EMS = 60*Vtip/(pi*NemS_max);
                JemS = (1/2)*MemS*(Drotor_EMS/2)^2;
                Drotor_EMR = 60*Vtip/(pi*NemR_max);
                JemR = (1/2)*MemR*(Drotor_EMR/2)^2;
                Drotor_EMC = 60*Vtip/(pi*NemC_max);
                JemC = (1/2)*MemC*(Drotor_EMC/2)^2;
                JS = JSo + Jlps/abs(NemSqNlps*NSqNemS)^2 + numEMS*JemS/abs(NSqNemS)^2; 
                JC = JCo + numEMC*JemC/abs(NCqNemC)^2;
                JR = JRo + Jhps/abs(NemRqNhps*NRqNemR)^2 + numEMR*JemR/abs(NRqNemR)^2; 
                JErr = max(abs(([JS JR JC]-[JSold JRold JCold])./[JSold JRold JCold]));
            end
            iterCnt = iterCnt + 1;
        end 
        JlpsEff = JEff.Sun*abs(NemSqNlps*NSqNemS)^2;
        JhpsEff = JEff.Ring*abs(NemRqNhps*NRqNemR)^2;
    
    elseif GBConfig == 4 % Sun - Coupling, Ring - HP, Carrier - LP
    
        % determine speeds of components
        % -- can infer the speed ratio of the components connected to the LPS
        %    and HPS relative to the LPS and HPS
        NemRqNhps = NemR_max/Nhps_max;
        NemCqNlps = NemC_max/Nlps_max;
        NemqNhps = NemRqNhps;
        NemqNlps = NemCqNlps;
        NRqNemR = NHPGBqNem;
        NCqNemC = NLPGBqNem;
        % -- planetary gearbox speeds
        NemR = Nhps*NemRqNhps;
        NR = NemR*NRqNemR;
        NemC = Nlps*NemCqNlps;
        NC = NemC*NCqNemC;
        NS = (1+1/rSqrR)*NC - (1/rSqrR)*NR;
        NSqNemS = max(abs(NS))/NemS_max;
        NemS = NS./NSqNemS;
        NemR_EPT = Nhps_EPT*NemRqNhps;
        NR_EPT = NemR_EPT*NRqNemR;
        NemC_EPT = Nlps_EPT*NemCqNlps;
        NC_EPT = NemC_EPT*NCqNemC;
        NS_EPT = (1+1/rSqrR)*NC_EPT - (1/rSqrR)*NR_EPT;
        NemS_EPT = NS_EPT./NSqNemS;
        % Cast EM speeds for PGB_PwrTrq.m
        NemH = NemR;
        NemL = NemC;
        NemCoup = NemS;
        NemH_EPT = NemR_EPT;
        NemCoup_EPT = NemS_EPT;
        % -- planetary gearbox inertias
        JS = JSo; % does not include EM(s) yet
        JC = JCo + Jlps/abs(NemCqNlps*NCqNemC)^2; % does not include EM(s) yet
        JR = JRo + Jhps/abs(NemRqNhps*NRqNemR)^2; % does not include EM(s) yet
        JErr = inf;
        iterCnt = 0;
        while JErr > 0.02 && iterCnt <= 7
            % Gearbox effect
            [TrqCoeff,JEff] = PGBimpact(rS,rR,nP,mP,JS,JR,JC,JP);
            PXeffH = TrqCoeff.CTrqRS*NR./NS;
            PXeffL = TrqCoeff.CTrqCS*NC./NS;
            PXeffH_EPT = TrqCoeff.CTrqRS*NR_EPT./NS_EPT;
            PXeffL_EPT = TrqCoeff.CTrqCS*NC_EPT./NS_EPT;
            PXeffH_bins = [trapz(Nlps(idx1:idx2),PXeffH(idx1:idx2))/dNlps, trapz(Nlps(idx2:idx3),PXeffH(idx2:idx3))/dNlps, ...
                trapz(Nlps(idx3:idx4),PXeffH(idx3:idx4))/dNlps, trapz(Nlps(idx4:idx5),PXeffH(idx4:idx5))/dNlps, ...
                trapz(Nlps(idx5:idx6),PXeffH(idx5:idx6))/dNlps];
            % -- EM power
            if max(PXeffH) >= 0 
                MassE = 1000 + (mS + mR1 + mR2 + mC + mP); %slug
                disp('WARNING: The gearbox design does not favor the LPS across the entire engine power range and will not be consistent with the control algorithm')
            else
                % Powers, Torques, Energy (TEEM)
                [PemS_nom_MA, PemS_EPT_MA, PemS_TEEMa_MA, PemS_TEEMd_MA, ...
                  TrqemS_nom_MA, TrqemS_EPT_MA, TrqemS_TEEMa_MA, TrqemS_TEEMd_MA, ...
                  PemR_nom_MA, PemR_EPT_MA, PemR_TEEMa_MA, PemR_TEEMd_MA, ...
                  TrqemR_nom_MA, TrqemR_EPT_MA, TrqemR_TEEMa_MA, TrqemR_TEEMd_MA, ...
                  PemC_nom_MA, PemC_EPT_MA, PemC_TEEMa_MA, PemC_TEEMd_MA, ...
                  TrqemC_nom_MA, TrqemC_EPT_MA, TrqemC_TEEMa_MA, TrqemC_TEEMd_MA, ...
                  PemCoup_TEEMonlyA_MA, PemCoup_TEEMonlyD_MA, PemH_TEEMonly, EnergyUse_TEEMa] ...
                  = PGB_PwrTrq(HybConfig,GBConfig,NemH,NemL,NemCoup,NemH_EPT,NemCoup_EPT, ...
                    PXeffH,PXeffL,PXeffH_EPT,PXeffL_EPT,Plps_nom,Phps_nom,Plps_nom_EPT, ...
                    Phps_nom_EPT,PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,...
                    TEEMaWeights,TEEMdWeights,PXeffH_bins,transientDur,Psf,Trqsf);
                % Mass
                PemS_max = max(UseToggle.*[PemS_nom_MA PemS_EPT_MA PemS_TEEMa_MA PemS_TEEMd_MA]);%Psf;
                TrqemS_max = max(UseToggle.*[TrqemS_nom_MA TrqemS_EPT_MA TrqemS_TEEMa_MA TrqemS_TEEMd_MA]);%*Trqsf;
                PemR_max = max(UseToggle.*[PemR_nom_MA PemR_EPT_MA PemR_TEEMa_MA PemR_TEEMd_MA]);%*Psf;
                TrqemR_max = max(UseToggle.*[TrqemR_nom_MA TrqemR_EPT_MA TrqemR_TEEMa_MA TrqemR_TEEMd_MA]);%*Trqsf;
                PemC_max = max(UseToggle.*[PemC_nom_MA PemC_EPT_MA PemC_TEEMa_MA PemC_TEEMd_MA]);%*Psf;
                TrqemC_max = max(UseToggle.*[TrqemC_nom_MA TrqemC_EPT_MA TrqemC_TEEMa_MA TrqemC_TEEMd_MA]);%*Trqsf;
                numEMS = ceil(PemS_max/PmaxPerEM);
                MemS = numEMS*1.2866*0.0685218*0.3462*(TrqemS_max/max([numEMS, 1]))^0.7486;
                numEMR = ceil(PemR_max/PmaxPerEM);
                MemR = numEMR*1.2866*0.0685218*0.3462*(TrqemR_max/max([numEMR, 1]))^0.7486;
                numEMC = ceil(PemC_max/PmaxPerEM);
                MemC = numEMC*1.2866*0.0685218*0.3462*(TrqemC_max/max([numEMC, 1]))^0.7486;
                if HybConfig == 2
                    Mesd = 0; %boost leverages large existing battery (not included)
                else
                    Mesd = max([EnergyUse_TEEMa 0])/energyDenESD * Esf;
                end
                MpeS = PemS_max/PwrElecPwrDen;
                MpeR = PemR_max/PwrElecPwrDen;
                MpeC = PemC_max/PwrElecPwrDen;
                if UseToggle(3) ~= 1
                    Pesd = 0;
                else
                    Pesd = (EnergyUse_TEEMa/(transientDur*745.7*2.77778e-7));
                end
                if HybConfig == 2
                    MpeConv = max(Plps_nom)/PwrElecPwrDen; %boost leverages converter sized for boost (assumed to handle TEEM for short durations)
                else
                    MpeConv = Pesd/PwrElecPwrDen;
                end
                MassE = MemS + MemR + MemC + Mesd + MpeS + MpeR + MpeC + MpeConv; %mass of electrical system, slug
                % Inertia update
                JSold = JS;
                JRold = JR;
                JCold = JC;
                Drotor_EMS = 60*Vtip/(pi*NemS_max);
                JemS = (1/2)*MemS*(Drotor_EMS/2)^2;
                Drotor_EMR = 60*Vtip/(pi*NemR_max);
                JemR = (1/2)*MemR*(Drotor_EMR/2)^2;
                Drotor_EMC = 60*Vtip/(pi*NemC_max);
                JemC = (1/2)*MemC*(Drotor_EMC/2)^2;
                JS = JSo + numEMS*JemS/abs(NSqNemS)^2; 
                JC = JCo + Jlps/abs(NemCqNlps*NCqNemC)^2 + numEMC*JemC/abs(NCqNemC)^2;
                JR = JRo + Jhps/abs(NemRqNhps*NRqNemR)^2 + numEMR*JemR/abs(NRqNemR)^2; 
                JErr = max(abs(([JS JR JC]-[JSold JRold JCold])./[JSold JRold JCold]));
            end
            iterCnt = iterCnt + 1;
        end 
        JlpsEff = JEff.Carrier*abs(NemCqNlps*NCqNemC)^2;
        JhpsEff = JEff.Ring*abs(NemRqNhps*NRqNemR)^2;
    
    elseif GBConfig == 5 % Sun - LP, Ring - Coupling, Carrier - HP
    
        % determine speeds of components
        % -- can infer the speed ratio of the components connected to the LPS
        %    and HPS relative to the LPS and HPS
        NemCqNhps = NemC_max/Nhps_max;
        NemSqNlps = NemS_max/Nlps_max;
        NemqNhps = NemCqNhps;
        NemqNlps = NemSqNlps;
        NCqNemC = NHPGBqNem;
        NSqNemS = NLPGBqNem;
        % -- planetary gearbox speeds
        NemC = Nhps*NemCqNhps;
        NC = NemC*NCqNemC;
        NemS = Nlps*NemSqNlps;
        NS = NemS*NSqNemS;
        NR = (rSqrR+1)*NC-rSqrR*NS;
        NRqNemR = max(abs(NR))/NemR_max;
        NemR = NR./NRqNemR;
        NemC_EPT = Nhps_EPT*NemCqNhps;
        NC_EPT = NemC_EPT*NCqNemC;
        NemS_EPT = Nlps_EPT*NemSqNlps;
        NS_EPT = NemS_EPT*NSqNemS;
        NR_EPT = (rSqrR+1)*NC_EPT-rSqrR*NS_EPT;
        NemR_EPT = NR_EPT./NRqNemR;
        % Cast EM speeds for PGB_PwrTrq.m
        NemH = NemC;
        NemL = NemS;
        NemCoup = NemR;
        NemH_EPT = NemC_EPT;
        NemCoup_EPT = NemR_EPT;
        % -- planetary gearbox inertias
        JS = JSo + Jlps/abs(NemSqNlps*NSqNemS)^2; % does not include EM(s) yet
        JC = JCo + Jhps/abs(NemCqNhps*NCqNemC)^2; % does not include EM(s) yet
        JR = JRo; % does not include EM(s) yet
        JErr = inf;
        iterCnt = 0;
        while JErr > 0.02 && iterCnt <= 7
            % Gearbox effect
            [TrqCoeff,JEff] = PGBimpact(rS,rR,nP,mP,JS,JR,JC,JP);
            PXeffH = TrqCoeff.CTrqCR*NC./NR;
            PXeffL = TrqCoeff.CTrqSR*NS./NR;
            PXeffH_EPT = TrqCoeff.CTrqSR*NC_EPT./NR_EPT;
            PXeffL_EPT = TrqCoeff.CTrqCR*NS_EPT./NR_EPT;
            PXeffH_bins = [trapz(Nlps(idx1:idx2),PXeffH(idx1:idx2))/dNlps, trapz(Nlps(idx2:idx3),PXeffH(idx2:idx3))/dNlps, ...
                trapz(Nlps(idx3:idx4),PXeffH(idx3:idx4))/dNlps, trapz(Nlps(idx4:idx5),PXeffH(idx4:idx5))/dNlps, ...
                trapz(Nlps(idx5:idx6),PXeffH(idx5:idx6))/dNlps];
            % -- EM power
            if max(PXeffH) >= 0 
                MassE = 1000 + (mS + mR1 + mR2 + mC + mP); %slug
                disp('WARNING: The gearbox design does not favor the LPS across the entire engine power range and will not be consistent with the control algorithm')
            else
                % Powers, Torques, Energy (TEEM)
                [PemS_nom_MA, PemS_EPT_MA, PemS_TEEMa_MA, PemS_TEEMd_MA, ...
                  TrqemS_nom_MA, TrqemS_EPT_MA, TrqemS_TEEMa_MA, TrqemS_TEEMd_MA, ...
                  PemR_nom_MA, PemR_EPT_MA, PemR_TEEMa_MA, PemR_TEEMd_MA, ...
                  TrqemR_nom_MA, TrqemR_EPT_MA, TrqemR_TEEMa_MA, TrqemR_TEEMd_MA, ...
                  PemC_nom_MA, PemC_EPT_MA, PemC_TEEMa_MA, PemC_TEEMd_MA, ...
                  TrqemC_nom_MA, TrqemC_EPT_MA, TrqemC_TEEMa_MA, TrqemC_TEEMd_MA, ...
                  PemCoup_TEEMonlyA_MA, PemCoup_TEEMonlyD_MA, PemH_TEEMonly, EnergyUse_TEEMa] ...
                  = PGB_PwrTrq(HybConfig,GBConfig,NemH,NemL,NemCoup,NemH_EPT,NemCoup_EPT, ...
                    PXeffH,PXeffL,PXeffH_EPT,PXeffL_EPT,Plps_nom,Phps_nom,Plps_nom_EPT, ...
                    Phps_nom_EPT,PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,...
                    TEEMaWeights,TEEMdWeights,PXeffH_bins,transientDur,Psf,Trqsf);
                % Mass
                PemS_max = max(UseToggle.*[PemS_nom_MA PemS_EPT_MA PemS_TEEMa_MA PemS_TEEMd_MA]);%Psf;
                TrqemS_max = max(UseToggle.*[TrqemS_nom_MA TrqemS_EPT_MA TrqemS_TEEMa_MA TrqemS_TEEMd_MA]);%*Trqsf;
                PemR_max = max(UseToggle.*[PemR_nom_MA PemR_EPT_MA PemR_TEEMa_MA PemR_TEEMd_MA]);%*Psf;
                TrqemR_max = max(UseToggle.*[TrqemR_nom_MA TrqemR_EPT_MA TrqemR_TEEMa_MA TrqemR_TEEMd_MA]);%*Trqsf;
                PemC_max = max(UseToggle.*[PemC_nom_MA PemC_EPT_MA PemC_TEEMa_MA PemC_TEEMd_MA]);%*Psf;
                TrqemC_max = max(UseToggle.*[TrqemC_nom_MA TrqemC_EPT_MA TrqemC_TEEMa_MA TrqemC_TEEMd_MA]);%*Trqsf;
                numEMS = ceil(PemS_max/PmaxPerEM);
                MemS = numEMS*1.2866*0.0685218*0.3462*(TrqemS_max/max([numEMS, 1]))^0.7486;
                numEMR = ceil(PemR_max/PmaxPerEM);
                MemR = numEMR*1.2866*0.0685218*0.3462*(TrqemR_max/max([numEMR, 1]))^0.7486;
                numEMC = ceil(PemC_max/PmaxPerEM);
                MemC = numEMC*1.2866*0.0685218*0.3462*(TrqemC_max/max([numEMC, 1]))^0.7486;
                if HybConfig == 2
                    Mesd = 0; %boost leverages large existing battery (not included)
                else
                    Mesd = max([EnergyUse_TEEMa 0])/energyDenESD * Esf;
                end
                MpeS = PemS_max/PwrElecPwrDen;
                MpeR = PemR_max/PwrElecPwrDen;
                MpeC = PemC_max/PwrElecPwrDen;
                if UseToggle(3) ~= 1
                    Pesd = 0;
                else
                    Pesd = (EnergyUse_TEEMa/(transientDur*745.7*2.77778e-7));
                end
                if HybConfig == 2
                    MpeConv = max(Plps_nom)/PwrElecPwrDen; %boost leverages converter sized for boost (assumed to handle TEEM for short durations)
                else
                    MpeConv = Pesd/PwrElecPwrDen;
                end
                MassE = MemS + MemR + MemC + Mesd + MpeS + MpeR + MpeC + MpeConv; %mass of electrical system, slug
                % Inertia update
                JSold = JS;
                JRold = JR;
                JCold = JC;
                Drotor_EMS = 60*Vtip/(pi*NemS_max);
                JemS = (1/2)*MemS*(Drotor_EMS/2)^2;
                Drotor_EMR = 60*Vtip/(pi*NemR_max);
                JemR = (1/2)*MemR*(Drotor_EMR/2)^2;
                Drotor_EMC = 60*Vtip/(pi*NemC_max);
                JemC = (1/2)*MemC*(Drotor_EMC/2)^2;
                JS = JSo + Jlps/abs(NemSqNlps*NSqNemS)^2 + numEMS*JemS/abs(NSqNemS)^2; 
                JC = JCo + Jhps/abs(NemCqNhps*NCqNemC)^2 + numEMC*JemC/abs(NCqNemC)^2;
                JR = JRo + numEMR*JemR/abs(NRqNemR)^2; 
                JErr = max(abs(([JS JR JC]-[JSold JRold JCold])./[JSold JRold JCold]));
            end
            iterCnt = iterCnt + 1;
        end 
        JlpsEff = JEff.Sun*abs(NemSqNlps*NSqNemS)^2;
        JhpsEff = JEff.Carrier*abs(NemCqNhps*NCqNemC)^2;
    
    elseif GBConfig == 6 % Sun - Coupling, Ring - LP, Carrier - HP
    
        % determine speeds of components
        % -- can infer the speed ratio of the components connected to the LPS
        %    and HPS relative to the LPS and HPS
        NemCqNhps = NemC_max/Nhps_max;
        NemRqNlps = NemR_max/Nlps_max;
        NemqNhps = NemCqNhps;
        NemqNlps = NemRqNlps;
        NCqNemC = NHPGBqNem;
        NRqNemR = NLPGBqNem;
        % -- planetary gearbox speeds
        NemC = Nhps*NemCqNhps;
        NC = NemC*NCqNemC;
        NemR = Nlps*NemRqNlps;
        NR = NemR*NRqNemR;
        NS = (1+1/rSqrR)*NC - (1/rSqrR)*NR;
        NSqNemS = max(abs(NS))/NemS_max;
        NemS = NS./NSqNemS;
        NemC_EPT = Nhps_EPT*NemCqNhps;
        NC_EPT = NemC_EPT*NCqNemC;
        NemR_EPT = Nlps_EPT*NemRqNlps;
        NR_EPT = NemR_EPT*NRqNemR;
        NS_EPT = (1+1/rSqrR)*NC_EPT - (1/rSqrR)*NR_EPT;
        NemS_EPT = NS_EPT./NSqNemS;
        % Cast EM speeds for PGB_PwrTrq.m
        NemH = NemC;
        NemL = NemR;
        NemCoup = NemS;
        NemH_EPT = NemC_EPT;
        NemCoup_EPT = NemS_EPT;
        % -- planetary gearbox inertias
        JS = JSo; % does not include EM(s) yet
        JC = JCo + Jhps/abs(NemCqNhps*NCqNemC)^2; % does not include EM(s) yet
        JR = JRo + Jlps/abs(NemRqNlps*NRqNemR)^2; % does not include EM(s) yet
        JErr = inf;
        iterCnt = 0;
        while JErr > 0.02 && iterCnt <= 7
            % Gearbox effect
            [TrqCoeff,JEff] = PGBimpact(rS,rR,nP,mP,JS,JR,JC,JP);
            PXeffH = TrqCoeff.CTrqCS*NC./NS;
            PXeffL = TrqCoeff.CTrqRS*NR./NS;
            PXeffH_EPT = TrqCoeff.CTrqCS*NC_EPT./NS_EPT;
            PXeffL_EPT = TrqCoeff.CTrqRS*NR_EPT./NS_EPT;
            PXeffH_bins = [trapz(Nlps(idx1:idx2),PXeffH(idx1:idx2))/dNlps, trapz(Nlps(idx2:idx3),PXeffH(idx2:idx3))/dNlps, ...
                trapz(Nlps(idx3:idx4),PXeffH(idx3:idx4))/dNlps, trapz(Nlps(idx4:idx5),PXeffH(idx4:idx5))/dNlps, ...
                trapz(Nlps(idx5:idx6),PXeffH(idx5:idx6))/dNlps];
            % -- EM power
            if max(PXeffH) >= 0 
                MassE = 1000 + (mS + mR1 + mR2 + mC + mP); %slug
                disp('WARNING: The gearbox design does not favor the LPS across the entire engine power range and will not be consistent with the control algorithm')
            else
                % Powers, Torques, Energy (TEEM)
                [PemS_nom_MA, PemS_EPT_MA, PemS_TEEMa_MA, PemS_TEEMd_MA, ...
                  TrqemS_nom_MA, TrqemS_EPT_MA, TrqemS_TEEMa_MA, TrqemS_TEEMd_MA, ...
                  PemR_nom_MA, PemR_EPT_MA, PemR_TEEMa_MA, PemR_TEEMd_MA, ...
                  TrqemR_nom_MA, TrqemR_EPT_MA, TrqemR_TEEMa_MA, TrqemR_TEEMd_MA, ...
                  PemC_nom_MA, PemC_EPT_MA, PemC_TEEMa_MA, PemC_TEEMd_MA, ...
                  TrqemC_nom_MA, TrqemC_EPT_MA, TrqemC_TEEMa_MA, TrqemC_TEEMd_MA, ...
                  PemCoup_TEEMonlyA_MA, PemCoup_TEEMonlyD_MA, PemH_TEEMonly, EnergyUse_TEEMa] ...
                  = PGB_PwrTrq(HybConfig,GBConfig,NemH,NemL,NemCoup,NemH_EPT,NemCoup_EPT, ...
                    PXeffH,PXeffL,PXeffH_EPT,PXeffL_EPT,Plps_nom,Phps_nom,Plps_nom_EPT, ...
                    Phps_nom_EPT,PEx_ACsys,PX_EPT,Pin_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,...
                    TEEMaWeights,TEEMdWeights,PXeffH_bins,transientDur,Psf,Trqsf);
                % Mass
                PemS_max = max(UseToggle.*[PemS_nom_MA PemS_EPT_MA PemS_TEEMa_MA PemS_TEEMd_MA]);%Psf;
                TrqemS_max = max(UseToggle.*[TrqemS_nom_MA TrqemS_EPT_MA TrqemS_TEEMa_MA TrqemS_TEEMd_MA]);%*Trqsf;
                PemR_max = max(UseToggle.*[PemR_nom_MA PemR_EPT_MA PemR_TEEMa_MA PemR_TEEMd_MA]);%*Psf;
                TrqemR_max = max(UseToggle.*[TrqemR_nom_MA TrqemR_EPT_MA TrqemR_TEEMa_MA TrqemR_TEEMd_MA]);%*Trqsf;
                PemC_max = max(UseToggle.*[PemC_nom_MA PemC_EPT_MA PemC_TEEMa_MA PemC_TEEMd_MA]);%*Psf;
                TrqemC_max = max(UseToggle.*[TrqemC_nom_MA TrqemC_EPT_MA TrqemC_TEEMa_MA TrqemC_TEEMd_MA]);%*Trqsf;
                numEMS = ceil(PemS_max/PmaxPerEM);
                MemS = numEMS*1.2866*0.0685218*0.3462*(TrqemS_max/max([numEMS, 1]))^0.7486;
                numEMR = ceil(PemR_max/PmaxPerEM);
                MemR = numEMR*1.2866*0.0685218*0.3462*(TrqemR_max/max([numEMR, 1]))^0.7486;
                numEMC = ceil(PemC_max/PmaxPerEM);
                MemC = numEMC*1.2866*0.0685218*0.3462*(TrqemC_max/max([numEMC, 1]))^0.7486;
                if HybConfig == 2
                    Mesd = 0; %boost leverages large existing battery (not included)
                else
                    Mesd = max([EnergyUse_TEEMa 0])/energyDenESD * Esf;
                end
                MpeS = PemS_max/PwrElecPwrDen;
                MpeR = PemR_max/PwrElecPwrDen;
                MpeC = PemC_max/PwrElecPwrDen;
                if UseToggle(3) ~= 1
                    Pesd = 0;
                else
                    Pesd = (EnergyUse_TEEMa/(transientDur*745.7*2.77778e-7));
                end
                if HybConfig == 2
                    MpeConv = max(Plps_nom)/PwrElecPwrDen; %boost leverages converter sized for boost (assumed to handle TEEM for short durations)
                else
                    MpeConv = Pesd/PwrElecPwrDen;
                end
                MassE = MemS + MemR + MemC + Mesd + MpeS + MpeR + MpeC + MpeConv; %mass of electrical system, slug
                % Inertia update
                JSold = JS;
                JRold = JR;
                JCold = JC;
                Drotor_EMS = 60*Vtip/(pi*NemS_max);
                JemS = (1/2)*MemS*(Drotor_EMS/2)^2;
                Drotor_EMR = 60*Vtip/(pi*NemR_max);
                JemR = (1/2)*MemR*(Drotor_EMR/2)^2;
                Drotor_EMC = 60*Vtip/(pi*NemC_max);
                JemC = (1/2)*MemC*(Drotor_EMC/2)^2;
                JS = JSo + numEMS*JemS/abs(NSqNemS)^2; 
                JC = JCo + Jhps/abs(NemCqNhps*NCqNemC)^2 + numEMC*JemC/abs(NCqNemC)^2;
                JR = JRo + Jlps/abs(NemRqNlps*NRqNemR)^2 + numEMR*JemR/abs(NRqNemR)^2; 
                JErr = max(abs(([JS JR JC]-[JSold JRold JCold])./[JSold JRold JCold]));
            end
            iterCnt = iterCnt + 1;
        end 
        JlpsEff = JEff.Ring*abs(NemRqNlps*NRqNemR)^2;
        JhpsEff = JEff.Carrier*abs(NemCqNhps*NCqNemC)^2;

    else
        fprintf('Error: Inproper definition of Config. It must have an integer >= 1 and <= 6')
    end

    %iteration check
    iter = iter + 1;
    k = k*0.95;
    if iter >= iter_max
        disp('Unable to converge in Nem while-loop')
        break;
    end

end

if iterCnt >= 7

    Mass = 1000 + (mS + mR1 + mR2 + mC + mP); %slug
    disp('WARNING: The gearbox design could not be evaluated')

    NP = (rS/rP)*(NC-NS);

    info.NemSmax = NemS_max;
    info.NemRmax = NemR_max;
    info.NemCmax = NemC_max;
    info.NSrange = [min(NS) max(NS)];
    info.NRrange = [min(NR) max(NR)];
    info.NCrange = [min(NC) max(NC)];
    info.NPrange = [min(NP) max(NP)];
    info.NemSrange = [min(NemS) max(NemS)];
    info.NemRrange = [min(NemR) max(NemR)];
    info.NemCrange = [min(NemC) max(NemC)];
    info.NSqNemS = NSqNemS;
    info.NRqNemR = NRqNemR;
    info.NCqNemC = NCqNemC;
    info.NemqNhps = NemqNhps;
    info.NemqNlps = NemqNlps;
    info.PemS = PemS_max;
    info.PemR = PemR_max;
    info.PemC = PemC_max;
    info.Pesd = 0; %Pesd;
    info.TrqemS = 0; %TrqemS_max;
    info.TrqemR = 0; %TrqemR_max;
    info.TrqemC = 0; %TrqemC_max;
    info.PemCoup_TEEMonlyA = 0; %PemCoup_TEEMonlyA_MA;
    info.PemCoup_TEEMonlyD = 0; %PemCoup_TEEMonlyD_MA;
    info.PemH_TEEMonly = 0; %PemH_TEEMonly;
    info.MemS = 0; %MemS;
    info.MemR = 0; %MemR;
    info.MemC = 0; %MemC;
    info.Mesd = 0; %Mesd;
    info.MpeS = 0; %MpeS;
    info.MpeR = 0; %MpeR;
    info.MpeC = 0; %MpeC;
    info.MpeConv = 0; %MpeConv;
    info.Mps = 0; %MassE;
    % info.MS = mS;
    % info.MR = mR1 + mR2;
    % info.MC = mC;
    % info.MP = mP; 
    % info.Mgb = mS + mR1 + mR2 + mC + nP*mP; %very rough estimate
    info.JSo = JSo;
    info.JRo = JRo;
    info.JCo = JCo;
    info.JemS = 0; %JemS;
    info.JemR = 0; %JemR;
    info.JemC = 0; %JemC;
    info.JS = 0; %JS;
    info.JR = 0; %JR;
    info.JC = 0; %JC;
    info.JP = 0; %JP;
    info.nP = nP;
    info.JSEff = 0; %JEff.Sun;
    info.JREff = 0; %JEff.Ring;
    info.JCEff = 0; %JEff.Carrier;
    info.JlpsEff = 0; %JlpsEff;
    info.JhpsEff = 0; %JhpsEff;
    info.MassDriver_emS = 0;
    info.MassDriver_emR = 0;
    info.MassDriver_emC = 0;
    info.MassDriver_peS = 0;
    info.MassDriver_peR = 0;
    info.MassDriver_peC = 0;
    info.EnergyStorage = 0; %EnergyUse_TEEMa*Esf;
    info.TrqCoeff = TrqCoeff;
    info.PXeff_H = PXeffH;
    info.PXeff_L = PXeffL;
    info.Mts = 0;
    info.Mbg = 0;
    info.Msg = 0;
    info.Mpgb = 0;
    info.Mmech = 0;
    info.Mtot = Mass;
    info.rSqrR = rSqrR;
    info.NLPGBqNem = NLPGBqNem;
    info.NtsqNhps = NtsqNhps;
    info.NtsqNlps = NtsqNlps;
    if iter >= iter_max
        info.NemMaxLoopPass = 0;
    else
        info.NemMaxLoopPass = 1;
    end

else

    % penalize speed contraint violations
    if max(abs(NS)) > NS_max
        MassPenaltyS = (max(abs(NS))-NS_max)*0.1; %artificially adds 0.1 slug per rpm over the limit (1 slug / 10 rpm)
    else
        MassPenaltyS = 0;
    end
    if max(abs(NR)) > NR_max
        MassPenaltyR = (max(abs(NR))-NR_max)*0.1; %artificially adds 0.1 slug per rpm over the limit (1 slug / 10 rpm)
    else
        MassPenaltyR = 0;
    end
    if max(abs(NC)) > NC_max
        MassPenaltyC = (max(abs(NC))-NC_max)*0.1; %artificially adds 0.1 slug per rpm over the limit (1 slug / 10 rpm)
    else
        MassPenaltyC = 0;
    end
    NP = (rS/rP)*(NC-NS);
    if max(abs(NP)) > NP_max
        MassPenaltyP = (max(abs(NP))-NP_max)*0.1; %artificially adds 0.1 slug per rpm over the limit (1 slug / 10 rpm)
    else
        MassPenaltyP = 0;
    end
    MassE = MassE + MassPenaltyS + MassPenaltyR + MassPenaltyC + MassPenaltyP;

    % Mechanical system weight approximation
    %   Approach: approximate torques from the engine shaft requirements
    %   and applie saftey factors. This will be good working back as far as
    %   the PGB components. The PGB weight is estimated based on ring gear
    %   torque which will fall out of back propagating the engine shaft
    %   requirements if the ring gear is interfaced with one of the shaft.
    %   Otherwise it is the coupling component and the maximum torque is
    %   approximated as the maximum torque of the coupling component
    %   electric machine.

    % Properties
    E = 200*145038; %Young's Modulus, GPa->psi
    fv = 0.9; %1-v^2 where v is Poisson's Ratio
    sigy = 290*145.038; %yield stress, MPa->psi
    rho = 7900*1.12287e-6; %density, kg/m3->slug/in3
    L = 60; %length of shaft, in
    minSWT = 0.125; %minimum shaft wall thickness, in
    alpha = 10; %gear tooth contact stress angle
    sigc = 1380*145.038; %contact strength, MPa->psi
    Xsg = 1; %fill factor for spur gears
    Xc = 0.75; %fill factor for carrier
    Xp = 1; %fill factor for planets gears
    Xr = 0.3225; %fill factor for ring gear
    Xs = 1; %fill factor for sun gear

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
    NemLqNts = abs(NemqNlps/NtsqNlps);
    NemHqNts = abs(NemqNhps/NtsqNhps);
    
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
    Trqmax = max(abs([Trqhps_max/NtsqNhps, Trqhps_max/NemqNhps]));
    Q_emh = (Trqmax/5252.113)*(nbg_emh+1)^3/nbg_emh^2;
    Mbgemh = (55/32.174)*Q_emh^0.91;
    % --LP shaft - spool side
    nbg_l = min(NtsqNlps,1/NtsqNlps);
    Trqmax = max(abs([Trqlps_max/NtsqNlps, Trqlps_max]));
    Q_l = (Trqmax/5252.113)*(nbg_l+1)^3/nbg_l^2;
    Mbgl = (55/32.174)*Q_l^0.91;
    % --LP shaft - EM side
    nbg_eml = min(NemLqNts,1/NemLqNts);
    Trqmax = max(abs([Trqlps_max/NtsqNlps, Trqlps_max/NemqNlps]));
    Q_eml = (Trqmax/5252.113)*(nbg_eml+1)^3/nbg_eml^2;
    Mbgeml = (55/32.174)*Q_eml^0.91;
    % --Total
    Mbg = Mbgemh + Mbgh + Mbgeml + Mbgl; %slug

    % Spur Gears
    % --LP shaft to PGB
    NsgL = abs(NLPGBqNem);
    TrqsgL2_max = Trqlps_max/NemqNlps; %Trq on EM side
    Vsg2L = (E*144/fv)*(TrqsgL2_max/sind(2*alpha))*(GBsf/(sigc*144))^2*(1+NsgL);
    Vsg1L = Vsg2L/NsgL^2;
    if NLPGBqNem > 0 % need to switch direction with another spur gear
        if NLPGBqNem >= 1
            VsgSD = min([Vsg1L Vsg2L]);
        else
            VsgSD = max([Vsg1L Vsg2L]);
        end
    else
        VsgSD = 0;
    end
    MsgL = GBsf*(rho*12^3)*(Vsg1L + Vsg2L + VsgSD)*Xsg;
    % --Coupling EM to PGB
    if GBConfig == 1 || GBConfig == 3
        NsgC = abs(NCqNemC);
        TrqsgC2_max = TrqemC_max; %Trq on EM side
    elseif GBConfig == 2 || GBConfig == 5
        NsgC = abs(NRqNemR);
        TrqsgC2_max = TrqemR_max; %Trq on EM side
    else
        NsgC = abs(NSqNemS);
        TrqsgC2_max = TrqemS_max; %Trq on EM side
    end
    Vsg2C = (E*144/fv)*(TrqsgC2_max/sind(2*alpha))*(GBsf/(sigc*144))^2*(1+NsgC);
    Vsg1C = Vsg2C/NsgC^2;
    MsgC = GBsf*(rho*12^3)*(Vsg1C + Vsg2C)*Xsg;
    % --Total
    Msg = MsgL + MsgC;

    % Planetary Gearbox 
    if GBConfig == 1 || GBConfig == 6 % ring on LPS
        NRqNlps = abs(NRqNemR*NemqNlps);
        TrqR_max = Trqlps_max/abs(NRqNlps);
    elseif GBConfig == 2 || GBConfig == 5 % ring on coupling
        TrqR_max = TrqsgC2_max/NsgC;
    else
        NRqNhps = abs(NRqNemR*NemqNhps);
        TrqR_max = Trqhps_max/abs(NRqNhps);
    end
    ns = rSqrR;
    np = (1-ns)/2;
    K = (rho*12^3)*(E*144/fv)*(TrqR_max/sind(2*alpha))*(GBsf/(sigc*144))^2*(1/nP)*(1+ns)/(2*ns);
    % --Carrier
    Mc = Xc*K*np;
    % --Planet
    Mp = Xp*K/np;
    % --Ring
    Mr = Xr*K/np;
    % --Sun
    Ms = Xs*K*np;
    % --Total
    Mpgb = Mc + nP*Mp + Mr + Ms;

    % Total Mech System Mass
    Mmech = Mts + Mbg + Msg + Mpgb; %slug

    % Overall System Mass
    Mass = MassE + Mmech; %slug

    % END of Mechanical system weight approximation

    % pick off mass drivers for the motors and power electronics
    %   1 - nominial use 
    %   2 - EPT
    %   3 - TEEM-accel
    %   4 - TEEM-decel
    [dumbie,MassDriver_emS] = max(UseToggle.*[TrqemS_nom_MA TrqemS_EPT_MA TrqemS_TEEMa_MA TrqemS_TEEMd_MA]);
    [dumbie,MassDriver_emR] = max(UseToggle.*[TrqemR_nom_MA TrqemR_EPT_MA TrqemR_TEEMa_MA TrqemR_TEEMd_MA]);
    [dumbie,MassDriver_emC] = max(UseToggle.*[TrqemC_nom_MA TrqemC_EPT_MA TrqemC_TEEMa_MA TrqemC_TEEMd_MA]);
    [dumbie,MassDriver_peS] = max(UseToggle.*[PemS_nom_MA PemS_EPT_MA PemS_TEEMa_MA PemS_TEEMd_MA]);
    [dumbie,MassDriver_peR] = max(UseToggle.*[PemR_nom_MA PemR_EPT_MA PemR_TEEMa_MA PemR_TEEMd_MA]);
    [dumbie,MassDriver_peC] = max(UseToggle.*[PemC_nom_MA PemC_EPT_MA PemC_TEEMa_MA PemC_TEEMd_MA]);
    
    info.NemSmax = NemS_max;
    info.NemRmax = NemR_max;
    info.NemCmax = NemC_max;
    info.NSrange = [min(NS) max(NS)];
    info.NRrange = [min(NR) max(NR)];
    info.NCrange = [min(NC) max(NC)];
    info.NPrange = [min(NP) max(NP)];
    info.NemSrange = [min(NemS) max(NemS)];
    info.NemRrange = [min(NemR) max(NemR)];
    info.NemCrange = [min(NemC) max(NemC)];
    info.NSqNemS = NSqNemS;
    info.NRqNemR = NRqNemR;
    info.NCqNemC = NCqNemC;
    info.NemqNhps = NemqNhps;
    info.NemqNlps = NemqNlps;
    info.PemS = PemS_max;
    info.PemR = PemR_max;
    info.PemC = PemC_max;
    info.Pesd = Pesd;
    info.TrqemS = TrqemS_max;
    info.TrqemR = TrqemR_max;
    info.TrqemC = TrqemC_max;
    info.PemCoup_TEEMonlyA = PemCoup_TEEMonlyA_MA;
    info.PemCoup_TEEMonlyD = PemCoup_TEEMonlyD_MA;
    info.PemH_TEEMonly = PemH_TEEMonly;
    info.MemS = MemS;
    info.MemR = MemR;
    info.MemC = MemC;
    info.Mesd = Mesd;
    info.MpeS = MpeS;
    info.MpeR = MpeR;
    info.MpeC = MpeC;
    info.MpeConv = MpeConv;
    info.Mps = MassE;
    % info.MS = mS;
    % info.MR = mR1 + mR2;
    % info.MC = mC;
    % info.MP = mP; 
    % info.Mgb = mS + mR1 + mR2 + mC + nP*mP; %very rough estimate
    info.JSo = JSo;
    info.JRo = JRo;
    info.JCo = JCo;
    info.JemS = JemS;
    info.JemR = JemR;
    info.JemC = JemC;
    info.JS = JS;
    info.JR = JR;
    info.JC = JC;
    info.JP = JP;
    info.nP = nP;
    info.mP = Mp;
    info.JSEff = JEff.Sun;
    info.JREff = JEff.Ring;
    info.JCEff = JEff.Carrier;
    info.JlpsEff = JlpsEff;
    info.JhpsEff = JhpsEff;
    info.MassDriver_emS = MassDriver_emS;
    info.MassDriver_emR = MassDriver_emR;
    info.MassDriver_emC = MassDriver_emC;
    info.MassDriver_peS = MassDriver_peS;
    info.MassDriver_peR = MassDriver_peR;
    info.MassDriver_peC = MassDriver_peC;
    info.EnergyStorage = EnergyUse_TEEMa*Esf;
    info.TrqCoeff = TrqCoeff;
    info.PXeff_H = PXeffH;
    info.PXeff_L = PXeffL;
    info.Mts = Mts;
    info.Mbg = Mbg;
    info.Msg = Msg;
    info.Mpgb = Mpgb;
    info.Mmech = Mmech;
    info.Mtot = Mass;
    info.rSqrR = rSqrR;
    info.NLPGBqNem = NLPGBqNem;
    info.NtsqNhps = NtsqNhps;
    info.NtsqNlps = NtsqNlps;
    if iter >= iter_max
        info.NemMaxLoopPass = 0;
    else
        info.NemMaxLoopPass = 1;
    end

end

end