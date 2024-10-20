function [PemS_nom_MA, PemS_EPT_MA, PemS_TEEMa_MA, PemS_TEEMd_MA, ...
          TrqemS_nom_MA, TrqemS_EPT_MA, TrqemS_TEEMa_MA, TrqemS_TEEMd_MA, ...
          PemR_nom_MA, PemR_EPT_MA, PemR_TEEMa_MA, PemR_TEEMd_MA, ...
          TrqemR_nom_MA, TrqemR_EPT_MA, TrqemR_TEEMa_MA, TrqemR_TEEMd_MA, ...
          PemC_nom_MA, PemC_EPT_MA, PemC_TEEMa_MA, PemC_TEEMd_MA, ...
          TrqemC_nom_MA, TrqemC_EPT_MA, TrqemC_TEEMa_MA, TrqemC_TEEMd_MA, ...
          PemCoup_TEEMonlyA_MA, PemCoup_TEEMonlyD_MA, PemH_TEEMonly, EnergyUse_TEEMa] ...
          = PGB_PwrTrq(HybConfig,GBConfig,NemH,NemL,NemCoup,NemH_EPT,NemCoup_EPT, ...
            PXeffH,PXeffL,PXeffH_EPT,PXeffL_EPT,Plps_nom,Phps_nom,Plps_nom_EPT, ...
            Phps_nom_EPT,PEx_ACsys,PX_EPT,Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,PemCoup_maxTEEM,PemCoup_maxNom,...
            TEEMaWeights,TEEMdWeights,PXeffH_bins,transientDur,Psf,Trqsf)

if HybConfig == 1 % standard engine
    % nominal use req.
    PemH_nom_MA = PEx_ACsys;
    PemCoup_nom_MA = 0;
    PemL_nom_MA = 0;
    TrqemH_nom_MA = max(abs(PemH_nom_MA./NemH * 5252.113));
    TrqemCoup_nom_MA = 0;
    TrqemL_nom_MA = 0;
    % EPT req.
    PemCoup_EPT = -PX_EPT./PXeffL_EPT;
    PemH_EPT = (-PXeffL_EPT-PXeffH_EPT).*PemCoup_EPT; %should be about equal to -PemCoup_EPT
    PemCoup_EPT_MA = max(abs(PemCoup_EPT));
    PemH_EPT_MA = max(abs(PemH_EPT));
    PemL_EPT_MA = 0;
    TrqemCoup_EPT = PemCoup_EPT./NemCoup_EPT * 5252.113;
    TrqemH_EPT = PemH_EPT./NemH_EPT * 5252.113;
    TrqemCoup_EPT_MA = max(abs(TrqemCoup_EPT));
    TrqemH_EPT_MA = max(abs(TrqemH_EPT));
    TrqemL_EPT_MA = 0;
    % TEEM req.
    % -- accel
    PXeffH_TEEMa = sum(TEEMaWeights.*PXeffH_bins); %weighted average 
    PHcoup_TEEMa = min([-PemCoup_maxTEEM*PXeffH_TEEMa Pin_TEEMa PExLMax_TEEMa]);
    PemCoup_TEEMa_MA = min([-PHcoup_TEEMa/PXeffH_TEEMa Pin_TEEMa/(1-PXeffH_TEEMa)])/Psf;
    PHcoup_TEEMa = -(PemCoup_TEEMa_MA*Psf)*PXeffH_TEEMa;
    TrqemCoup_TEEMa_MA = max(abs(PemCoup_TEEMa_MA*Psf./NemCoup * 5252.113))/Trqsf;
    PemH_TEEMa_MA = (Pin_TEEMa - PHcoup_TEEMa)/Psf;
    TrqemH_TEEMa_MA = max(abs(PemH_TEEMa_MA*Psf./NemH * 5252.113))/Trqsf;
    PemL_TEEMa_MA = 0;
    TrqemL_TEEMa_MA = 0;
    EnergyUse_TEEMa = (7/10)*745.7*sum([PemH_TEEMa_MA -PemCoup_TEEMa_MA PemL_TEEMa_MA])*transientDur*2.77778e-7*Psf; %kW-hr
        %NOTE: by observation standard engine use ~7/10 less energy
    PemCoup_TEEMonlyA_MA = PemCoup_TEEMa_MA*Psf;
    % -- decel
    PXeffH_TEEMd = sum(TEEMdWeights.*PXeffH_bins); %weighted average
    PHcoup_TEEMd = min(abs([PemCoup_maxTEEM*PXeffH_TEEMd PX_TEEMd]));
    PemCoup_TEEMd_MA = abs(PHcoup_TEEMd/PXeffH_TEEMd)/Psf;
    TrqemCoup_TEEMd_MA = max(abs(PemCoup_TEEMd_MA*Psf./NemCoup * 5252.113))/Trqsf;
    PemH_TEEMd_MA = PemCoup_TEEMd_MA/Psf;
    TrqemH_TEEMd_MA = max(abs(PemH_TEEMd_MA*Psf./NemH * 5252.113))/Trqsf;
    PemL_TEEMd_MA = 0;
    TrqemL_TEEMd_MA = 0;
    PemCoup_TEEMonlyD_MA = PemCoup_TEEMd_MA*Psf;
    PemH_TEEMonly = max([PemH_TEEMa_MA PemCoup_TEEMd_MA])*Psf;
elseif HybConfig == 2 % boost
    % nominal use req.
    PL_nom = Plps_nom;
    PH_nom = Phps_nom;
    for i = 1:length(PXeffL)
        PLcoup_nom(i) = min([PXeffL(i)*PemCoup_maxNom PL_nom(i)]);
        PemCoup_nom(i) = PLcoup_nom(i)/PXeffL(i);
        PemH_nom(i) = PH_nom(i)-PemCoup_nom(i)*PXeffH(i); %leaves out PEx_ACsys as it could reduce
        PemL_nom(i) = PL_nom(i)-PemCoup_nom(i)*PXeffL(i);
    end
    PemCoup_nom_MA = max(abs(PemCoup_nom));
    PemH_nom_MA = max(abs(PemH_nom));
    PemL_nom_MA = max(abs(PemL_nom));
    TrqemCoup_nom = PemCoup_nom./NemCoup * 5252.113;
    TrqemH_nom = PemH_nom./NemH * 5252.113;
    TrqemL_nom = PemL_nom./NemL * 5252.113;
    TrqemCoup_nom_MA = max(abs(TrqemCoup_nom));
    TrqemH_nom_MA = max(abs(TrqemH_nom));
    TrqemL_nom_MA = max(abs(TrqemL_nom));
    % EPT req. (same process as standard engine)
    PemCoup_EPT = -PX_EPT./PXeffL_EPT;
    PemH_EPT = (-PXeffL_EPT-PXeffH_EPT).*PemCoup_EPT; %should be about equal to -PemCoup_EPT
    PemCoup_EPT_MA = max(abs(PemCoup_EPT));
    PemH_EPT_MA = max(abs(PemH_EPT));
    PemL_EPT_MA = 0;
    TrqemCoup_EPT = PemCoup_EPT./NemCoup_EPT * 5252.113;
    TrqemH_EPT = PemH_EPT./NemH_EPT * 5252.113;
    TrqemCoup_EPT_MA = max(abs(TrqemCoup_EPT));
    TrqemH_EPT_MA = max(abs(TrqemH_EPT));
    TrqemL_EPT_MA = 0;
    % TEEM req.
    % -- accel
    PXeffH_TEEMa = sum(TEEMaWeights.*PXeffH_bins); %weighted average 
    PHcoup_TEEMa = min([-PemCoup_maxTEEM*PXeffH_TEEMa Pin_TEEMa PExLMax_TEEMa]);
    PemCoup_TEEMa_TPonly = min([-PHcoup_TEEMa/PXeffH_TEEMa Pin_TEEMa/(1-PXeffH_TEEMa)])/Psf; %transient portion only
    PHcoup_TEEMa = -(PemCoup_TEEMa_TPonly*Psf)*PXeffH_TEEMa;
    PemCoup_TEEMa = PHcoup_TEEMa/PXeffH_TEEMa + PemCoup_nom;
    PemCoup_TEEMa_MA = max(abs(PemCoup_TEEMa))/Psf;
    TrqemCoup_TEEMa_MA = max(abs(PemCoup_TEEMa*Psf./NemCoup * 5252.113))/Trqsf;
    PemH_TEEMa = PemH_nom + Pin_TEEMa - PHcoup_TEEMa;
    PemH_TEEMa_MA = max(abs(PemH_TEEMa))/Psf;
    TrqemH_TEEMa_MA = max(abs(PemH_TEEMa*Psf./NemH * 5252.113))/Trqsf;
    PemL_TEEMa_MA = 0;
    TrqemL_TEEMa_MA = 0;
    EnergyUse_TEEMa = (7/10)*745.7*sum([max(PemH_TEEMa-PemH_nom) -max(abs(PemCoup_TEEMa-PemCoup_nom)) PemL_TEEMa_MA])*transientDur*2.77778e-7*Psf; %kW-hr
        %NOTE: by observation standard engine use ~7/10 less energy
    PemCoup_TEEMonlyA_MA = max(abs(PemCoup_TEEMa-PemCoup_nom))*Psf;
    % -- decel
    PXeffH_TEEMd = sum(TEEMdWeights.*PXeffH_bins); %weighted average
    PHcoup_TEEMd = min(abs([PemCoup_maxTEEM*PXeffH_TEEMd PX_TEEMd]));
    PemCoup_TEEMd = PHcoup_TEEMd/PXeffH_TEEMd + PemCoup_nom;
    PemCoup_TEEMd_MA = max(abs(PemCoup_TEEMd))/Psf;
    TrqemCoup_TEEMd_MA = max(abs(PemCoup_TEEMd_MA*Psf./NemCoup * 5252.113))/Trqsf;
    PemH_TEEMd = PemH_nom - PHcoup_TEEMd/PXeffH_TEEMd;
    PemH_TEEMd_MA = max(PemH_TEEMd)/Psf;
    TrqemH_TEEMd_MA = max(abs(PemH_TEEMd_MA*Psf./NemH * 5252.113))/Trqsf;
    PemL_TEEMd_MA = 0;
    TrqemL_TEEMd_MA = 0;
    PemCoup_TEEMonlyD_MA = abs(PHcoup_TEEMd/PXeffH_TEEMd);
    PemH_TEEMonly = max([max(PemH_TEEMa-PemH_nom)*Psf PemCoup_TEEMonlyD_MA]);
else % Power extraction (PEx)
    % nominal use req.
    PL_nom = Plps_nom;
    PH_nom = Phps_nom;
    for i = 1:length(PXeffL)
        PLcoup_nom(i) = max([-PXeffL(i)*PemCoup_maxNom PL_nom(i)]);
        %PHcoup_nom(i) = max([PXeffH(i)*PemCoup_maxNom PH_nom(i)]);
        PemCoup_nom(i) = PLcoup_nom(i)/PXeffL(i); %max([PLcoup_nom(i)/PXeffL(i) PHcoup_nom(i)/PXeffH(i)]);
        PemH_nom(i) = PH_nom(i)-PemCoup_nom(i)*PXeffH(i)-PEx_ACsys;
        PemL_nom(i) = PL_nom(i)-PemCoup_nom(i)*PXeffL(i);
    end
    PemCoup_nom_MA = max(abs(PemCoup_nom));
    PemH_nom_MA = max(abs(PemH_nom));
    PemL_nom_MA = max(abs(PemL_nom));
    TrqemCoup_nom = PemCoup_nom./NemCoup * 5252.113;
    TrqemH_nom = PemH_nom./NemH * 5252.113;
    TrqemL_nom = PemL_nom./NemL * 5252.113;
    TrqemCoup_nom_MA = max(abs(TrqemCoup_nom));
    TrqemH_nom_MA = max(abs(TrqemH_nom));
    TrqemL_nom_MA = max(abs(TrqemL_nom));
    % EPT req.
    PL_nom_EPT = Plps_nom_EPT;
    PH_nom_EPT = Phps_nom_EPT;
    for i = 1:length(PXeffL_EPT)
        PLcoup_nom_EPT(i) = max([-PXeffL_EPT(i)*PemCoup_maxNom PL_nom_EPT(i)]);
        %PHcoup_nom_EPT(i) = max([PXeffH_EPT(i)*PemCoup_maxNom PH_nom_EPT(i)]);
        PemCoup_nom_EPT(i) = PLcoup_nom_EPT(i)/PXeffL_EPT(i); %max([-PLcoup_nom_EPT(i)/PXeffL_EPT(i) PHcoup_nom_EPT(i)/PXeffH_EPT(i)]);
        PemH_nom_EPT(i) = PH_nom_EPT(i)-PemCoup_nom_EPT(i)*PXeffH_EPT(i)-PEx_ACsys;
        PemL_nom_EPT(i) = PL_nom_EPT(i)-PemCoup_nom_EPT(i)*PXeffL_EPT(i);
    end
    PemCoup_EPT = -PX_EPT./PXeffL_EPT + PemCoup_nom_EPT;
    PemH_EPT = (-PXeffL_EPT-PXeffH_EPT).*PemCoup_EPT + PemH_nom_EPT; %should be about equal to -PemCoup_EPT
    PemCoup_EPT_MA = max(abs(PemCoup_EPT));
    PemH_EPT_MA = max(abs(PemH_EPT));
    PemL_EPT_MA = 0;
    TrqemCoup_EPT = PemCoup_EPT./NemCoup_EPT * 5252.113;
    TrqemH_EPT = PemH_EPT./NemH_EPT * 5252.113;
    TrqemCoup_EPT_MA = max(abs(TrqemCoup_EPT));
    TrqemH_EPT_MA = max(abs(TrqemH_EPT));
    TrqemL_EPT_MA = 0;
    % TEEM req.
    % -- accel
    PXeffH_TEEMa = sum(TEEMaWeights.*PXeffH_bins); %weighted average 
    PHcoup_TEEMa = min([-PemCoup_maxTEEM*PXeffH_TEEMa Pin_TEEMa PExLMax_TEEMa]);
    PemCoup_TEEMa_TPonly = min([-PHcoup_TEEMa/PXeffH_TEEMa Pin_TEEMa/(1-PXeffH_TEEMa)])/Psf; %transient portion only
    PHcoup_TEEMa = -(PemCoup_TEEMa_TPonly*Psf)*PXeffH_TEEMa;
    PemCoup_TEEMa = PHcoup_TEEMa/PXeffH_TEEMa + PemCoup_nom;
    PemCoup_TEEMa_MA = max(abs(PemCoup_TEEMa))/Psf;
    TrqemCoup_TEEMa_MA = max(abs(PemCoup_TEEMa*Psf./NemCoup * 5252.113))/Trqsf;
    PemH_TEEMa = PemH_nom + Pin_TEEMa - PHcoup_TEEMa;
    PemH_TEEMa_MA = max(abs(PemH_TEEMa))/Psf;
    TrqemH_TEEMa_MA = max(abs(PemH_TEEMa*Psf./NemH * 5252.113))/Trqsf;
    PemL_TEEMa_MA = 0;
    TrqemL_TEEMa_MA = 0;
    EnergyUse_TEEMa = (8/9)*745.7*sum([max(PemH_TEEMa-PemH_nom) -max(abs(PemCoup_TEEMa-PemCoup_nom)) PemL_TEEMa_MA])*transientDur*2.77778e-7; %kW-hr
        %NOTE: by observation Boost and PEx use ~8/9 less energy
    PemCoup_TEEMonlyA_MA = max(abs(PemCoup_TEEMa-PemCoup_nom));
    % -- decel
    PXeffH_TEEMd = sum(TEEMdWeights.*PXeffH_bins); %weighted average
    PHcoup_TEEMd = min(abs([PemCoup_maxTEEM*PXeffH_TEEMd PX_TEEMd]));
    PemCoup_TEEMd = PHcoup_TEEMd/PXeffH_TEEMd + PemCoup_nom;
    PemCoup_TEEMd_MA = max(abs(PemCoup_TEEMd))/Psf;
    TrqemCoup_TEEMd_MA = max(abs(PemCoup_TEEMd_MA*Psf./NemCoup * 5252.113))/Trqsf;
    PemH_TEEMd = PemH_nom - PHcoup_TEEMd/PXeffH_TEEMd;
    PemH_TEEMd_MA = max(abs(PemH_TEEMd))/Psf;
    TrqemH_TEEMd_MA = max(abs(PemH_TEEMd_MA*Psf./NemH * 5252.113))/Trqsf;
    PemL_TEEMd_MA = 0;
    TrqemL_TEEMd_MA = 0;
    PemCoup_TEEMonlyD_MA = abs(PHcoup_TEEMd/PXeffH_TEEMd);
    PemH_TEEMonly = max([max(PemH_TEEMa-PemH_nom) PemCoup_TEEMonlyD_MA]);
end

if GBConfig == 1 % HP-Sun, LP-Ring, Coup-Carrier

    PemS_nom_MA = PemH_nom_MA;
    PemS_EPT_MA = PemH_EPT_MA;
    PemS_TEEMa_MA = PemH_TEEMa_MA;
    PemS_TEEMd_MA = PemH_TEEMd_MA;
    TrqemS_nom_MA = TrqemH_nom_MA;
    TrqemS_EPT_MA = TrqemH_EPT_MA;
    TrqemS_TEEMa_MA = TrqemH_TEEMa_MA;
    TrqemS_TEEMd_MA = TrqemH_TEEMd_MA;
    PemR_nom_MA = PemL_nom_MA;
    PemR_EPT_MA = PemL_EPT_MA;
    PemR_TEEMa_MA = PemL_TEEMa_MA;
    PemR_TEEMd_MA = PemL_TEEMd_MA;
    TrqemR_nom_MA = TrqemL_nom_MA;
    TrqemR_EPT_MA = TrqemL_EPT_MA;
    TrqemR_TEEMa_MA = TrqemL_TEEMa_MA;
    TrqemR_TEEMd_MA = TrqemL_TEEMd_MA;
    PemC_nom_MA = PemCoup_nom_MA;
    PemC_EPT_MA = PemCoup_EPT_MA;
    PemC_TEEMa_MA = PemCoup_TEEMa_MA;
    PemC_TEEMd_MA = PemCoup_TEEMd_MA;
    TrqemC_nom_MA = TrqemCoup_nom_MA;
    TrqemC_EPT_MA = TrqemCoup_EPT_MA;
    TrqemC_TEEMa_MA = TrqemCoup_TEEMa_MA;
    TrqemC_TEEMd_MA = TrqemCoup_TEEMd_MA;

elseif GBConfig == 2 % HP-Sun, LP-Carrier, Coup-Ring

    PemS_nom_MA = PemH_nom_MA;
    PemS_EPT_MA = PemH_EPT_MA;
    PemS_TEEMa_MA = PemH_TEEMa_MA;
    PemS_TEEMd_MA = PemH_TEEMd_MA;
    TrqemS_nom_MA = TrqemH_nom_MA;
    TrqemS_EPT_MA = TrqemH_EPT_MA;
    TrqemS_TEEMa_MA = TrqemH_TEEMa_MA;
    TrqemS_TEEMd_MA = TrqemH_TEEMd_MA;
    PemR_nom_MA = PemCoup_nom_MA;
    PemR_EPT_MA = PemCoup_EPT_MA;
    PemR_TEEMa_MA = PemCoup_TEEMa_MA;
    PemR_TEEMd_MA = PemCoup_TEEMd_MA;
    TrqemR_nom_MA = TrqemCoup_nom_MA;
    TrqemR_EPT_MA = TrqemCoup_EPT_MA;
    TrqemR_TEEMa_MA = TrqemCoup_TEEMa_MA;
    TrqemR_TEEMd_MA = TrqemCoup_TEEMd_MA;
    PemC_nom_MA = PemL_nom_MA;
    PemC_EPT_MA = PemL_EPT_MA;
    PemC_TEEMa_MA = PemL_TEEMa_MA;
    PemC_TEEMd_MA = PemL_TEEMd_MA;
    TrqemC_nom_MA = TrqemL_nom_MA;
    TrqemC_EPT_MA = TrqemL_EPT_MA;
    TrqemC_TEEMa_MA = TrqemL_TEEMa_MA;
    TrqemC_TEEMd_MA = TrqemL_TEEMd_MA;

elseif GBConfig == 3 % HP-Ring, LP-Sun, Coup-Carrier

    PemS_nom_MA = PemL_nom_MA;
    PemS_EPT_MA = PemL_EPT_MA;
    PemS_TEEMa_MA = PemL_TEEMa_MA;
    PemS_TEEMd_MA = PemL_TEEMd_MA;
    TrqemS_nom_MA = TrqemL_nom_MA;
    TrqemS_EPT_MA = TrqemL_EPT_MA;
    TrqemS_TEEMa_MA = TrqemL_TEEMa_MA;
    TrqemS_TEEMd_MA = TrqemL_TEEMd_MA;
    PemR_nom_MA = PemH_nom_MA;
    PemR_EPT_MA = PemH_EPT_MA;
    PemR_TEEMa_MA = PemH_TEEMa_MA;
    PemR_TEEMd_MA = PemH_TEEMd_MA;
    TrqemR_nom_MA = TrqemH_nom_MA;
    TrqemR_EPT_MA = TrqemH_EPT_MA;
    TrqemR_TEEMa_MA = TrqemH_TEEMa_MA;
    TrqemR_TEEMd_MA = TrqemH_TEEMd_MA;
    PemC_nom_MA = PemCoup_nom_MA;
    PemC_EPT_MA = PemCoup_EPT_MA;
    PemC_TEEMa_MA = PemCoup_TEEMa_MA;
    PemC_TEEMd_MA = PemCoup_TEEMd_MA;
    TrqemC_nom_MA = TrqemCoup_nom_MA;
    TrqemC_EPT_MA = TrqemCoup_EPT_MA;
    TrqemC_TEEMa_MA = TrqemCoup_TEEMa_MA;
    TrqemC_TEEMd_MA = TrqemCoup_TEEMd_MA;

elseif GBConfig == 4 % HP-Ring, LP-Carrier, Coup-Sun

    PemS_nom_MA = PemCoup_nom_MA;
    PemS_EPT_MA = PemCoup_EPT_MA;
    PemS_TEEMa_MA = PemCoup_TEEMa_MA;
    PemS_TEEMd_MA = PemCoup_TEEMd_MA;
    TrqemS_nom_MA = TrqemCoup_nom_MA;
    TrqemS_EPT_MA = TrqemCoup_EPT_MA;
    TrqemS_TEEMa_MA = TrqemCoup_TEEMa_MA;
    TrqemS_TEEMd_MA = TrqemCoup_TEEMd_MA;
    PemR_nom_MA = PemH_nom_MA;
    PemR_EPT_MA = PemH_EPT_MA;
    PemR_TEEMa_MA = PemH_TEEMa_MA;
    PemR_TEEMd_MA = PemH_TEEMd_MA;
    TrqemR_nom_MA = TrqemH_nom_MA;
    TrqemR_EPT_MA = TrqemH_EPT_MA;
    TrqemR_TEEMa_MA = TrqemH_TEEMa_MA;
    TrqemR_TEEMd_MA = TrqemH_TEEMd_MA;
    PemC_nom_MA = PemL_nom_MA;
    PemC_EPT_MA = PemL_EPT_MA;
    PemC_TEEMa_MA = PemL_TEEMa_MA;
    PemC_TEEMd_MA = PemL_TEEMd_MA;
    TrqemC_nom_MA = TrqemL_nom_MA;
    TrqemC_EPT_MA = TrqemL_EPT_MA;
    TrqemC_TEEMa_MA = TrqemL_TEEMa_MA;
    TrqemC_TEEMd_MA = TrqemL_TEEMd_MA;

elseif GBConfig == 5 % HP-Carrier, LP-Sun, Coup-Ring

    PemS_nom_MA = PemL_nom_MA;
    PemS_EPT_MA = PemL_EPT_MA;
    PemS_TEEMa_MA = PemL_TEEMa_MA;
    PemS_TEEMd_MA = PemL_TEEMd_MA;
    TrqemS_nom_MA = TrqemL_nom_MA;
    TrqemS_EPT_MA = TrqemL_EPT_MA;
    TrqemS_TEEMa_MA = TrqemL_TEEMa_MA;
    TrqemS_TEEMd_MA = TrqemL_TEEMd_MA;
    PemR_nom_MA = PemCoup_nom_MA;
    PemR_EPT_MA = PemCoup_EPT_MA;
    PemR_TEEMa_MA = PemCoup_TEEMa_MA;
    PemR_TEEMd_MA = PemCoup_TEEMd_MA;
    TrqemR_nom_MA = TrqemCoup_nom_MA;
    TrqemR_EPT_MA = TrqemCoup_EPT_MA;
    TrqemR_TEEMa_MA = TrqemCoup_TEEMa_MA;
    TrqemR_TEEMd_MA = TrqemCoup_TEEMd_MA;
    PemC_nom_MA = PemH_nom_MA;
    PemC_EPT_MA = PemH_EPT_MA;
    PemC_TEEMa_MA = PemH_TEEMa_MA;
    PemC_TEEMd_MA = PemH_TEEMd_MA;
    TrqemC_nom_MA = TrqemH_nom_MA;
    TrqemC_EPT_MA = TrqemH_EPT_MA;
    TrqemC_TEEMa_MA = TrqemH_TEEMa_MA;
    TrqemC_TEEMd_MA = TrqemH_TEEMd_MA;

else % HP-Carrier, LP-Ring, Coup-Sun

    PemS_nom_MA = PemCoup_nom_MA;
    PemS_EPT_MA = PemCoup_EPT_MA;
    PemS_TEEMa_MA = PemCoup_TEEMa_MA;
    PemS_TEEMd_MA = PemCoup_TEEMd_MA;
    TrqemS_nom_MA = TrqemCoup_nom_MA;
    TrqemS_EPT_MA = TrqemCoup_EPT_MA;
    TrqemS_TEEMa_MA = TrqemCoup_TEEMa_MA;
    TrqemS_TEEMd_MA = TrqemCoup_TEEMd_MA;
    PemR_nom_MA = PemL_nom_MA;
    PemR_EPT_MA = PemL_EPT_MA;
    PemR_TEEMa_MA = PemL_TEEMa_MA;
    PemR_TEEMd_MA = PemL_TEEMd_MA;
    TrqemR_nom_MA = TrqemL_nom_MA;
    TrqemR_EPT_MA = TrqemL_EPT_MA;
    TrqemR_TEEMa_MA = TrqemL_TEEMa_MA;
    TrqemR_TEEMd_MA = TrqemL_TEEMd_MA;
    PemC_nom_MA = PemH_nom_MA;
    PemC_EPT_MA = PemH_EPT_MA;
    PemC_TEEMa_MA = PemH_TEEMa_MA;
    PemC_TEEMd_MA = PemH_TEEMd_MA;
    TrqemC_nom_MA = TrqemH_nom_MA;
    TrqemC_EPT_MA = TrqemH_EPT_MA;
    TrqemC_TEEMa_MA = TrqemH_TEEMa_MA;
    TrqemC_TEEMd_MA = TrqemH_TEEMd_MA;

end

% PemCoup_TEEMonlyA_MA = 0;
% PemCoup_TEEMonlyD_MA = 0;
% PemH_TEEMonly = 0;

end