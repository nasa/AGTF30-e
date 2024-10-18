function [PemH_nom_MA, PemH_EPT_MA, PemH_TEEMa_MA, PemH_TEEMd_MA, ...
          TrqemH_nom_MA, TrqemH_EPT_MA, TrqemH_TEEMa_MA, TrqemH_TEEMd_MA, ...
          PemL_nom_MA, PemL_EPT_MA, PemL_TEEMa_MA, PemL_TEEMd_MA, ...
          TrqemL_nom_MA, TrqemL_EPT_MA, TrqemL_TEEMa_MA, TrqemL_TEEMd_MA, ...
          EnergyUse_TEEMa] ...
          = DEM_PwrTrq(HybConfig,NemH,NemL,NemH_EPT,NemL_EPT,...
            Plps_nom,Phps_nom,Plps_nom_EPT,Phps_nom_EPT,PEx_ACsys,PX_EPT, ...
            Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,transientDur,Psf,Trqsf)

if HybConfig == 1 % standard engine
    % nominal use req.
    PemH_nom_MA = PEx_ACsys;
    PemL_nom_MA = 0;
    TrqemH_nom_MA = max(abs(PemH_nom_MA./NemH * 5252.113));
    TrqemL_nom_MA = 0;
    % EPT req.
    PemH_EPT = PX_EPT;
    PemL_EPT = -PX_EPT;
    PemH_EPT_MA = max(abs(PemH_EPT));
    PemL_EPT_MA = max(abs(PemL_EPT));
    TrqemH_EPT = PemH_EPT./NemH_EPT * 5252.113;
    TrqemL_EPT = PemL_EPT./NemL_EPT * 5252.113;
    TrqemH_EPT_MA = max(abs(TrqemH_EPT));
    TrqemL_EPT_MA = max(abs(TrqemL_EPT));
    % TEEM req.
    % -- accel
    PemH_TEEMa_MA = Pin_TEEMa;
    PemL_TEEMa_MA = PExLMax_TEEMa;
    TrqemH_TEEMa_MA = max(abs(PemH_TEEMa_MA*Psf./NemH * 5252.113))/Trqsf;
    TrqemL_TEEMa_MA = max(abs(PemL_TEEMa_MA*Psf./NemL * 5252.113))/Trqsf;
    EnergyUse_TEEMa = 745.7*sum([PemH_TEEMa_MA -PemL_TEEMa_MA])*transientDur*2.77778e-7*Psf; %kW-hr
    % -- decel
    PemH_TEEMd_MA = PX_TEEMd;
    PemL_TEEMd_MA = PX_TEEMd;
    TrqemH_TEEMd_MA = max(abs(PemH_TEEMa_MA*Psf./NemH * 5252.113))/Trqsf;
    TrqemL_TEEMd_MA = max(abs(PemL_TEEMa_MA*Psf./NemL * 5252.113))/Trqsf;
elseif HybConfig == 2 % boost
    % nominal use req.
    PemL_nom = Plps_nom;
    PemH_nom = Phps_nom;
    PemH_nom_MA = max(abs(PemH_nom));
    PemL_nom_MA = max(abs(PemL_nom));
    TrqemH_nom = PemH_nom./NemH * 5252.113;
    TrqemL_nom = PemL_nom./NemL * 5252.113;
    TrqemH_nom_MA = max(abs(TrqemH_nom));
    TrqemL_nom_MA = max(abs(TrqemL_nom));
    % EPT req. (same process as standard engine)
    PemH_EPT = PX_EPT;
    PemL_EPT = -PX_EPT;
    PemH_EPT_MA = max(abs(PemH_EPT));
    PemL_EPT_MA = max(abs(PemL_EPT));
    TrqemH_EPT = PemH_EPT./NemH_EPT * 5252.113;
    TrqemL_EPT = PemL_EPT./NemL_EPT * 5252.113;
    TrqemH_EPT_MA = max(abs(TrqemH_EPT));
    TrqemL_EPT_MA = max(abs(TrqemL_EPT));
    % TEEM req.
    % -- accel
    PemH_TEEMa_MA = max(abs(Pin_TEEMa + PemH_nom));
    PemL_TEEMa_MA = max(abs(-PExLMax_TEEMa + PemL_nom));
    TrqemH_TEEMa_MA = max(abs(PemH_TEEMa_MA*Psf./NemH * 5252.113))/Trqsf;
    TrqemL_TEEMa_MA = max(abs(PemL_TEEMa_MA*Psf./NemL * 5252.113))/Trqsf;
    EnergyUse_TEEMa = 745.7*sum([Pin_TEEMa -PExLMax_TEEMa])*transientDur*2.77778e-7*Psf; %kW-hr
    % -- decel
    PemH_TEEMd_MA = max(abs(PX_TEEMd + PemH_nom));
    PemL_TEEMd_MA = max(abs(PX_TEEMd + PemL_nom));
    TrqemH_TEEMd_MA = max(abs(PemH_TEEMa_MA*Psf./NemH * 5252.113))/Trqsf;
    TrqemL_TEEMd_MA = max(abs(PemL_TEEMa_MA*Psf./NemL * 5252.113))/Trqsf;
else % Power extraction (PEx)
    % nominal use req.
    PemL_nom = Plps_nom;
    PemH_nom = Phps_nom - PEx_ACsys;
    PemH_nom_MA = max(abs(PemH_nom));
    PemL_nom_MA = max(abs(PemL_nom));
    TrqemH_nom = PemH_nom./NemH * 5252.113;
    TrqemL_nom = PemL_nom./NemL * 5252.113;
    TrqemH_nom_MA = max(abs(TrqemH_nom));
    TrqemL_nom_MA = max(abs(TrqemL_nom));
    % EPT req. 
    PemH_EPT = PX_EPT + Phps_nom_EPT;
    PemL_EPT = -PX_EPT + Plps_nom_EPT;
    PemH_EPT_MA = max(abs(PemH_EPT));
    PemL_EPT_MA = max(abs(PemL_EPT));
    TrqemH_EPT = PemH_EPT./NemH_EPT * 5252.113;
    TrqemL_EPT = PemL_EPT./NemL_EPT * 5252.113;
    TrqemH_EPT_MA = max(abs(TrqemH_EPT));
    TrqemL_EPT_MA = max(abs(TrqemL_EPT));
    % TEEM req.
    % -- accel
    PemH_TEEMa_MA = max(abs(Pin_TEEMa + PemH_nom));
    PemL_TEEMa_MA = max(abs(-PExLMax_TEEMa + PemL_nom));
    TrqemH_TEEMa_MA = max(abs(PemH_TEEMa_MA*Psf./NemH * 5252.113))/Trqsf;
    TrqemL_TEEMa_MA = max(abs(PemL_TEEMa_MA*Psf./NemL * 5252.113))/Trqsf;
    EnergyUse_TEEMa = 745.7*sum([Pin_TEEMa -PExLMax_TEEMa])*transientDur*2.77778e-7*Psf; %kW-hr
    % -- decel
    PemH_TEEMd_MA = max(abs(PX_TEEMd + PemH_nom));
    PemL_TEEMd_MA = max(abs(PX_TEEMd + PemL_nom));
    TrqemH_TEEMd_MA = max(abs(PemH_TEEMa_MA*Psf./NemH * 5252.113))/Trqsf;
    TrqemL_TEEMd_MA = max(abs(PemL_TEEMa_MA*Psf./NemL * 5252.113))/Trqsf;
end

end