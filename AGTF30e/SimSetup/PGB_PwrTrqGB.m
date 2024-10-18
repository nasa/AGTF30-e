function [Phps_nom_MA, Phps_EPT_MA, Phps_TEEMa_MA, Phps_TEEMd_MA, ...
          Trqhps_nom_MA, Trqhps_EPT_MA, Trqhps_TEEMa_MA, Trqhps_TEEMd_MA, ...
          Plps_nom_MA, Plps_EPT_MA, Plps_TEEMa_MA, Plps_TEEMd_MA, ...
          Trqlps_nom_MA, Trqlps_EPT_MA, Trqlps_TEEMa_MA, Trqlps_TEEMd_MA] ...
          = PGB_PwrTrqGB(HybConfig,Nhps,Nlps,Nhps_EPT,Nlps_EPT, ...
            Plps_nom,Phps_nom,Plps_nom_EPT,Phps_nom_EPT,PEx_ACsys,PX_EPT, ...
            Pin_TEEMa,PExLMax_TEEMa,PX_TEEMd,Psf,Trqsf)

if HybConfig == 1 % standard engine
    % nominal use req.
    Phps_nom_MA = PEx_ACsys;
    Trqhps_nom_MA = max(abs(Phps_nom_MA./Nhps * 5252.113));
    Plps_nom_MA = 0;
    Trqlps_nom_MA = 0;
    % EPT req.
    Plps_EPT_MA = PX_EPT; 
    Trqlps_EPT_MA = max(abs(Plps_EPT_MA./Nlps_EPT));
    Phps_EPT_MA = PX_EPT; 
    Trqhps_EPT_MA = max(abs(Phps_EPT_MA./Nhps_EPT));
    % TEEM req.
    % -- accel
    Phps_TEEMa_MA = Pin_TEEMa;
    Plps_TEEMa_MA = PExLMax_TEEMa;
    Trqhps_TEEMa_MA = max(abs(Phps_TEEMa_MA*Psf./Nhps * 5252.113))/Trqsf;
    Trqlps_TEEMa_MA = max(abs(Plps_TEEMa_MA*Psf./Nlps * 5252.113))/Trqsf;
    % -- decel
    Plps_TEEMd_MA = PX_TEEMd;
    Phps_TEEMd_MA = PX_TEEMd;
    Trqlps_TEEMd_MA = max(abs(Plps_TEEMd_MA*Psf./Nlps * 5252.113))/Trqsf;
    Trqhps_TEEMd_MA = max(abs(Phps_TEEMd_MA*Psf./Nhps * 5252.113))/Trqsf;
elseif HybConfig == 2 % boost
    % nominal use req.
    Plps_nom_MA = max(abs(Plps_nom));
    Phps_nom_MA = max(abs(Phps_nom));
    Trqlps_nom_MA = max(abs(Plps_nom./Nlps * 5252.113));
    Trqhps_nom_MA = max(abs(Phps_nom./Nhps * 5252.113));
    % EPT req. (same process as standard engine)
    Plps_EPT_MA = PX_EPT; 
    Trqlps_EPT_MA = max(abs(Plps_EPT_MA./Nlps_EPT));
    Phps_EPT_MA = PX_EPT; 
    Trqhps_EPT_MA = max(abs(Phps_EPT_MA./Nhps_EPT));
    % TEEM req.
    % -- accel
    Phps_TEEMa_MA = max(abs(Phps_nom + Pin_TEEMa));
    Trqhps_TEEMa_MA = max(abs((Phps_nom + Pin_TEEMa)*Psf./Nhps * 5252.113))/Trqsf;
    Plps_TEEMa_MA = max(abs(Plps_nom - PExLMax_TEEMa));
    Trqlps_TEEMa_MA = max(abs((Plps_nom - PExLMax_TEEMa)*Psf./Nlps * 5252.113))/Trqsf;
    % -- decel
    Plps_TEEMd_MA = max(abs(Phps_nom + PX_TEEMd));
    Phps_TEEMd_MA = max(abs(Plps_nom - PX_TEEMd));
    Trqhps_TEEMd_MA = max(abs((Phps_nom + PX_TEEMd)*Psf./Nhps * 5252.113))/Trqsf;
    Trqlps_TEEMd_MA = max(abs((Plps_nom - PX_TEEMd)*Psf./Nlps * 5252.113))/Trqsf;
else % Power extraction (PEx)
    % nominal use req.
    Plps_nom_MA = max(abs(Plps_nom));
    Phps_nom_MA = max(abs(Phps_nom));
    Trqlps_nom_MA = max(abs(Plps_nom./Nlps * 5252.113));
    Trqhps_nom_MA = max(abs(Phps_nom./Nhps * 5252.113));
    % EPT req.
    Plps_EPT_MA = max(abs(Plps_nom_EPT - PX_EPT)); 
    Trqlps_EPT_MA = max(abs((Plps_nom_EPT - PX_EPT)./Nlps_EPT));
    Phps_EPT_MA = max(abs(Plps_nom_EPT + PX_EPT)); 
    Trqhps_EPT_MA = max(abs((Phps_nom_EPT + PX_EPT)./Nhps_EPT));
    % TEEM req.
    % -- accel
    Phps_TEEMa_MA = max(abs(Phps_nom + Pin_TEEMa));
    Trqhps_TEEMa_MA = max(abs((Phps_nom + Pin_TEEMa)*Psf./Nhps * 5252.113))/Trqsf;
    Plps_TEEMa_MA = max(abs(Plps_nom - PExLMax_TEEMa));
    Trqlps_TEEMa_MA = max(abs((Plps_nom - PExLMax_TEEMa)*Psf./Nlps * 5252.113))/Trqsf;
    % -- decel
    Plps_TEEMd_MA = max(abs(Phps_nom + PX_TEEMd));
    Phps_TEEMd_MA = max(abs(Plps_nom - PX_TEEMd));
    Trqhps_TEEMd_MA = max(abs((Phps_nom + PX_TEEMd)*Psf./Nhps * 5252.113))/Trqsf;
    Trqlps_TEEMd_MA = max(abs((Plps_nom - PX_TEEMd)*Psf./Nlps * 5252.113))/Trqsf;
end

end