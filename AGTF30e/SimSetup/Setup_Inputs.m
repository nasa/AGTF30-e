function [MWS] = Setup_Inputs(MWS,method,In,filename)
%Setup_Inputs.m Summary
%   Written By: Jonathan Kratz
%   Date: June, 2023
%   Description: This function is ran prior to execuation of the AGTF30-e
%   engine model to setup the model inputs. This mainly consists
%   of environmental variables the define the flight profiles as well as
%   control input setting and hybrid configuration. 
%
%   Inputs:
%       method - 1: use Input structure In
%                2: use data in structured excel file given by filename
%                any other entry: use default inputs
%       In - input structure that defines the inputs
%       filename - string with the name (and file extension if necessary) 
%                of the file that contains the data to be used for defining
%                the inputs
%   Outputs:
%       MWS structure with "In" element that defines the inputs

% Simulation Time Step & Control Time Step
MWS.In.Ts = 0.015;

if method == 1 % use input struction In
    
    % Enable Limit Violtation Simulation Halts
    MWS.In.EnLimitStop_Prim = In.EnLimitStop_Prim;
    MWS.In.EnLimitStop_VBV = In.EnLimitStop_VBV;
    MWS.In.EnLimitStop_VAFN = In.EnLimitStop_VAFN;

    % Options
    % -- Hybrid Option (1-Standard Engine, 2-Boost, 3-PEx)
    MWS.In.Options.HybridConfig = In.Options.HybridConfig; 
    % -- Electric Power System Loss Option (1-No losses (idealisitc), 2-Losses (realistic))
    MWS.In.Options.EPSLosses = In.Options.EPSLosses;
    % -- Engine-EM Integration Option (1-Dedicated EM Approach, 2-VEATE Gearbox Approach)
    MWS.In.Options.EngineEMInt = In.Options.EngineEMInt;
    % -- VEATE Gearbox Configuration Option
    % ---- 1 - HP coupled to sun gear, LP coupled to ring gear
    % ---- 2 - HP coupled to sun gear, LP coupled to carrier
    % ---- 3 - HP coupled to ring gear, LP coupled to sun gear
    % ---- 4 - HP coupled to ring gear, LP coupled to carrier
    % ---- 5 - HP coupled to carrier, LP coupled to sun gear
    % ---- 6 - HP coupled to carrier, LP coupled to ring gear
    MWS.In.Options.PGBConfig = In.Options.PGBConfig;
    % -- Fuel Flow Transient Logic Option (1-Generic limiter, 2-Ratio Unit Limiter)
    MWS.In.Options.WfTransientLogic = In.Options.WfTransientLogic;
    % -- Turbine Electrified Energy Managment Option (0-Disable, 1-Enable)
    MWS.In.Options.TEEM = In.Options.TEEM;

    % Engine Input Variables
    % -- Altitude, ft
    MWS.In.t_Alt = In.t_Alt;
    MWS.In.Alt = In.Alt;
    % -- Ambient Temperature Difference from Standard Day, degR or degF
    MWS.In.t_dTamb = In.t_dTamb;
    MWS.In.dTamb = In.dTamb;
    % -- Mach Number
    MWS.In.t_MN = In.t_MN;
    MWS.In.MN = In.MN;
    % -- Aircraft Power Load
    MWS.In.t_PExAC = In.t_PExAC;
    MWS.In.PExAC = In.PExAC;

    % Control Inputs
    % -- Power Lever Angle, deg
    MWS.In.t_PLA = In.t_PLA;
    MWS.In.PLA = In.PLA;
    % -- Boost Toggle
    MWS.In.t_Boost = In.t_Boost;
    MWS.In.Boost = In.Boost;
    % -- EPT Toggle
    MWS.In.t_EPT = In.t_EPT;
    MWS.In.EPT = In.EPT;
    % --- Charging Toggle
    MWS.In.t_Charge = In.t_Charge;
    MWS.In.Charge = In.Charge;
    % -- Manual Control Inputs
    % ---- N1c Set-Point
    MWS.In.N1cManEn = In.N1cManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_N1cMan = In.t_N1cMan;
    MWS.In.N1cMan = In.N1cMan;
    % ---- Fuel Flow, lbm/s
    MWS.In.WfManEn = In.WfManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_WfMan = In.t_WfMan; 
    MWS.In.WfMan = In.WfMan;
    % ---- Variable bleed valve, frac. open
    MWS.In.VBVManEn = In.VBVManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_VBVMan = In.t_VBVMan; 
    MWS.In.VBVMan = In.VBVMan;
    % ---- Variable area fan nozzle, sqin
    MWS.In.VAFNManEn = In.VAFNManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_VAFNMan = In.t_VAFNMan; 
    MWS.In.VAFNMan = In.VAFNMan;
    % ---- LPS Additional Power (nominal), hp
    MWS.In.PwrInLPNomManEn = In.PwrInLPNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPNomMan = In.t_PwrInLPNomMan; 
    MWS.In.PwrInLPNomMan = In.PwrInLPNomMan;
    % ---- LPS Additional Power (off-nominal), hp
    MWS.In.PwrInLPOffNomManEn = In.PwrInLPOffNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPOffNomMan = In.t_PwrInLPOffNomMan; 
    MWS.In.PwrInLPOffNomMan = In.PwrInLPOffNomMan;
    % ---- HPS Additional Power (nominal), hp
    MWS.In.PwrInHPNomManEn = In.PwrInHPNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPNomMan = In.t_PwrInHPNomMan; 
    MWS.In.PwrInHPNomMan = In.PwrInHPNomMan;
    % ---- HPS Additional Power (off-nominal), hp
    MWS.In.PwrInHPOffNomManEn = In.PwrInHPOffNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPOffNomMan = In.t_PwrInHPOffNomMan; 
    MWS.In.PwrInHPOffNomMan = In.PwrInHPOffNomMan;
    % ---- LPS Electric Machine Power (nominal), hp
    MWS.In.PwrInLPEMNomManEn = In.PwrInLPEMNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPEMNomMan = In.t_PwrInLPEMNomMan; 
    MWS.In.PwrInLPEMNomMan = In.PwrInLPEMNomMan;
    % ---- LPS Electric Machine Power (off-nominal), hp
    MWS.In.PwrInLPEMOffNomManEn = In.PwrInLPEMOffNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPEMOffNomMan = In.t_PwrInLPEMOffNomMan; 
    MWS.In.PwrInLPEMOffNomMan = In.PwrInLPEMOffNomMan;
    % ---- HPS Electric Machine Power (nominal), hp
    MWS.In.PwrInHPEMNomManEn = In.PwrInHPEMNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPEMNomMan = In.t_PwrInHPEMNomMan; 
    MWS.In.PwrInHPEMNomMan = In.PwrInHPEMNomMan;
    % ---- HPS Electric Machine Power (off-nominal), hp
    MWS.In.PwrInHPEMOffNomManEn = In.PwrInHPEMOffNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPEMOffNomMan = In.t_PwrInHPEMOffNomMan; 
    MWS.In.PwrInHPEMOffNomMan = In.PwrInHPEMOffNomMan;
    % ---- Sun Gear Electric Machine Power (nominal), hp
    MWS.In.PwrInSEMNomManEn = In.PwrInSEMNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInSEMNomMan = In.t_PwrInSEMNomMan; 
    MWS.In.PwrInSEMNomMan = In.PwrInSEMNomMan;
    % ---- Sun Gear Electric Machine Power (off-nominal), hp
    MWS.In.PwrInSEMOffNomManEn = In.PwrInSEMOffNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInSEMOffNomMan = In.t_PwrInSEMOffNomMan; 
    MWS.In.PwrInSEMOffNomMan = In.PwrInSEMOffNomMan;
    % ---- Ring Gear Electric Machine Power (nominal), hp
    MWS.In.PwrInREMNomManEn = In.PwrInREMNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInREMNomMan = In.t_PwrInREMNomMan; 
    MWS.In.PwrInREMNomMan = In.PwrInREMNomMan;
    % ---- Ring Gear Electric Machine Power (off-nominal), hp
    MWS.In.PwrInREMOffNomManEn = In.PwrInREMOffNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInREMOffNomMan = In.t_PwrInREMOffNomMan; 
    MWS.In.PwrInREMOffNomMan = In.PwrInREMOffNomMan;
    % ---- Carrier Electric Machine Power (nominal), hp
    MWS.In.PwrInCEMNomManEn = In.PwrInCEMNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInCEMNomMan = In.t_PwrInCEMNomMan; 
    MWS.In.PwrInCEMNomMan = In.PwrInCEMNomMan;
    % ---- Carrier Electric Machine Power (off-nominal), hp
    MWS.In.PwrInCEMOffNomManEn = In.PwrInCEMOffNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInCEMOffNomMan = In.t_PwrInCEMOffNomMan; 
    MWS.In.PwrInCEMOffNomMan = In.PwrInCEMOffNomMan;
    % ---- Planet Gear Electric Machine Power (nominal), hp
    MWS.In.PwrInPEMNomManEn = In.PwrInPEMNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInPEMNomMan = In.t_PwrInPEMNomMan; 
    MWS.In.PwrInPEMNomMan = In.PwrInPEMNomMan;
    % ---- Planet Gear Electric Machine Power (off-nominal), hp
    MWS.In.PwrInPEMOffNomManEn = In.PwrInPEMOffNomManEn; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInPEMOffNomMan = In.t_PwrInPEMOffNomMan; 
    MWS.In.PwrInPEMOffNomMan = In.PwrInPEMOffNomMan;
    
elseif method == 2 % use input from formated excel file
    
    %Note that the excel file has to be formated appropriately and the
    %inputs you wish to use should be on the first sheet in the
    %spreadsheet. The file should also be saved in the active directory or 
    %the full file extension should be provided in the file name. The data
    %set length is limited to 10000 time instances per variable. The excel
    %file should be formated as follows (note x is the row number 
    %corresponding to the last entry for each individual data set which 
    %must have time and value vectors of the same length):
    %   Option.HybridConfig: B12
    %   Option.EPSLosses: B13
    %   Option.EngineEMInt: B14
    %   Option.PGBConfig: B15
    %   Option.WfTransientLogic: B16
    %   Option.TEEM: B17
    %   t_Alt: A23:Ax
    %   Alt: B23:Bx
    %   t_dTamb: C23:Cx
    %   dTamb: D23:Dx
    %   t_MN: E23:Ex
    %   MN: F23:Fx
    %   t_PExAC: G23:Gx
    %   PExAC: H23:Hx
    %   t_PLA: J23:Jx
    %   PLA: K23:Kx
    %   t_Boost: L23:Lx
    %   Boost: M23:Mx
    %   t_EPT: N23:Nx
    %   EPT: O23:Ox
    %   t_Charge: P23:Px
    %   Charge: Q23:Qx
    %   N1cManEn = BG21
    %   t_N1cMan: R23:Rx
    %   N1cMan: S23:SSx
    %   WfManEn = BG22
    %   t_WfMan: T23:Tx
    %   WfMan: U23:Ux
    %   VBVManEn = BG23
    %   t_VBVMan: V23:Vx
    %   VBVMan: W23:Wx
    %   VAFNManEn = BG24
    %   t_VAFNMan: X23:Xx
    %   VAFNMan: Y23:Yx
    %   PwrInLPNomManEn = BG25
    %   t_PwrInLPNomMan: Z23:Zx
    %   PwrInLPNomMan: AA23:AAx
    %   PwrInLPOffNomManEn = BG26
    %   t_PwrInLPOffNomMan: AB23:ABx
    %   PwrInLPOffNomMan: AC23:ACx
    %   PwrInHPNomManEn = BG27
    %   t_PwrInHPNomMan: AD23:ADx
    %   PwrInHPNomMan: AE23:AEx
    %   PwrInHPOffNomManEn = BG28
    %   t_PwrInHPOffNomMan: AF23:AFx
    %   PwrInHPOffNomMan: AG23:AGx
    %   PwrInLPEMNomManEn = BG29
    %   t_PwrInLPEMNomMan: AH23:AHx
    %   PwrInLPEMNomMan: AI23:AIx
    %   PwrInLPEMOffNomManEn = BG30
    %   t_PwrInLPEMOffNomMan: AJ23:AJx
    %   PwrInLPEMOffNomMan: AK23:AKx
    %   PwrInHPEMNomManEn = BG31
    %   t_PwrInHPEMNomMan: AL23:ALx
    %   PwrInHPEMNomMan: AM23:AMx
    %   PwrInHPEMOffNomManEn = BG32
    %   t_PwrInHPEMOffNomMan: AN23:ANx
    %   PwrInHPEMOffNomMan: AO23:AOx
    %   PwrInSEMNomManEn = BG33
    %   t_PwrInSEMNomMan: AP23:APx
    %   PwrInSEMNomMan: AQ23:AQx
    %   PwrInSEMOffNomManEn = BG34
    %   t_PwrInSEMOffNomMan: AR23:ARx
    %   PwrInSEMOffNomMan: AS23:ASx
    %   PwrInREMNomManEn = BG35
    %   t_PwrInREMNomMan: AT23:ATx
    %   PwrInREMNomMan: AU23:AUx
    %   PwrInREMOffNomManEn = BG36
    %   t_PwrInREMOffNomMan: AV23:AVx
    %   PwrInREMOffNomMan: AW23:AWx
    %   PwrInCEMNomManEn = BG37
    %   t_PwrInCEMNomMan: AX23:AXx
    %   PwrInCEMNomMan: AY23:AYx
    %   PwrInCEMOffNomManEn = BG38
    %   t_PwrInCEMOffNomMan: AZ23:AZx
    %   PwrInCEMOffNomMan: BA23:BAx
    %   PwrInPEMNomManEn = BG39
    %   t_PwrInPEMNomMan: BB23:BBx
    %   PwrInCEMNomMan: BC23:BCx
    %   PwrInPEMOffNomManEn = BG40
    %   t_PwrInPEMOffNomMan: BD23:BDx
    %   PwrInPEMOffNomMan: BE23:BEx
    %   EnLimitStop_Prim: BJ21:BJ49
    %   EnLimitStop_VBV: BM21:BM22
    %   EnLimitStop_VAFN: BM26:BM27
    
    % Enable Limit Violtation Simulation Halts
    MWS.In.EnLimitStop_Prim = readmatrix(filename,'Range','BJ21:BJ49');
    MWS.In.EnLimitStop_VBV = readmatrix(filename,'Range','BM21:BM22');
    MWS.In.EnLimitStop_VAFN = readmatrix(filename,'Range','BM26:BM27');

    % Options
    % -- Hybrid Option (1-Standard Engine, 2-Boost, 3-PEx)
    MWS.In.Options.HybridConfig = readmatrix(filename,'Range','B12:B12'); 
    % -- Electric Power System Loss Option (1-No losses (idealisitc), 2-Losses (realistic))
    MWS.In.Options.EPSLosses = readmatrix(filename,'Range','B13:B13');
    % -- Engine-EM Integration Option (1-Dedicated EM Approach, 2-VEATE Gearbox Approach)
    MWS.In.Options.EngineEMInt = readmatrix(filename,'Range','B14:B14'); 
    % -- VEATE Gearbox Configuration Option
    % ---- 1 - HP coupled to sun gear, LP coupled to ring gear
    % ---- 2 - HP coupled to sun gear, LP coupled to carrier
    % ---- 3 - HP coupled to ring gear, LP coupled to sun gear
    % ---- 4 - HP coupled to ring gear, LP coupled to carrier
    % ---- 5 - HP coupled to carrier, LP coupled to sun gear
    % ---- 6 - HP coupled to carrier, LP coupled to ring gear
    MWS.In.Options.PGBConfig = readmatrix(filename,'Range','B15:B15'); 
    % -- Fuel Flow Transient Logic Option (1-Generic limiter, 2-Ratio Unit Limiter)
    MWS.In.Options.WfTransientLogic = readmatrix(filename,'Range','B16:B16'); 
    % -- Turbine Electrified Energy Managment Option (0-Disable, 1-Enable)
    MWS.In.Options.TEEM = readmatrix(filename,'Range','B17:B17'); 

    % Engine Input Variables
    % -- Altitude, ft
    MWS.In.t_Alt = readmatrix(filename,'Range','A23:A10023');
    MWS.In.t_Alt = MWS.In.t_Alt(~isnan(MWS.In.t_Alt));
    MWS.In.Alt = readmatrix(filename,'Range','B23:B10023');
    MWS.In.Alt = MWS.In.Alt(~isnan(MWS.In.Alt));
    % -- Ambient Temperature Difference from Standard Day, degR or degF
    MWS.In.t_dTamb = readmatrix(filename,'Range','C23:C10023');
    MWS.In.t_dTamb = MWS.In.t_dTamb(~isnan(MWS.In.t_dTamb));
    MWS.In.dTamb = readmatrix(filename,'Range','D23:D10023');
    MWS.In.dTamb = MWS.In.dTamb(~isnan(MWS.In.dTamb));
    % -- Mach Number
    MWS.In.t_MN = readmatrix(filename,'Range','E23:E10023');
    MWS.In.t_MN = MWS.In.t_MN(~isnan(MWS.In.t_MN));
    MWS.In.MN = readmatrix(filename,'Range','F23:F10023');
    MWS.In.MN = MWS.In.MN(~isnan(MWS.In.MN));
    % -- Aircraft Power Load
    MWS.In.t_PExAC = readmatrix(filename,'Range','G23:G10023');
    MWS.In.t_PExAC = MWS.In.t_PExAC(~isnan(MWS.In.t_PExAC));
    MWS.In.PExAC = readmatrix(filename,'Range','H23:H10023');
    MWS.In.PExAC = MWS.In.PExAC(~isnan(MWS.In.PExAC));

    % Control Inputs
    % -- Power Lever Angle, deg
    MWS.In.t_PLA = readmatrix(filename,'Range','J23:J10023');
    MWS.In.t_PLA = MWS.In.t_PLA(~isnan(MWS.In.t_PLA));
    MWS.In.PLA = readmatrix(filename,'Range','K23:K10023');
    MWS.In.PLA = MWS.In.PLA(~isnan(MWS.In.PLA));
    % -- Boost Toggle
    MWS.In.t_Boost = readmatrix(filename,'Range','L23:L10023');
    MWS.In.t_Boost = MWS.In.t_Boost(~isnan(MWS.In.t_Boost));
    MWS.In.Boost = readmatrix(filename,'Range','M23:M10023');
    MWS.In.Boost = MWS.In.Boost(~isnan(MWS.In.Boost));
    % -- EPT Toggle
    MWS.In.t_EPT = readmatrix(filename,'Range','N23:N10023');
    MWS.In.t_EPT = MWS.In.t_EPT(~isnan(MWS.In.t_EPT));
    MWS.In.EPT = readmatrix(filename,'Range','O23:O10023');
    MWS.In.EPT = MWS.In.EPT(~isnan(MWS.In.EPT));
    % --- Charging Toggle
    MWS.In.t_Charge = readmatrix(filename,'Range','P23:P10023');
    MWS.In.t_Charge = MWS.In.t_Charge(~isnan(MWS.In.t_Charge));
    MWS.In.Charge = readmatrix(filename,'Range','Q23:Q10023');
    MWS.In.Charge = MWS.In.Charge(~isnan(MWS.In.Charge));
    % -- Manual Control Inputs
    % ---- Corrected Fan Speed Command, rpm
    MWS.In.N1cManEn = readmatrix(filename,'Range','BG21:BG21'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_N1cMan = readmatrix(filename,'Range','R23:R10023'); 
    MWS.In.t_N1cMan = MWS.In.t_N1cMan(~isnan(MWS.In.t_N1cMan));
    MWS.In.N1cMan = readmatrix(filename,'Range','S23:S10023');
    MWS.In.N1cMan = MWS.In.N1cMan(~isnan(MWS.In.N1cMan));
    % ---- Fuel Flow, lbm/s
    MWS.In.WfManEn = readmatrix(filename,'Range','BG22:BG22'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_WfMan = readmatrix(filename,'Range','T23:T10023'); 
    MWS.In.t_WfMan = MWS.In.t_WfMan(~isnan(MWS.In.t_WfMan));
    MWS.In.WfMan = readmatrix(filename,'Range','U23:U10023');
    MWS.In.WfMan = MWS.In.WfMan(~isnan(MWS.In.WfMan));
    % ---- Variable bleed valve, frac. open
    MWS.In.VBVManEn = readmatrix(filename,'Range','BG23:BG23'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_VBVMan = readmatrix(filename,'Range','V23:V10023');
    MWS.In.t_VBVMan = MWS.In.t_VBVMan(~isnan(MWS.In.t_VBVMan));
    MWS.In.VBVMan = readmatrix(filename,'Range','W23:W10023');
    MWS.In.VBVMan = MWS.In.VBVMan(~isnan(MWS.In.VBVMan));
    % ---- Variable area fan nozzle, sqin
    MWS.In.VAFNManEn = readmatrix(filename,'Range','BG24:BG24'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_VAFNMan = readmatrix(filename,'Range','X23:X10023');
    MWS.In.t_VAFNMan = MWS.In.t_VAFNMan(~isnan(MWS.In.t_VAFNMan));
    MWS.In.VAFNMan = readmatrix(filename,'Range','Y23:Y10023');
    MWS.In.VAFNMan = MWS.In.VAFNMan(~isnan(MWS.In.VAFNMan));
    % ---- LPS Additional Power (nominal), hp
    MWS.In.PwrInLPNomManEn = readmatrix(filename,'Range','BG25:BG25'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPNomMan = readmatrix(filename,'Range','Z23:Z10023');
    MWS.In.t_PwrInLPNomMan = MWS.In.t_PwrInLPNomMan(~isnan(MWS.In.t_PwrInLPNomMan));
    MWS.In.PwrInLPNomMan = readmatrix(filename,'Range','AA23:AA10023');
    MWS.In.PwrInLPNomMan = MWS.In.PwrInLPNomMan(~isnan(MWS.In.PwrInLPNomMan));
    % ---- LPS Additional Power (off-nominal), hp
    MWS.In.PwrInLPOffNomManEn = readmatrix(filename,'Range','BG26:BG26'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPOffNomMan = readmatrix(filename,'Range','AB23:AB10023');
    MWS.In.t_PwrInLPOffNomMan = MWS.In.t_PwrInLPOffNomMan(~isnan(MWS.In.t_PwrInLPOffNomMan));
    MWS.In.PwrInLPOffNomMan = readmatrix(filename,'Range','AC23:AC10023');
    MWS.In.PwrInLPOffNomMan = MWS.In.PwrInLPOffNomMan(~isnan(MWS.In.PwrInLPOffNomMan));
    % ---- HPS Additional Power (nominal), hp
    MWS.In.PwrInHPNomManEn = readmatrix(filename,'Range','BG27:BG27'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPNomMan = readmatrix(filename,'Range','AD23:AD10023'); 
    MWS.In.t_PwrInHPNomMan = MWS.In.t_PwrInHPNomMan(~isnan(MWS.In.t_PwrInHPNomMan));
    MWS.In.PwrInHPNomMan = readmatrix(filename,'Range','AE23:AE10023');
    MWS.In.PwrInHPNomMan = MWS.In.PwrInHPNomMan(~isnan(MWS.In.PwrInHPNomMan));
    % ---- HPS Additional Power (off-nominal), hp
    MWS.In.PwrInHPOffNomManEn = readmatrix(filename,'Range','BG28:BG28'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPOffNomMan = readmatrix(filename,'Range','AF23:AF10023'); 
    MWS.In.t_PwrInHPOffNomMan = MWS.In.t_PwrInHPOffNomMan(~isnan(MWS.In.t_PwrInHPOffNomMan));
    MWS.In.PwrInHPOffNomMan = readmatrix(filename,'Range','AG23:AG10023');
    MWS.In.PwrInHPOffNomMan = MWS.In.PwrInHPOffNomMan(~isnan(MWS.In.PwrInHPOffNomMan));
    % ---- LPS Electric Machine (nominal), hp
    MWS.In.PwrInLPEMNomManEn = readmatrix(filename,'Range','BG29:BG29'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPEMNomMan = readmatrix(filename,'Range','AH23:AH10023');
    MWS.In.t_PwrInLPEMNomMan = MWS.In.t_PwrInLPEMNomMan(~isnan(MWS.In.t_PwrInLPEMNomMan));
    MWS.In.PwrInLPEMNomMan = readmatrix(filename,'Range','AI23:AI10023');
    MWS.In.PwrInLPEMNomMan = MWS.In.PwrInLPEMNomMan(~isnan(MWS.In.PwrInLPEMNomMan));
    % ---- LPS Electric Machine (off-nominal), hp
    MWS.In.PwrInLPEMOffNomManEn = readmatrix(filename,'Range','BG30:BG30'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPEMOffNomMan = readmatrix(filename,'Range','AJ23:AJ10023');
    MWS.In.t_PwrInLPEMOffNomMan = MWS.In.t_PwrInLPEMOffNomMan(~isnan(MWS.In.t_PwrInLPEMOffNomMan));
    MWS.In.PwrInLPEMOffNomMan = readmatrix(filename,'Range','AK23:AK10023');
    MWS.In.PwrInLPEMOffNomMan = MWS.In.PwrInLPEMOffNomMan(~isnan(MWS.In.PwrInLPEMOffNomMan));
    % ---- HPS Electric Machine (nominal), hp
    MWS.In.PwrInHPEMNomManEn = readmatrix(filename,'Range','BG31:BG31'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPEMNomMan = readmatrix(filename,'Range','AL23:AL10023'); 
    MWS.In.t_PwrInHPEMNomMan = MWS.In.t_PwrInHPEMNomMan(~isnan(MWS.In.t_PwrInHPEMNomMan));
    MWS.In.PwrInHPEMNomMan = readmatrix(filename,'Range','AM23:AM10023');
    MWS.In.PwrInHPEMNomMan = MWS.In.PwrInHPEMNomMan(~isnan(MWS.In.PwrInHPEMNomMan));
    % ---- HPS Electric Machine (off-nominal), hp
    MWS.In.PwrInHPEMOffNomManEn = readmatrix(filename,'Range','BG32:BG32'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPEMOffNomMan = readmatrix(filename,'Range','AN23:AN10023'); 
    MWS.In.t_PwrInHPEMOffNomMan = MWS.In.t_PwrInHPEMOffNomMan(~isnan(MWS.In.t_PwrInHPEMOffNomMan));
    MWS.In.PwrInHPEMOffNomMan = readmatrix(filename,'Range','AO23:AO10023');
    MWS.In.PwrInHPEMOffNomMan = MWS.In.PwrInHPEMOffNomMan(~isnan(MWS.In.PwrInHPEMOffNomMan));
    % ---- Sun Gear Electric Machine (nominal), hp
    MWS.In.PwrInSEMNomManEn = readmatrix(filename,'Range','BG33:BG33'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInSEMNomMan = readmatrix(filename,'Range','AP23:AP10023'); 
    MWS.In.t_PwrInSEMNomMan = MWS.In.t_PwrInSEMNomMan(~isnan(MWS.In.t_PwrInSEMNomMan));
    MWS.In.PwrInSEMNomMan = readmatrix(filename,'Range','AQ23:AQ10023');
    MWS.In.PwrInSEMNomMan = MWS.In.PwrInSEMNomMan(~isnan(MWS.In.PwrInSEMNomMan));
    % ---- Sun Gear Electric Machine (off-nominal), hp
    MWS.In.PwrInSEMOffNomManEn = readmatrix(filename,'Range','BG34:BG34'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInSEMOffNomMan = readmatrix(filename,'Range','AR23:AR10023'); 
    MWS.In.t_PwrInSEMOffNomMan = MWS.In.t_PwrInSEMOffNomMan(~isnan(MWS.In.t_PwrInSEMOffNomMan));
    MWS.In.PwrInSEMOffNomMan = readmatrix(filename,'Range','AS23:AS10023');
    MWS.In.PwrInSEMOffNomMan = MWS.In.PwrInSEMOffNomMan(~isnan(MWS.In.PwrInSEMOffNomMan));
    % ---- Ring Gear Electric Machine (nominal), hp
    MWS.In.PwrInREMNomManEn = readmatrix(filename,'Range','BG35:BG35'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInREMNomMan = readmatrix(filename,'Range','AT23:AT10023'); 
    MWS.In.t_PwrInREMNomMan = MWS.In.t_PwrInREMNomMan(~isnan(MWS.In.t_PwrInREMNomMan));
    MWS.In.PwrInREMNomMan = readmatrix(filename,'Range','AU23:AU10023');
    MWS.In.PwrInREMNomMan = MWS.In.PwrInREMNomMan(~isnan(MWS.In.PwrInREMNomMan));
    % ---- Ring Gear Electric Machine (off-nominal), hp
    MWS.In.PwrInREMOffNomManEn = readmatrix(filename,'Range','BG36:BG36'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInREMOffNomMan = readmatrix(filename,'Range','AV23:AV10023'); 
    MWS.In.t_PwrInREMOffNomMan = MWS.In.t_PwrInREMOffNomMan(~isnan(MWS.In.t_PwrInREMOffNomMan));
    MWS.In.PwrInREMOffNomMan = readmatrix(filename,'Range','AW23:AW10023');
    MWS.In.PwrInREMOffNomMan = MWS.In.PwrInREMOffNomMan(~isnan(MWS.In.PwrInREMOffNomMan));
    % ---- Carrier Electric Machine (nominal), hp
    MWS.In.PwrInCEMNomManEn = readmatrix(filename,'Range','BG37:BG37'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInCEMNomMan = readmatrix(filename,'Range','AX23:AX10023'); 
    MWS.In.t_PwrInCEMNomMan = MWS.In.t_PwrInCEMNomMan(~isnan(MWS.In.t_PwrInCEMNomMan));
    MWS.In.PwrInCEMNomMan = readmatrix(filename,'Range','AY23:AY10023');
    MWS.In.PwrInCEMNomMan = MWS.In.PwrInCEMNomMan(~isnan(MWS.In.PwrInCEMNomMan));
    % ---- Carrier Electric Machine (off-nominal), hp
    MWS.In.PwrInCEMOffNomManEn = readmatrix(filename,'Range','BG38:BG38'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInCEMOffNomMan = readmatrix(filename,'Range','AZ23:AZ10023'); 
    MWS.In.t_PwrInCEMOffNomMan = MWS.In.t_PwrInCEMOffNomMan(~isnan(MWS.In.t_PwrInCEMOffNomMan));
    MWS.In.PwrInCEMOffNomMan = readmatrix(filename,'Range','BA23:BA10023');
    MWS.In.PwrInCEMOffNomMan = MWS.In.PwrInCEMOffNomMan(~isnan(MWS.In.PwrInCEMOffNomMan));
    % ---- Planet Gear Electric Machine (nominal), hp
    MWS.In.PwrInPEMNomManEn = readmatrix(filename,'Range','BG39:BG39'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInPEMNomMan = readmatrix(filename,'Range','BB23:BB10023'); 
    MWS.In.t_PwrInPEMNomMan = MWS.In.t_PwrInPEMNomMan(~isnan(MWS.In.t_PwrInPEMNomMan));
    MWS.In.PwrInPEMNomMan = readmatrix(filename,'Range','BC23:BC10023');
    MWS.In.PwrInPEMNomMan = MWS.In.PwrInPEMNomMan(~isnan(MWS.In.PwrInPEMNomMan));
    % ---- Planet Gear Electric Machine (off-nominal), hp
    MWS.In.PwrInPEMOffNomManEn = readmatrix(filename,'Range','BG40:BG40'); % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInPEMOffNomMan = readmatrix(filename,'Range','BD23:BD10023'); 
    MWS.In.t_PwrInPEMOffNomMan = MWS.In.t_PwrInPEMOffNomMan(~isnan(MWS.In.t_PwrInPEMOffNomMan));
    MWS.In.PwrInPEMOffNomMan = readmatrix(filename,'Range','BE23:BE10023');
    MWS.In.PwrInPEMOffNomMan = MWS.In.PwrInPEMOffNomMan(~isnan(MWS.In.PwrInPEMOffNomMan));

else % use default inputs
    
    % Enable Limit Violtation Simulation Halts
    MWS.In.EnLimitStop_Prim = zeros(29,1);
    MWS.In.EnLimitStop_VBV = [0; 0];
    MWS.In.EnLimitStop_VAFN = [0; 0];

    % Options
    % -- Hybrid Option (1-Standard Engine, 2-Boost, 3-PEx)
    MWS.In.Options.HybridConfig = 1; 
    
    % -- Engine-EM Integration Option (1-Dedicated EM Approach, 2-VEATE Gearbox Approach)
    MWS.In.Options.EngineEMInt = 1;
    % -- VEATE Gearbox Configuration Option
    % ---- 1 - HP coupled to sun gear, LP coupled to ring gear
    % ---- 2 - HP coupled to sun gear, LP coupled to carrier
    % ---- 3 - HP coupled to ring gear, LP coupled to sun gear
    % ---- 4 - HP coupled to ring gear, LP coupled to carrier
    % ---- 5 - HP coupled to carrier, LP coupled to sun gear
    % ---- 6 - HP coupled to carrier, LP coupled to ring gear
    MWS.In.Options.PGBConfig = 2;
    % -- Fuel Flow Transient Logic Option (1-Generic limiter, 2-Ratio Unit Limiter)
    MWS.In.Options.WfTransientLogic = 1;
    % -- Turbine Electrified Energy Managment Option (0-Disable, 1-Enable)
    MWS.In.Options.TEEM = 0;

    % Engine Input Variables
    % -- Altitude, ft
    MWS.In.t_Alt = [0 10];
    MWS.In.Alt = In.Alt;
    % -- Ambient Temperature Difference from Standard Day, degR or degF
    MWS.In.t_dTamb = [0 10];
    MWS.In.dTamb = In.dTamb;
    % -- Mach Number
    MWS.In.t_MN = [0 10];
    MWS.In.MN = In.MN;
    % -- Aircraft Power Load
    MWS.In.t_PExAC = [0 10];
    MWS.In.PExAC = In.PExAC;

    % Control Inputs
    % -- Power Lever Angle, deg
    MWS.In.t_PLA = [0 10];
    MWS.In.PLA = [80 80];
    % -- Boost Toggle
    MWS.In.t_Boost = [0 10];
    MWS.In.Boost = [0 0];
    % -- EPT Toggle
    MWS.In.t_EPT = [0 10];
    MWS.In.EPT = [0 0];
    % --- Charging Toggle
    MWS.In.t_Charge = [0 10];
    MWS.In.Charge = [0 0];
    % -- Manual Control Inputs
    % ---- N1c Set-Point
    MWS.In.N1cManEn = [0 10]; % 0-use control loop, 1-use manual prescription
    MWS.In.t_N1cMan = [0 10];
    MWS.In.N1cMan = [2003.2 2003.2];
    % ---- Fuel Flow, lbm/s
    MWS.In.WfManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_WfMan = [0 10]; 
    MWS.In.WfMan = [1.5379 1.5379];
    % ---- Variable bleed valve, frac. open
    MWS.In.VBVManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_VBVMan = [0 10]; 
    MWS.In.VBVMan = [0 0];
    % ---- Variable area fan nozzle, sqin
    MWS.In.VAFNManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_VAFNMan = [0 10];
    MWS.In.VAFNMan = [6105.5 6105.5];
    % ---- LPS Additional Power (nominal), hp
    MWS.In.PwrInLPNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPNomMan = [0 10];
    MWS.In.PwrInLPNomMan = [0 0];
    % ---- LPS Additional Power (off-nominal), hp
    MWS.In.PwrInLPOffNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPOffNomMan = [0 10]; 
    MWS.In.PwrInLPOffNomMan = [0 0];
    % ---- HPS Additional Power (nominal), hp
    MWS.In.PwrInHPNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPNomMan = [0 10];
    MWS.In.PwrInHPNomMan = [0 0];
    % ---- HPS Additional Power (off-nominal), hp
    MWS.In.PwrInHPOffNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPOffNomMan = [0 10];
    MWS.In.PwrInHPOffNomMan = [0 0];
    % ---- LPS Electric Machine Power (nominal), hp
    MWS.In.PwrInLPEMNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPEMNomMan = [0 10];
    MWS.In.PwrInLPEMNomMan = [0 0];
    % ---- LPS Electric Machine Power (off-nominal), hp
    MWS.In.PwrInLPEMOffNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInLPEMOffNomMan = [0 10];
    MWS.In.PwrInLPEMOffNomMan = [0 0];
    % ---- HPS Electric Machine Power (nominal), hp
    MWS.In.PwrInHPEMNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPEMNomMan = [0 10];
    MWS.In.PwrInHPEMNomMan = [0 0];
    % ---- HPS Electric Machine Power (off-nominal), hp
    MWS.In.PwrInHPEMOffNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInHPEMOffNomMan = [0 10];
    MWS.In.PwrInHPEMOffNomMan = [0 0];
    % ---- Sun Gear Electric Machine Power (nominal), hp
    MWS.In.PwrInSEMNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInSEMNomMan = [0 10];
    MWS.In.PwrInSEMNomMan = [0 0];
    % ---- Sun Gear Electric Machine Power (off-nominal), hp
    MWS.In.PwrInSEMOffNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInSEMOffNomMan = [0 10];
    MWS.In.PwrInSEMOffNomMan = [0 0];
    % ---- Ring Gear Electric Machine Power (nominal), hp
    MWS.In.PwrInREMNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInREMNomMan = [0 10]; 
    MWS.In.PwrInREMNomMan = [0 0];
    % ---- Ring Gear Electric Machine Power (off-nominal), hp
    MWS.In.PwrInREMOffNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInREMOffNomMan = [0 10];
    MWS.In.PwrInREMOffNomMan = [0 0];
    % ---- Carrier Electric Machine Power (nominal), hp
    MWS.In.PwrInCEMNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInCEMNomMan = [0 10]; 
    MWS.In.PwrInCEMNomMan = [0 0];
    % ---- Carrier Electric Machine Power (off-nominal), hp
    MWS.In.PwrInCEMOffNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInCEMOffNomMan = [0 10]; 
    MWS.In.PwrInCEMOffNomMan = [0 0];
    % ---- Planet Gear Electric Machine Power (nominal), hp
    MWS.In.PwrInPEMNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInPEMNomMan = [0 10];
    MWS.In.PwrInPEMNomMan = [0 0];
    % ---- Planet Gear Electric Machine Power (off-nominal), hp
    MWS.In.PwrInPEMOffNomManEn = 0; % 0-use control loop, 1-use manual prescription
    MWS.In.t_PwrInPEMOffNomMan = [0 10];
    MWS.In.PwrInPEMOffNomMan = [0 0];
    
end

% Simulation End Time
MWS.In.tend = max([MWS.In.t_Alt(end) MWS.In.t_dTamb(end) MWS.In.t_MN(end) ...
    MWS.In.t_PExAC(end) MWS.In.t_PLA(end) MWS.In.t_Boost(end) MWS.In.t_EPT(end) ...
    MWS.In.t_Charge(end) MWS.In.t_N1cMan(end) MWS.In.t_WfMan(end) MWS.In.t_VBVMan(end) ...
    MWS.In.t_VAFNMan(end) MWS.In.t_PwrInLPNomMan(end) MWS.In.t_PwrInLPOffNomMan(end) ...
    MWS.In.t_PwrInHPNomMan(end) MWS.In.t_PwrInHPOffNomMan(end) MWS.In.t_PwrInLPEMNomMan(end) ...
    MWS.In.t_PwrInLPEMOffNomMan(end) MWS.In.t_PwrInHPEMNomMan(end) MWS.In.t_PwrInHPEMOffNomMan(end) ...
    MWS.In.t_PwrInSEMNomMan(end) MWS.In.t_PwrInSEMOffNomMan(end) MWS.In.t_PwrInREMNomMan(end) ...
    MWS.In.t_PwrInREMOffNomMan(end) MWS.In.t_PwrInCEMNomMan(end) MWS.In.t_PwrInCEMOffNomMan(end) ...
    MWS.In.t_PwrInPEMNomMan(end) MWS.In.t_PwrInPEMOffNomMan(end)]);

end

