%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jackson Steiner    %     
% 2024-06-12                 %
% NASA Glenn Research Center %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup script function for AGTF30-e Electrical Power System
% EMTAT v1.3.2 

function [MWS] = Setup_EMTAT_EPS(MWS)

% Energy Storage Device Selection
if MWS.In.Options.HybridConfig == 2
    MWS.EPS.Options.ESD = 0;
else
    MWS.EPS.Options.ESD = 1;    % 0 : battery
                                % 1 : supercapacitor
end

% DC-DC Converter Selection
MWS.EPS.Options.DCDC = 0;   % 0 : regulator
                            % 1 : set point

% Battery Scaling Method (only applicable to battery)
MWS.EPS.Options.BattScale = 0;  % 0 : usable energy
                                % 1 : total energy

%%%% Battery %%%%
Battery_SoC_IC = 85; % Initial State of Charge (%)
V_nom = 2e3; %(0.13457/2); % Nominal Voltage (V)
E_chem = MWS.PowerSystem.ESD.EnergyCapacity * 10^3; % Stored Chemical Energy (Wh)
Q_min = 20; % Minimum usable state of charge (%)
Q_max = 90; % Maximum usable state of charge (%)

% Equivalent Series Resistance (Ohm)
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.Batt.ESR = 0;
else
    MWS.EPS.Batt.ESR = 0.005;
end

% Cell Voltage Charge Curve (1x110 array, V)
Unit_VCC = [3.60022000000000  3.47014090756723  3.40534522081438  3.36797345615108  3.34637300833001  3.33384164719904  3.32652442405484  3.32220404518678  3.31960513667983  3.31799408159663  3.31694886010866  3.31622666940565  3.31568778312557  3.31525175779256  3.31487224411964  3.31452249917241  3.31418705328483  3.31385691696240  3.31352682385014  3.31319364476023  3.31285547523761  3.31251111050485  3.31215974319668  3.31180078921639  3.31143378726315  3.31105834070993  3.31067408381608  3.31028066191027  3.30987771958033  3.30946489343758  3.30904180747880  3.30860806990567  3.30816327074273  3.30770697987113  3.30723874525401  3.30675809122011  3.30626451672401  3.30575749353022  3.30523646428523  3.30470084044935  3.30415000006513  3.30358328534007  3.30300000002157  3.30239940654061  3.30178072289873  3.30114311927020  3.30048571428812  3.29980757097933  3.29910769230854  3.29838501628717  3.29763841059637  3.29686666666691  3.29606849315088  3.29524250871097  3.29438723404271  3.29350108303264  3.29258235294130  3.29162921348325  3.29063969465657  3.28961167315181  3.28854285714290  3.28743076923078  3.28627272727272  3.28506582278477  3.28380689655166  3.28249251101312  3.28111891891878  3.27968202764960  3.27817735849035  3.27659999999974  3.27494455445513  3.27320507614176  3.27137499999957  3.26944705882304  3.26741318681262  3.26526440677901  3.26299069767368  3.26058083832252  3.25802222222128  3.25530063694161  3.25239999999879  3.24930204081497  3.24598591549143  3.24242773722455  3.23859999999805  3.23447086613954  3.23000327868604  3.22515384615104  3.21987142856824  3.21409532709918  3.20775294117234  3.20075670102620  3.19299999999457  3.18435172413165  3.17464878048051  3.16368571427719  3.15119999998996  3.13685074625673  3.12018709675987  3.10059999998259  3.07724615382467  3.04892340422834  3.01385714282250  2.96931351346777  2.91084999993737  2.83073333324329  2.71419999986125  2.52911764682123  2.18979999951267  1.36574285567986]; 
% Cell Voltage Discharge Curve (1x110 array, V)
Unit_VDC = [3.60022000000000  3.47014090756723  3.40534522081438  3.36797345615108  3.34637300833001  3.33384164719904  3.32652442405484  3.32220404518678  3.31960513667983  3.31799408159663  3.31694886010866  3.31622666940565  3.31568778312557  3.31525175779256  3.31487224411964  3.31452249917241  3.31418705328483  3.31385691696240  3.31352682385014  3.31319364476023  3.31285547523761  3.31251111050485  3.31215974319668  3.31180078921639  3.31143378726315  3.31105834070993  3.31067408381608  3.31028066191027  3.30987771958033  3.30946489343758  3.30904180747880  3.30860806990567  3.30816327074273  3.30770697987113  3.30723874525401  3.30675809122011  3.30626451672401  3.30575749353022  3.30523646428523  3.30470084044935  3.30415000006513  3.30358328534007  3.30300000002157  3.30239940654061  3.30178072289873  3.30114311927020  3.30048571428812  3.29980757097933  3.29910769230854  3.29838501628717  3.29763841059637  3.29686666666691  3.29606849315088  3.29524250871097  3.29438723404271  3.29350108303264  3.29258235294130  3.29162921348325  3.29063969465657  3.28961167315181  3.28854285714290  3.28743076923078  3.28627272727272  3.28506582278477  3.28380689655166  3.28249251101312  3.28111891891878  3.27968202764960  3.27817735849035  3.27659999999974  3.27494455445513  3.27320507614176  3.27137499999957  3.26944705882304  3.26741318681262  3.26526440677901  3.26299069767368  3.26058083832252  3.25802222222128  3.25530063694161  3.25239999999879  3.24930204081497  3.24598591549143  3.24242773722455  3.23859999999805  3.23447086613954  3.23000327868604  3.22515384615104  3.21987142856824  3.21409532709918  3.20775294117234  3.20075670102620  3.19299999999457  3.18435172413165  3.17464878048051  3.16368571427719  3.15119999998996  3.13685074625673  3.12018709675987  3.10059999998259  3.07724615382467  3.04892340422834  3.01385714282250  2.96931351346777  2.91084999993737  2.83073333324329  2.71419999986125  2.52911764682123  2.18979999951267  1.36574285567986]; 
% Cell Depth of Discharge (1x110 array, Ah)
Unit_DoD = [0  0.0209908023779841  0.0419816047559682  0.0629724071339524  0.0839632095119365  0.104954011889921  0.125944814267905  0.146935616645889  0.167926419023873  0.188917221401857  0.209908023779841  0.230898826157825  0.251889628535809  0.272880430913794  0.293871233291778  0.314862035669762  0.335852838047746  0.356843640425730  0.377834442803714  0.398825245181698  0.419816047559682  0.440806849937666  0.461797652315651  0.482788454693635  0.503779257071619  0.524770059449603  0.545760861827587  0.566751664205571  0.587742466583555  0.608733268961539  0.629724071339524  0.650714873717508  0.671705676095492  0.692696478473476  0.713687280851460  0.734678083229444  0.755668885607428  0.776659687985412  0.797650490363396  0.818641292741380  0.839632095119365  0.860622897497349  0.881613699875333  0.902604502253317  0.923595304631301  0.944586107009285  0.965576909387269  0.986567711765253  1.00755851414324  1.02854931652122  1.04954011889921  1.07053092127719  1.09152172365517  1.11251252603316  1.13350332841114  1.15449413078913  1.17548493316711  1.19647573554509  1.21746653792308  1.23845734030106  1.25944814267905  1.28043894505703  1.30142974743502  1.32242054981300  1.34341135219098  1.36440215456897  1.38539295694695  1.40638375932494  1.42737456170292  1.44836536408090  1.46935616645889  1.49034696883687  1.51133777121486  1.53232857359284  1.55331937597082  1.57431017834881  1.59530098072679  1.61629178310478  1.63728258548276  1.65827338786075  1.67926419023873  1.70025499261671  1.72124579499470  1.74223659737268  1.76322739975067  1.78421820212865  1.80520900450663  1.82619980688462  1.84719060926260  1.86818141164059  1.88917221401857  1.91016301639655  1.93115381877454  1.95214462115252  1.97313542353051  1.99412622590849  2.01511702828648  2.03610783066446  2.05709863304244  2.07808943542043  2.09908023779841  2.12007104017640  2.14106184255438  2.16205264493236  2.18304344731035  2.20403424968833  2.22502505206632  2.24601585444430  2.26700665682228  2.28799745920027];

% Identify DoD indicies correspondent to minimum and maximum usable states of charge
[~, min_idx] = min(abs(Unit_DoD - ((1 - Q_min / 100) * max(Unit_DoD))));
[~, max_idx] = min(abs(Unit_DoD - ((1 - Q_max / 100) * max(Unit_DoD))));

% Cell nominal voltage approximated as voltage at minimum DoD index
V_cell = Unit_VDC(min_idx);

% Battery Pack Voltage Charge Curve (V) (scaled by V_nom / cell nominal voltage)
MWS.EPS.Batt.VCC = V_nom / V_cell * Unit_VCC;
% Battery Pack Voltage Discharge Curve (V) (scaled by V_nom / cell nominal voltage)
MWS.EPS.Batt.VDC = V_nom / V_cell * Unit_VDC;

% Depth of Discharge Scaling Method
if MWS.EPS.Options.BattScale == 0
    % Usable battery capacity is the length between SoC limits
    usable_capacity = Unit_DoD(min_idx) - Unit_DoD(max_idx);

    % Capacity Table (Depth of Discharge) (Ah) (scaled by E_chem / V_nom / usable_capacity)
    MWS.EPS.Batt.DoD = E_chem / V_nom / usable_capacity * Unit_DoD;
else
    % Total energy scaling method
    MWS.EPS.Batt.DoD = Unit_DoD; % Initalize DoD vector 

    % Integrate VDC over DoD vector to find initial total energy
    E_total = trapz(MWS.EPS.Batt.DoD, MWS.EPS.Batt.VDC);

    a = 1; % DoD scalar

    % Iterate until E_total ~ E_chem
    while (E_total - E_chem)^2 > E_chem/1000000
   
        % Increment or decrement DoD scalar based on sign of error
        if E_total - E_chem > 0
            a = a - 0.001;
        else
            a = a + 0.001;
        end
   
   % Scale DoD vector by adjusted value
   MWS.EPS.Batt.DoD = a * MWS.EPS.Batt.DoD;

   % Integrate to find energy over scaled DoD vector
   E_total = trapz(MWS.EPS.Batt.DoD, MWS.EPS.Batt.VDC);

   end
end

% Depth of Discharge Initial Condition (C)
MWS.EPS.Batt.DoD_IC = 3600 * max(MWS.EPS.Batt.DoD) * (1 - Battery_SoC_IC / 100);

%%%% Supercapacitor %%%%
SC_SoC_IC = 95; % Initial State of Charge (%)

MWS.EPS.SC.Vout_IC = 2e3 + 60; %sqrt((2e3 + 60)^2*85/95); % Output Voltage Initial Condition (V)

% Equivalent Series Resistance (Ohm)
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.SC.ESR = 0;
else
    MWS.EPS.SC.ESR = 0.02;
end

MWS.EPS.SC.TSE = MWS.PowerSystem.ESD.EnergyCapacity * 10^3; % Max Stored Energy (W-h)
MWS.EPS.SC.Iout_IC = 0; % Output Current Initial Condition (A)

% Capacitance (F)
MWS.EPS.SC.C = (2 * SC_SoC_IC * MWS.EPS.SC.TSE) / (100 * MWS.EPS.SC.Vout_IC^2 *0.000277778);

% W(t) [Wh] = V(t)^2 * C / 2 [J] * 0.000277778 [Wh/J] ; SoC = 100 * W(t) / TSE [%]
% Output Voltage Initial Condition (V)
%MWS.EPS.SC.Vout_IC = sqrt(2 * SC_SoC_IC * MWS.EPS.SC.TSE / (100 * MWS.EPS.SC.C * 0.000277778));

%%%% DC-DC Converter (Regulator and Set Point) %%%%
MWS.EPS.DCDC.Vset = 2000; % DC Bus Voltage Set Point (V)
MWS.EPS.DCDC.Vvec = [0 100 1000 1e4]; % Efficiency Table Voltage In Vector
MWS.EPS.DCDC.Ivec = [0 10 100 1e3]; % Efficiency Table Current Out Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.DCDC.Eff = 1*ones(4); % Power Loss Table Values
else
    MWS.EPS.DCDC.Eff = 0.98*ones(4); % Power Loss Table Values
end
MWS.EPS.DCDC.Iout_IC = 0; % Output Current Initial Condition (A)

% DCDC Regulator PID Controller gains NEED TO RETUNE!!!!!!!!
% Optimized for Boost, TEEM On, DEM, Ratio Unit Schedule
if MWS.EPS.Options.ESD == 1     % gains for supercapacitor
    MWS.EPS.DCDC.PID.kp = 0.47;
    MWS.EPS.DCDC.PID.ki = 50; %215;
    MWS.EPS.DCDC.PID.kd = 0;
else                            % gains for battery
    MWS.EPS.DCDC.PID.kp = 0.51;
    MWS.EPS.DCDC.PID.ki = 45; %195;
    MWS.EPS.DCDC.PID.kd = 0;
end

%%%% Cable %%%%
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.Cable.R = 0; % Cable Resistance (Ohm)
else
    MWS.EPS.Cable.R = 0.2; % Cable Resistance (Ohm)
end
MWS.EPS.Cable.Iout_IC = 0; % Output Current Initial Condition (A)

%%%% DC-AC Inverter %%%%
% InvL
MWS.EPS.InvL.Vvec = [0 1e4]; % DC Voltage Vector
MWS.EPS.InvL.Ivec = [0 2e3]; % RMS Line Current Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.InvL.Eff = 1*ones(2); % Power Loss Table Values
else
    MWS.EPS.InvL.Eff = 0.98*ones(2); % Power Loss Table Values
end
MWS.EPS.InvL.K_emf = 0.0503; % Back EMF Constant, V_LL_peak/RPM
% InvH
MWS.EPS.InvH.Vvec = [0 1e4]; % DC Voltage Vector
MWS.EPS.InvH.Ivec = [0 2e3]; % RMS Line Current Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.InvH.Eff = 1*ones(2); % Power Loss Table Values
else
    MWS.EPS.InvH.Eff = 0.98*ones(2); % Power Loss Table Values
end
MWS.EPS.InvH.K_emf = 0.0503; % Back EMF Constant, V_LL_peak/RPM
% InvS
MWS.EPS.InvS.Vvec = [0 1e4]; % DC Voltage Vector
MWS.EPS.InvS.Ivec = [0 2e3]; % RMS Line Current Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.InvS.Eff = 1*ones(2); % Power Loss Table Values
else
    MWS.EPS.InvS.Eff = 0.98*ones(2); % Power Loss Table Values
end
MWS.EPS.InvS.K_emf = 0.0503; % Back EMF Constant, V_LL_peak/RPM
% InvR
MWS.EPS.InvR.Vvec = [0 1e4]; % DC Voltage Vector
MWS.EPS.InvR.Ivec = [0 2e3]; % RMS Line Current Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.InvR.Eff = 1*ones(2); % Power Loss Table Values
else
    MWS.EPS.InvR.Eff = 0.98*ones(2); % Power Loss Table Values
end
MWS.EPS.InvR.K_emf = 0.0503; % Back EMF Constant, V_LL_peak/RPM
% InvC
MWS.EPS.InvC.Vvec = [0 1e4]; % DC Voltage Vector
MWS.EPS.InvC.Ivec = [0 2e3]; % RMS Line Current Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.InvC.Eff = 1*ones(2); % Power Loss Table Values
else
    MWS.EPS.InvC.Eff = 0.98*ones(2); % Power Loss Table Values
end
MWS.EPS.InvC.K_emf = 0.0503; % Back EMF Constant, V_LL_peak/RPM
% InvP
MWS.EPS.InvP.Vvec = [0 1e4]; % DC Voltage Vector
MWS.EPS.InvP.Ivec = [0 2e3]; % RMS Line Current Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.InvP.Eff = 1*ones(2); % Power Loss Table Values
else
    MWS.EPS.InvP.Eff = 0.98*ones(2); % Power Loss Table Values
end
MWS.EPS.InvP.K_emf = 0.0503; % Back EMF Constant, V_LL_peak/RPM

%%%% Electric Machines %%%%
% EML
MWS.EPS.EML.Nvec = [0 2.5e4]; % Efficiency Table Speed Vector
MWS.EPS.EML.Tvec = [-1e3 1e3]; % Efficiency Table Torque Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.EML.Eff =  1*ones(2); % Efficiency Table Values
else
    MWS.EPS.EML.Eff =  0.96*ones(2); %[0.97 0.95; 0.93 0.94]; % Efficiency Table Values
end
% EMH
MWS.EPS.EMH.Nvec = [0 2.5e4]; % Efficiency Table Speed Vector
MWS.EPS.EMH.Tvec = [-1e3 1e3]; % Efficiency Table Torque Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.EMH.Eff =  1*ones(2); % Efficiency Table Values
else
    MWS.EPS.EMH.Eff =  0.96*ones(2); %[0.97 0.95; 0.93 0.94]; % Efficiency Table Values
end
% EMS
MWS.EPS.EMS.Nvec = [0 2.5e4]; % Efficiency Table Speed Vector
MWS.EPS.EMS.Tvec = [-1e3 1e3]; % Efficiency Table Torque Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.EMS.Eff =  1*ones(2); % Efficiency Table Values
else
    MWS.EPS.EMS.Eff =  0.96*ones(2); %[0.97 0.95; 0.93 0.94]; % Efficiency Table Values
end
% EMR
MWS.EPS.EMR.Nvec = [0 2.5e4]; % Efficiency Table Speed Vector
MWS.EPS.EMR.Tvec = [-1e3 1e3]; % Efficiency Table Torque Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.EMR.Eff =  1*ones(2); % Efficiency Table Values
else
    MWS.EPS.EMR.Eff =  0.96*ones(2); %[0.97 0.95; 0.93 0.94]; % Efficiency Table Values
end
% EMC
MWS.EPS.EMC.Nvec = [0 2.5e4]; % Efficiency Table Speed Vector
MWS.EPS.EMC.Tvec = [-1e3 1e3]; % Efficiency Table Torque Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.EMC.Eff =  1*ones(2); % Efficiency Table Values
else
    MWS.EPS.EMC.Eff =  0.96*ones(2); %[0.97 0.95; 0.93 0.94]; % Efficiency Table Values
end
% EMP
MWS.EPS.EMP.Nvec = [0 2.5e4]; % Efficiency Table Speed Vector
MWS.EPS.EMP.Tvec = [-1e3 1e3]; % Efficiency Table Torque Vector
if MWS.In.Options.EPSLosses == 1
    MWS.EPS.EMP.Eff =  1*ones(2); % Efficiency Table Values
else
    MWS.EPS.EMP.Eff =  0.96*ones(2); %[0.97 0.95; 0.93 0.94]; % Efficiency Table Values
end

%%%% Electrical Bus %%%%
MWS.EPS.Bus.C = 6e-3; % Capacitance (F)
MWS.EPS.Bus.Vout_IC = MWS.EPS.DCDC.Vset; % Voltage Initial Condition (V)
MWS.EPS.Bus.R = 0.2; % Resistance (Ohm) (Unused)
MWS.EPS.Bus.L = 6e-3; % Inductance (H) (Unused)

% Initialize SoC sensor condition
if MWS.EPS.Options.ESD == 0
    MWS.Init.IC.SOC = Battery_SoC_IC;
else
    MWS.Init.IC.SOC = SC_SoC_IC;
end

end