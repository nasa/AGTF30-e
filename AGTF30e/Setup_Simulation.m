function MWS = Setup_Simulation(inputMethod,In,filename,Mdl,SSsolverSP)
%Setup_Simulation Summary
%   Written By: Jonathan Kratz
%   Date: 11/20/2019
%   Description: This function is ran prior to execuation of the AGTF30-e
%   engine model to setup the workspace. It calls all setup
%   functions 
%
%   Inputs:
%       inputMethd - input method: 
%           1: use Input structure In
%           2: use data in structured excel file given by filename
%              any other entry: use default inputs
%       In - input structure that defines the inputs
%       filename - string with the name (and file extension if necessary) 
%           of the file that contains the data to be used for defining
%           the inputs. (AGTF30e_Inputs.xlsx provides a template and can be
%           used to defined input profiles)
%       Mdl - Model Selection: 0 - Steady-State, 1 - Dynamic, 2 - Linearization
%       SSsolverSP - steady-state set-point (1 - drive to set N1c, 2 - drive to set T4, 3 - drive to set Fn)
%   Outputs:
%       MWS - Matlab Workspace Structure with all model parameters

% Add Paths --------------------------------------------------------------%

addpath(genpath([cd,'/Functions']));
addpath(genpath([cd,'/SimSetup']));
dir = cd;
cd ..
addpath(genpath([cd,'/TMATS']))
addpath(genpath([cd,'/EMTAT-1.3.2']))
cd(dir);

% Change directory to SimulationSetup folder -----------------------------%

cd([cd,'/SimSetup'])

% Load Engine Bus --------------------------------------------------------%

load('Eng_Bus.mat')
load('SuppSys_Bus.mat')

% Assign MWS, Eng_Bus, and SuppSys Elements to base workspace ------------%

MWS = [];
% MWS
assignin('base','MWS',MWS);
% ---- Amb
assignin('base','Amb',Amb);
% ---- DEM
assignin('base','DEM',DEM);
% ---- Data
assignin('base','Data',Data);
% ---- Eng_Bus
assignin('base','Eng_Bus',Eng_Bus);
% ---- FAN_Data
assignin('base','FAN_Data',FAN_Data);
% ---- GBData
assignin('base','GBData',GBData);
% ---- HPC_Data
assignin('base','HPC_Data',HPC_Data);
% ---- HPT_Data
assignin('base','HPT_Data',HPT_Data);
% ---- LPC_Data
assignin('base','LPC_Data',LPC_Data);
% ---- LPT_Data
assignin('base','LPT_Data',LPT_Data);
% ---- PGB
assignin('base','PGB',PGB);
% ---- PGB_Data
assignin('base','PGB_Data',PGB_Data);
% ---- Perf
assignin('base','Perf',Perf);
% ---- S0
assignin('base','S0',S0);
% ---- S13
assignin('base','S13',S13);
% ---- S17
assignin('base','S17',S17);
% ---- S2
assignin('base','S2',S2);
% ---- S21
assignin('base','S21',S21);
% ---- S22
assignin('base','S22',S22);
% ---- S23
assignin('base','S23',S23);
% ---- S24
assignin('base','S24',S24);
% ---- S25
assignin('base','S25',S25);
% ---- S36
assignin('base','S36',S36);
% ---- S4
assignin('base','S4',S4);
% ---- S45
assignin('base','S45',S45);
% ---- S48
assignin('base','S48',S48);
% ---- S5
assignin('base','S5',S5);
% ---- S7
assignin('base','S7',S7);
% -- Shaft
assignin('base','Shaft',Shaft)
% -- SM
assignin('base','SM',SM)
% -- SuppSys_Bus
assignin('base','SuppSys_Bus',SupplySys)
% % -- Bus_Pwr
% assignin('base','Bus_Pwr',Bus_Pwr)
% % -- DCDC_Pwr
% assignin('base','DCDC_Pwr',DCDC_Pwr)
% % -- SuppCable_Pwr
% assignin('base','SuppCable_Pwr',SuppCable_Pwr)
% % -- ESD_Pwr
% assignin('base','ESD_Pwr',ESD_Pwr)
% % -- V_Bus
% assignin('base','V_Bus',V_Bus)
% % -- V_ESD
% assignin('base','V_ESD',V_ESD)
% % -- I_Bus
% assignin('base','I_Bus',I_Bus)
% % -- I_Load
% assignin('base','I_Load',I_Load)
% % -- I_ESD
% assignin('base','I_ESD',I_ESD)

% Inputs -----------------------------------------------------------------%

MWS = Setup_Inputs(MWS,inputMethod,In,filename);
if inEnvelope(MWS) == 0
    disp('WARNING: Initial flight condition is outside the flight envelope. The model inputs must remain within the operable range.')
end

% Engine -----------------------------------------------------------------%

MWS = setup_AllEng(MWS);	   % develops inputs for engine simulation
                               % includes sensor and actuator setup

% Controller -------------------------------------------------------------%

MWS = setup_Controller(MWS);   % develops inputs for engine controller

% Initialization ---------------------------------------------------------%

startPt = 1; % set to 0 for Steady-State analysis, 1 after
MWS = setup_initialize(MWS,startPt,SSsolverSP);

% Open Model -------------------------------------------------------------%

cd(dir)
if Mdl == 0
    open('AGTF30SysSS_PExPIn.slx');
elseif Mdl == 1
    open('AGTF30SysDyn_PExPIn.slx');
elseif Mdl == 2
    open('AGTF30SysLin_PExPIn.slx');
else
    disp('Invalid entry for Mdl')
end

% Change directory back to model directory -------------------------------%

cd(dir);

% Outputs ----------------------------------------------------------------%

disp('Initialization of the AGTF30-e Engine Model has completed. The model is ready for execution.')

end

