% HBVModelSimulator
close all
clear all
%% Load data
% PatientDataStruc = load('Study2158Data'); % load('VL_CD4_DataFile_Matrix_Censored');
% PatientSpecificDataMatrix = cell2mat(struct2cell(PatientDataStruc));
%
%  PatientID = [1,2,4,6,7,9,11,12,13,19,20,22,23,24,25,26,28,30,31,32];
% PatientDataIndex = 1 + 4*[0:length(PatientID)];
% DataTimeIndex = [1;2;8;14];
% DNAIndex = 4;
% RNAIndex = 5;
%% Load parameters

ParameterStruc = load('H0731_IndPatientParameters_13June24'); %  load('Study101_IndPatientParameters_30Aug2023'); %  load('VL_CD4_DataFile_Matrix_Censored');
disp('Convert V_0 into copies/mL from IU/mL')
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
PatientID = [3;5;6;7;8;9;10;11;13;15;16;18;19;20;21;24;26;27;28;29;31;32;33;34;35;36;37];
% Group100mgID = [];
% Group300mgID = [];
% Group500mgID = [];

PatientNumber = 1;

PA.dT = 0.004;
PA.lambda = 13000000*PA.dT;

NSteps = 5;
%% Storage for sensitivities
SensitivityDay1 = zeros(2,length(PatientID));
SensitivityDay2 = zeros(2,length(PatientID));
SensitivityDay7 = zeros(2,length(PatientID));
SensitivityDay14 = zeros(2,length(PatientID));
SensitivityDay28 = zeros(2,length(PatientID));
SensitivityComparisonMatrix = zeros(5,length(PatientID));

  %% Load parameters 
    %Population estimate  parameters
%% Load parameters
%Population estimate  parameters
HBeAgPos = 1
if HBeAgPos == 1 
    PA.epsC = 0.907; % [0.907,0.97,0.984]
    PA.beta = 10^(-6.375-3.572); % 6.8e-7;
    PA.alpha = 10^(-0.701+3.293) ; %  0.203 
    PA.delta = 0.0245; % 0.062
    PA.pi =  10^(2.311);
    PA.rhoR = 2.478;
    PA.rhoV = 2.478;
    PA.cR = 1; %Fixed
    PA.cV = 1;
    PA.muR = 0;
    PA.muV = 0;
    PA.epsN =  0; %0.9962; %  X(8);
    PA.VIC =  6e8; %Conversion factor from IU/mL to copies/mL
    PA.AIC = 10^(1.584);
    PA.AN =  10^(1.358);
    PA.cA = 0.0565; %Fixed 
elseif HBeAgPos == 0
    PA.epsC = 0.907; % [0.907,0.97,0.984]
    PA.beta = 10^(-6.375 ); % 6.8e-7;
    PA.alpha = 10^(-0.701 ) ; %  0.203 
    PA.delta = 0.0698; % 0.062
    PA.pi =  10^(2.311);
    PA.rhoR = 2.478;
    PA.rhoV = 2.478;
    PA.cR = 1; %Fixed
    PA.cV = 1;
    PA.muR = 0;
    PA.muV = 0;
    PA.epsN =  0; %0.9962; %  X(8);
    PA.VIC =  6e8; %Conversion factor from IU/mL to copies/mL
    PA.AIC = 10^(1.584);
    PA.AN =  10^(1.358);
    PA.cA = 0.0565; %Fixed 
end
    
    %     PA.xiP-PA.xiQ
    
    EpsVec =  [0.919,0.963,0.988]; % [0.911,0.975,0.993]; %  
    %% Initial condtions
    
    M = (PA.pi*PA.alpha)/(PA.delta*(PA.muV + PA.delta + PA.rhoV)*(PA.muR + PA.delta + PA.pi + PA.rhoR) );
    TIC = PA.cV./(PA.rhoV.*PA.beta.*M) ; % PA.lambda/(PA.dT + PA.beta*PA.VIC);
    IIC = PA.lambda/PA.delta - (PA.cV*PA.dT)/(PA.beta*PA.delta*PA.rhoV*M); % PA.beta*TIC*PA.VIC/PA.delta;
    PIC = (PA.alpha/(PA.muR + PA.delta + PA.pi + PA.rhoR) )*IIC; % PA.alpha*IIC/(PA.rhoR+PA.pi+PA.delta+PA.muR);
    RIC = PA.rhoR.*PIC/PA.cR; % (PA.alpha.*PA.rhoR/(PA.cR.*(PA.muR + PA.delta + PA.pi + PA.rhoR)) )*IIC; %
    QIC = PA.lambda.*M-PA.cV.*PA.dT./(PA.beta.*PA.rhoV); % ( PA.pi*PIC )/(PA.rhoV+PA.delta+PA.muV);
    PA.VIC = PA.lambda.*PA.rhoV.*M./PA.cV-PA.dT/PA.beta;
    
    %Fixed parameters
    PA.s = 10^(PA.AN)*PA.cA;
    PA.alphaA = PA.cA*( PA.AIC-PA.s/PA.cA)./IIC;
    %% Simulation time
    t0 = 0;
    tf =  28; % 28;
    TotalTime = [t0 tf];
    
    %% Solve the ODE systems
    IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
    %% Solve the Adaptive Treated ODE system
    PA.epsC = EpsVec(1);
    %     PA.mu = MuVec(1);
    [sol] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PA.epsC = EpsVec(2);
    %     PA.mu = MuVec(2);
    [sol2] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PA.epsC = EpsVec(3);
    %     PA.mu = MuVec(3);
    [sol3] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
%     PA.epsC = EpsVec(4);
%     %     PA.mu = MuVec(4);
%     [sol4] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    % PA.epsC = EpsVec(5);
    % %     PA.mu = MuVec(5);
    % [sol5] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    %% Solve the approximate Treated ODE system
    PA.beta = 0;
    PA.epsC = EpsVec(1);
    [solApprox] =  HBVMultiscaleModel(IC,TotalTime,PA);

     Conv = 0.25; % Convex combination coefficient: 0 is closer to 2nd phase, 1 closer to 1st
    EquilibVec = log( PA.rhoR*PA.alpha.*(1-PA.epsC)*IIC./( PA.cR*(PA.rhoR+PA.pi+PA.delta) ) )./log(10);
    A =  RIC - PA.rhoR*(PIC/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ) - (PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ))) - PA.rhoR*(PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-PA.delta)); ;
    B = PA.rhoR*(PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta+PA.muR) - PA.delta) )*(1/(PA.cR-PA.delta));
    TstarVecV1 = log( ((1-Conv)/Conv )*(A/B) )./ ( PA.cR-PA.delta ) + TotalTime(1) 

    
    PA.epsC = EpsVec(2);
    [sol2Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
     
    EquilibVec = log( PA.rhoR*PA.alpha.*(1-PA.epsC)*IIC./( PA.cR*(PA.rhoR+PA.pi+PA.delta) ) )./log(10);
    A =  RIC - PA.rhoR*(PIC/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ) - (PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ))) - PA.rhoR*(PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-PA.delta)); ;
    B = PA.rhoR*(PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta+PA.muR) - PA.delta) )*(1/(PA.cR-PA.delta));
    TstarVecV1 = log( ((1-Conv)/Conv )*(A/B) )./ ( PA.cR-PA.delta ) + TotalTime(1) 

    PA.epsC = EpsVec(3);
    [sol3Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
     
    EquilibVec = log( PA.rhoR*PA.alpha.*(1-PA.epsC)*IIC./( PA.cR*(PA.rhoR+PA.pi+PA.delta) ) )./log(10);
    A =  RIC - PA.rhoR*(PIC/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ) - (PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ))) - PA.rhoR*(PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-PA.delta)); ;
    B = PA.rhoR*(PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta+PA.muR) - PA.delta) )*(1/(PA.cR-PA.delta));
    TstarVecV1 = log( ((1-Conv)/Conv )*(A/B) )./ ( PA.cR-PA.delta ) + TotalTime(1) 

%     PA.epsC = EpsVec(4);
%     [sol4Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    % PA.epsC = EpsVec(5);
    % [sol5Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    
    %% Figures    
    %% HBV approximate solutions Fig
    % Xi = 2;
    Fig = figure(PatientID(PatientNumber)); % PK data  'Color',[152 0 67]./255
    g3 =  plot(sol3Approx.x(:),log10(sol3Approx.y(5,:)),'LineWidth',1.75,'Color', [67,147,195]/256,'LineStyle','--'); %grey
    hold on
    g2 =  plot(sol2Approx.x(:),log10(sol2Approx.y(5,:)),'LineWidth',1.75,'Color', [239,138,98]/256,'LineStyle','--'); %grey
    hold on
    g1 =  plot(solApprox.x(:),log10(solApprox.y(5,:)),'LineWidth',1.75,'Color', [118,42,13]/255,'LineStyle','--'); %grey
    hold on 
%     g4 =  plot(sol4Approx.x(:),log10(sol4Approx.y(5,:)),'LineWidth',1.75,'Color', [67,147,195]/256,'LineStyle','--'); %grey
%     hold on
    % g5 =  plot(sol5Approx.x(:),log10(sol5Approx.y(5,:)),'LineWidth',1.5,'Color', [179,88,6]/256,'LineStyle','--'); %grey
    % hold on
    %% HBV true solutions Fig
    TransparencyValue = 0.5*255;
    g3 =  plot(sol3.x(:),log10(sol3.y(5,:)),'LineWidth',1.75,'Color', [67,147,195,TransparencyValue]/256,'LineStyle','-'); %grey
    hold on
    g2 =  plot(sol2.x(:),log10(sol2.y(5,:)),'LineWidth',1.75,'Color', [239,138,98,TransparencyValue]/256,'LineStyle','-'); %grey
    hold on
    g1 =  plot(sol.x(:),log10(sol.y(5,:)),'LineWidth',1.75,'Color', [118,42,131,TransparencyValue]/255,'LineStyle','-'); %grey
    hold on 
    %     g4 =  plot(sol4.x(:),log10(sol4.y(5,:)),'LineWidth',1.75,'Color', [67,147,195,TransparencyValue]/256,'LineStyle','-'); %grey
%     hold on
    % g5 =  plot(sol5.x(:),log10(sol5.y(5,:)),'LineWidth',1.5,'Color', [179,88,6,TransparencyValue]/256,'LineStyle','-'); %grey
    % hold on
    %    g12 =  plot([TstarVecV1(1) TstarVecV1(1)],[-6 0] ,'LineWidth',0.75,'Color',[171,217,233]/255,'LineStyle','--'); %grey
    %    hold on
    %    g22 =  plot([TstarVecV1(2) TstarVecV1(2)],[-6 0] ,'LineWidth',0.75,'Color',[239,138,98]/256,'LineStyle','--'); %grey
    %    hold on
    %    g32 =  plot([TstarVecV1(3) TstarVecV1(3)],[-6 0] ,'LineWidth',0.75,'Color',[118,42,131]/256,'LineStyle','--'); %grey
    %    hold on
    %    g42 =  plot([TstarVecV1(4) TstarVecV1(4)],[-6 0] ,'LineWidth',0.75,'Color',[67,147,195]/256,'LineStyle','--'); %grey
    %    hold on
    %    g52 =  plot([TstarVecV1(5) TstarVecV1(5)],[-6 0] ,'LineWidth',0.75,'Color',[179,88,6]/256,'LineStyle','--'); %grey
    %    hold on
    
%     ylabel('log_{10} HBV DNA','FontSize',13); % ,'Interpreter','latex','FontSize',15)
%     xlabel('Time (days)','FontSize',13)
     ylim([-0.0,10])
%     title(PatientID(PatientNumber))




function [sol] = HBVMultiscaleModel(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
sol = ode45(@HBVMultiscaleDynamics,IC,totaltime,opts);
    function dydt = HBVMultiscaleDynamics(t,y);   % [TIC,IIC,PIC,QIC,PA.VIC,PA.RIC,PA.AIC];
        dydt(1) = PA.lambda - PA.beta*y(1)*y(5) - PA.dT*y(1) ; %Differential equation for TIC
        dydt(2) = PA.beta*y(1)*y(5) - PA.delta*y(2);
        dydt(3) = PA.alpha*(1-PA.epsC)*y(2) - (PA.rhoR+PA.pi*(1-PA.epsN)+PA.delta+PA.muR)*y(3);
        dydt(4) = PA.pi*(1-PA.epsN)*y(3) - (PA.rhoV+PA.delta+PA.muV)*y(4); % PA.beta*y(1)*y(5)
        dydt(5) = PA.rhoV*y(4) - PA.cV*y(5);
        dydt(6) = PA.rhoR*y(3) - PA.cR*y(6);
        dydt(7) = PA.s + PA.alphaA*y(2) - PA.cA*y(7);
        dydt = dydt';
    end

end

