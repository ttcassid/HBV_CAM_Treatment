% Compare the solution of the full model with the solution of the model
% assuming no new infections after treatment initiation.
close all
clear all
%% Load data 
%% Load parameters

ParameterStruc = load('H0731_IndPatientParameters_13June24');   
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
PatientID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37]; 
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

for ii = 1: length(PatientID);
    %% Load parameters
    PatientNumber = ii;
    X = PatientParameterMatrix(PatientNumber,2:end); 

    % Patient specific parameters
    PA.epsC = X(1); 
    PA.beta = 10^(X(2));
    PA.alpha = X(3);
    PA.delta = X(4);
    PA.pi =  X(5);
    PA.rhoR = X(6);
    PA.rhoV =  X(7);
    PA.cR = X(8); %Fixed
    PA.cV = X(9);
    PA.muR = X(10);
    PA.muV = X(11);
    PA.epsN =  0; 
    PA.AIC = 10^(X(14));
    PA.AN =  X(15);
    PA.cA = X(16); %Fixed 
    
    PA.xiP = PA.pi+PA.rhoR+PA.delta+PA.muR;
    PA.xiQ = PA.rhoV+PA.delta+PA.muV; 
    EpsVec = [0.9,0.95,0.99,0.995];

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
    t0 = 1;
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
    
    PA.epsC = EpsVec(4);
    %     PA.mu = MuVec(4);
    [sol4] =  HBVMultiscaleModel(IC,TotalTime,PA); 
    
    %% Solve the approximate Treated ODE system for no new infections
    PA.beta = 0; 
    PA.epsC = EpsVec(1);
    [solApprox] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PA.epsC = EpsVec(2);
    [sol2Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PA.epsC = EpsVec(3);
    [sol3Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PA.epsC = EpsVec(4);
    [sol4Approx] =  HBVMultiscaleModel(IC,TotalTime,PA); 
     
    
     %% HBV approximate solutions Fig 
    CompareIndex = 5; % 5 for DNA, 6 for RNA 
    TVec = [TotalTime(1):0.05:TotalTime(2)];
    % Approx solution
    Sol1Approx = deval(solApprox,TVec,CompareIndex); 
    Sol2Approx = deval(sol2Approx,TVec,CompareIndex); 
    Sol3Approx = deval(sol3Approx,TVec,CompareIndex);
    Sol4Approx = deval(sol4Approx,TVec,CompareIndex);
    
    % True solution
    Sol1 = deval(sol,TVec,CompareIndex); 
    Sol2 = deval(sol2,TVec,CompareIndex); 
    Sol3 = deval(sol3,TVec,CompareIndex);
    Sol4 = deval(sol4,TVec,CompareIndex);


     Fig = figure(PatientID(PatientNumber)); % PK data  'Color',[152 0 67]./255
    g1 =  plot(TVec,log10(Sol1Approx) - log10(Sol1),'LineWidth',2.4,'Color', [171,217,233]/255,'LineStyle','-'); %grey
    hold on
    g2 =  plot(TVec,log10(Sol2Approx) - log10(Sol2),'LineWidth',2.4,'Color', [239,138,98]/256,'LineStyle','-'); %grey
    hold on
    g3 =  plot(TVec,log10(Sol3Approx) - log10(Sol3),'LineWidth',2.4,'Color', [118,42,131]/256,'LineStyle','-'); %grey
    hold on
    g4 =  plot(TVec,log10(Sol4Approx) - log10(Sol4),'LineWidth',2.4,'Color', [67,147,195]/256,'LineStyle','-'); %grey
    hold on  
    ax = gca(Fig); 
    ax.FontSize = 14;
  ylim([-0.3,0.12])
  yticks([-0.3 -0.2 -0.1 0 0.1 ])
  xticks([0 10 20])
  title(PatientID(PatientNumber))
    
    
end

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

