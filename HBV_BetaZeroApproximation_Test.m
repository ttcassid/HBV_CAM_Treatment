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

ParameterStruc = load('H0731_IndPatientParameters_12Dec24_Rebound'); % load('H0731_IndPatientParameters_13June24'); %
disp('Convert V_0 into copies/mL from IU/mL')
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
PatientID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37]; %  [3;5;6;7;8;9;10;11;13;15;16;18;19;20;21;24;26;27;28;29;31;32;33;34;35;36;37];
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
    PA.epsN =  X(12); %0.9962; %  X(8);
    % PA.VIC = 10^(X(13)).*5.82; %Conversion factor from IU/mL to copies/mL
    PA.AIC = 10^(X(14));
    PA.AN =  X(15);
    PA.cA = X(16); %Fixed

    
    
    PA.xiP = PA.pi+PA.rhoR+PA.delta+PA.muR;
    PA.xiQ = PA.rhoV+PA.delta+PA.muV;
    
    %     PA.xiP-PA.xiQ
    
    f = @(x) -(exp(PA.xiP*x(2))-1)./(1+(1-x(1))*(exp(PA.xiP*x(2))-1) ); % d log(P(t)/P_0)/d \eps_C as a function of time
    a = @(x) 1 -( PA.xiQ*PA.alpha/(PA.delta*PA.xiP+PA.pi ) )*(exp(-PA.pi*x(2))-1) ;
    b = @(x) ( PA.pi*PA.xiQ*PA.alpha/(PA.delta*PA.xiP+PA.pi ) )*( (exp( PA.xiQ*x(2) ) -1)/PA.xiQ + (exp(-PA.pi*x(2)) -1)/PA.pi );
    g = @(x) - b(x)./( a(x) + (1-x(1))*b(x) );
    RNASensitivity = f([PA.epsC,14]);
    DNASensitivity = g([PA.epsC,14]);
    %Day 1
    SensitivityDay1(1,ii) = f([PA.epsC,1]);
    SensitivityDay1(2,ii) = g([PA.epsC,1]);
    %Day 2
    SensitivityDay2(1,ii) = f([PA.epsC,2]);
    SensitivityDay2(2,ii) = g([PA.epsC,2]);
    %Day 7
    SensitivityDay7(1,ii) = f([PA.epsC,7]);
    SensitivityDay7(2,ii) = g([PA.epsC,7]);
    %Day 14
    SensitivityDay14(1,ii) = f([PA.epsC,14]);
    SensitivityDay14(2,ii) = g([PA.epsC,14]);
    %Day 28
    SensitivityDay28(1,ii) = f([PA.epsC,28]);
    SensitivityDay28(2,ii) = g([PA.epsC,28]);
    
    
    %Increased efficacy
    %     EpsFinal = 0.99; % 0.9999;
    %     X(1) = PA.epsC ; 0.5;
    %     EpsStep = ( log(EpsFinal)-log(X(1)) )./(NSteps-1);
    EpsVec = [0.9,0.95,0.99,0.995];
    % PA.epsC.*ones(1,5); % [0.8,0.95,0.99,0.999,0.9999]; %  [0,0.25,0.9,0.99,0.999]; %[0.98, 0.99, 0.995, 0.999, 0.9995 ]; %  exp( log(X(1)) + [0:NSteps-1].*EpsStep) ;
    
    % PA.mu = 0;
    % MuVec = [0,0.05,0.1,0.2,0.3,0.4];
    %% Initial condtions
    
    
    % T_e = c_v/(rho_v*beta*M);
    % I_e = lambda/delta - (c_v*dT)/(beta*delta*rho_v*M);
    % P_e = (alpha/(mu_r + delta + pi + rho_r) )*I_e;
    % C_e = lambda*M - (c_v*dT)/(beta*rho_v);
    % R_e = ( (alpha*rho_r)/( c_r*(mu_r + delta + pi + rho_r) ))*I_e;
    % V_e = lambda*rho_v*M/c_v  - dT/beta;
    
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
    tf =  27; % 28;
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
    
    % PA.epsC = EpsVec(5);
    % %     PA.mu = MuVec(5);
    % [sol5] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    %% Solve the approximate Treated ODE system
    PA.beta = 0; 
    PA.epsC = EpsVec(1);
    [solApprox] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PA.epsC = EpsVec(2);
    [sol2Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PA.epsC = EpsVec(3);
    [sol3Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PA.epsC = EpsVec(4);
    [sol4Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    % PA.epsC = EpsVec(5);
    % [sol5Approx] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    
    %% Figures    
    %% HBV approximate solutions Fig 
    
     %% HBV approximate solutions Fig 
    CompareIndex = 6; % 5 for DNA, 6 for RNA 
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

    % Fig = figure(PatientID(PatientNumber)); % PK data  'Color',[152 0 67]./255
    % g1 =  plot(solApprox.x(:),log10(solApprox.y(CompareIndex,:)) - log10(sol.y(CompareIndex,:) ),'LineWidth',2.4,'Color', [171,217,233]/255,'LineStyle','-'); %grey
    % hold on
    % g2 =  plot(sol2Approx.x(:),log10(sol2Approx.y(CompareIndex,:)) - log10(sol2.y(CompareIndex,:) ),'LineWidth',2.4,'Color', [239,138,98]/256,'LineStyle','-'); %grey
    % hold on
    % g3 =  plot(sol3Approx.x(:),log10(sol3Approx.y(CompareIndex,:)) - log10(sol3.y(CompareIndex,:)),'LineWidth',2.4,'Color', [118,42,131]/256,'LineStyle','-'); %grey
    % hold on
    % g4 =  plot(sol4Approx.x(:),log10(sol4Approx.y(CompareIndex,:)) - log10(sol4.y(CompareIndex,:)),'LineWidth',2.4,'Color', [67,147,195]/256,'LineStyle','-'); %grey
    % hold on 
    %% HBV true solutions Fig
%     TransparencyValue = 0.5*255;
%     g1 =  plot(sol.x(:),log10(sol.y(5,:)),'LineWidth',1.75,'Color', [171,217,233,TransparencyValue]/255,'LineStyle','-'); %grey
%     hold on
%     g2 =  plot(sol2.x(:),log10(sol2.y(5,:)),'LineWidth',1.75,'Color', [239,138,98,TransparencyValue]/256,'LineStyle','-'); %grey
%     hold on
%     g3 =  plot(sol3.x(:),log10(sol3.y(5,:)),'LineWidth',1.75,'Color', [118,42,131,TransparencyValue]/256,'LineStyle','-'); %grey
%     hold on
%     g4 =  plot(sol4.x(:),log10(sol4.y(5,:)),'LineWidth',1.75,'Color', [67,147,195,TransparencyValue]/256,'LineStyle','-'); %grey
%     hold on
    ax = gca(Fig); 
    ax.FontSize = 14;
  ylim([-0.3,0.12])
  yticks([-0.3 -0.2 -0.1 0 0.1 ])
  xticks([0 10 20])
  title(PatientID(PatientNumber))
    
    
end


SensitivityComparisonMatrix(1,:) = SensitivityDay1(1,:)./SensitivityDay1(2,:);
SensitivityComparisonMatrix(2,:) = SensitivityDay2(1,:)./SensitivityDay2(2,:);
SensitivityComparisonMatrix(3,:) = SensitivityDay7(1,:)./SensitivityDay7(2,:);
SensitivityComparisonMatrix(4,:) = SensitivityDay14(1,:)./SensitivityDay14(2,:);
SensitivityComparisonMatrix(5,:) = SensitivityDay28(1,:)./SensitivityDay28(2,:);

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

