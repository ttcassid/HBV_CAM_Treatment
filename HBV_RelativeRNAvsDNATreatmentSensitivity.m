% HBVModelSimulator
close all
clear all
%% Load data
 
ParameterStruc = load('H0731_IndPatientParameters_12Dec24_Rebound'); % load('H0731_IndPatientParameters_13June24'); %
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
PatientNumber = 1;
 
PatientParameterID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37]; 
disp('Convert V_0 into copies/mL from IU/mL')
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
PatientID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37]; 
 
PatientNumber = 1;

PA.dT = 0.004;
PA.lambda = 13000000*PA.dT;

NSteps = 3;
%% Storage for sensitivities
SensitivityDay1 = zeros(2,length(PatientID));
SensitivityDay2 = zeros(2,length(PatientID));
SensitivityDay7 = zeros(2,length(PatientID));
SensitivityDay14 = zeros(2,length(PatientID));
SensitivityDay28 = zeros(2,length(PatientID));
SensitivityComparisonMatrix = zeros(5,length(PatientID));

for ii = 1:  length(PatientID);
    %% Load parameters
    PatientNumber = ii
    X = PatientParameterMatrix(PatientNumber,2:end);
    %Patient specific parameters
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
    PA.VIC = 10^(X(13)).*5.82; %Conversion factor from IU/mL to copies/mL
    PA.AIC = 10^(X(14));
    PA.AN =  X(15);
    PA.cA = X(16); %Fixed
    
    PA.xiP = PA.pi+PA.rhoR+PA.delta+PA.muR;
    PA.xiQ = PA.rhoV+PA.delta+PA.muV;
    
    
    %Increased efficacy  EpsVec = [0.9,0.95,0.99,0.995];
    EpsVec = [PA.epsC ]; %  [0.9,0.95,0.99,0.995]; % [0.9,0.99  ,0.999,0.9999,0.99999];
    % PA.epsC.*ones(1,5); % [0.8,0.95,0.99,0.999,0.9999]; %  [0,0.25,0.9,0.99,0.999]; %[0.98, 0.99, 0.995, 0.999, 0.9995 ]; %  exp( log(X(1)) + [0:NSteps-1].*EpsStep) ;
   
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
    tf =  8; % 28;
    TotalTime = [t0 tf];
    
    %% Solve the ODE systems
    IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
    %% Solve the Adaptive Treated ODE system
%     PA.epsC = EpsVec(1);
    %     PA.mu = MuVec(1);
    [sol] =  HBVMultiscaleModel(IC,TotalTime,PA); 
    
    %% Relative change figures
    
    %% HBV DNA Fig
    % Xi = 2;
    Fig = figure(PatientID(PatientNumber)); % PK data  'Color',[152 0 67]./255  
    g1 =  plot(sol.x(:),log10(sol.y(5,:)./PA.VIC),'LineWidth',1.5,'Color', [118,42,131]/255,'LineStyle','--'); %grey [171,217,233]
    hold on
    %% HBV RNA Fig
    Fig = figure(PatientID(PatientNumber));   % figure(PatientID(end) + PatientID(PatientNumber)) ; %  
    EquilibVec = zeros(1,length(EpsVec));
    TstarVec =  zeros(1,length(EpsVec));
    TstarVecV1 =  zeros(1,length(EpsVec));
    for ii = 1:length(EpsVec)
        Conv = 0.25; % Convex combination coefficient: 0 is closer to 2nd phase, 1 closer to 1st
        EquilibVec(ii) = log( PA.rhoR*PA.alpha.*(1-EpsVec(ii))*IIC./( PA.cR*(PA.rhoR+PA.pi+PA.delta) ) )./log(10);
        % c to xi_p
        %         A =  RIC - PA.rhoR*(PIC/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ) - (PA.alpha*(1-EpsVec(ii))*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ))) - PA.rhoR*(PA.alpha*(1-EpsVec(ii))*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-PA.delta)); ;
        %         B = PA.rhoR*(PIC/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ) - (PA.alpha*(1-EpsVec(ii))*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ))) ;
        %         % xi_p to delta
        %         A = PA.rhoR*(PIC/(PA.cR-(PA.rhoR+PA.pi+PA.delta+PA.muR) ) - (PA.alpha*(1-EpsVec(ii))*IIC/((PA.rhoR+PA.pi+PA.delta+PA.muR) - PA.delta) )*(1/(PA.cR-(PA.rhoR+PA.pi+PA.delta+PA.muR) ))) ;
        %         B = PA.rhoR*(PA.alpha*(1-EpsVec(ii))*IIC/((PA.rhoR+PA.pi+PA.delta+PA.muR) - PA.delta) )*(1/(PA.cR-PA.delta));
        %         % c to delta
        A =  RIC - PA.rhoR*(PIC/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ) - (PA.alpha*(1-EpsVec(ii))*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ))) - PA.rhoR*(PA.alpha*(1-EpsVec(ii))*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-PA.delta)); ;
        B = PA.rhoR*(PA.alpha*(1-EpsVec(ii))*IIC/((PA.rhoR+PA.pi+PA.delta+PA.muR) - PA.delta) )*(1/(PA.cR-PA.delta));
        TstarVecV1(ii) = log( ((1-Conv)/Conv )*(A/B) )./ ( PA.cR-PA.delta ) +1;
        
        %1 plus as the simulation starts at day 1
        % TstarVec(ii) = 1+(log( PA.RIC - PA.rho*PIC*EpsVec(ii)/(PA.c-(PA.rho+PA.pi+PA.delta) ) - (1-EpsVec(ii))*PA.rho*PIC/PA.c ) - log( PA.rho*PIC*EpsVec(ii)/(PA.c-(PA.rho+PA.pi+PA.delta) ) ) )./(PA.c-(PA.rho+PA.pi+PA.delta) );
        %         TstarVecV1(ii) = log( ((1-Conv)/Conv )*(A/B) )./ ( 1*( PA.c-(PA.rho+PA.pi+PA.delta) ) )+1;
        %         TstarVecV1(ii) = log( ((1-Conv)/Conv )*(A/B) )./ ( (PA.rhoR+PA.pi+PA.delta+PA.muR)-PA.delta ) +1;
    end
    %
    %   Fig = figure(PatientID(PatientNumber)+PatientID(end) ); % PK data  'Color',[152 0 67]./25
    TransparencyValue = 0.8*255; 
    g1 =  plot(sol.x(:),log10(sol.y(6,:)./RIC),'LineWidth',1.5,'Color', [118,42,131,TransparencyValue]/255,'LineStyle','-'); %grey % 171,217,233,TransparencyValue
    hold on
    ylim([-3.0 0])
    yticks([-3:0.5:0])
    xticks([0 7])
    xlim([0 inf])
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

