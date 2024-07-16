% Script to plot individual fits for HBV RNA and DNA for ABI-H0731

close all
clear all
%% Load data 
DataTimeIndex = [0;1;7;14;21;28]; 
                
PrintFitFig = 1; % Print the individual fits or not

%% Load parameters from monolix fitting for all participants
PatientHBVDNADataStruc = load('H0731_Study101_HBVDNA_RowDataSet_FullPt');  
PatientSpecificHBVDNAData = cell2mat(struct2cell(PatientHBVDNADataStruc));  

PatientHBVRNADataStruc = load('H0731_Study101_HBVRNA_RowDataSet_FullPt');  
PatientSpecificHBVRNAData = cell2mat(struct2cell(PatientHBVRNADataStruc));  

PatientID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];  

ExcludedPartipants = 0;   
ParameterStruc = load('H0731_IndPatientParameters_13June24'); % load('Study101_IndPatientParameters_FullParticipant_07Feb24') ; % 
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
PatientNumber = 1;

PatientDataID = PatientID; 
PatientParameterID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37]; 
% Find patients that have both data and parameters
[IntersectionDataParameters,ia,ib] = intersect(PatientParameterID,PatientDataID);
PatientID = PatientParameterID(ia);
%% Set up simulation 
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
PatientNumber = 1;

PA.dT = 0.004;
PA.lambda = 13000000*PA.dT;
%% For correlations
    RNAFoldDecay =  zeros(1, length(PatientID));
    DNAFoldDecay =  zeros(1, length(PatientID));
    PredictedFoldDecay = zeros(1, length(PatientID));
    
for ii = 1: length(PatientID);
  
    if ExcludedPartipants ==1;
        %% Load parameters
        PatientNumber = ii;
        X = PatientParameterMatrix(PatientNumber,:);
        %% PatientData
        %Get patient specific data
        PatientDNADataVecPlot = PatientSpecificHBVDNAData(PatientNumber,2:end) ;
        PatientRNADataVecPlot = PatientSpecificHBVRNAData(PatientNumber,2:end) ;
    elseif ExcludedPartipants == 0;
        PatientNumberParameter = ia(ii); % PatientID(ii) ;
        PatientNumberData = ib(ii);
        PatientNumber = PatientParameterMatrix(PatientNumberParameter,1);
        X = PatientParameterMatrix(PatientNumberParameter,2:end);
        PatientDNADataVecPlot = PatientSpecificHBVDNAData(PatientNumberData,2:end) ; % only load treatment data
        PatientRNADataVecPlot = PatientSpecificHBVRNAData(PatientNumberData,2:end) ;
    end
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
    PA.epsN =  0; 
    PA.AIC = 10^(X(14));
    PA.AN =  X(15);
    PA.cA = X(16); %Fixed
    
    
    PA.xiP = PA.pi+PA.rhoR+PA.delta+PA.muR;
    PA.xiQ = PA.rhoV+PA.delta+PA.muV;
    
    
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
    tf =  28;
    TotalTime = [t0 tf];
    
    %% Solve the ODE systems
    IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
    %% Solve the Treated ODE system
    [sol] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    %% HBV DNA Fig
    if PrintFitFig ==1
        Fig = figure(PatientID(ii));
        set(Fig,'defaultAxesColorOrder',[ [118,42,131]./255 ; [33,102,172]./255 ]);
        yyaxis left
        % Shading for treatment duration
        Ylowlim = -1.0; %0.0
        YUpLim = 10; 
        hold on 
        g10 =  plot(sol.x(:),log(sol.y(5,:)./5.82)./log(10),'LineWidth',2,'Color', [153,112,171]/255,'LineStyle','-'); %grey 153,112,171
        hold on
        g60 = scatter(DataTimeIndex,PatientDNADataVecPlot,40,[118,42,131]/255,'o','filled'); %,'20', [33,102,172]/255,'*');
        hold on 
        ylim([Ylowlim,YUpLim])
        xlim([-0.5,inf])
        yticks([0:2:10]);
        ylabel('log_{10} HBV DNA (copies/mL) ','FontSize',15);
        %% HBV RNA Fig
        yyaxis right
        g1 =  plot(sol.x(:),log(sol.y(6,:))./log(10),'LineWidth',2,'Color', [67,147,195]/255,'LineStyle','-'); %grey
        hold on
        g6 = scatter(DataTimeIndex,PatientRNADataVecPlot,40,[33,102,172]/255,'o','filled'); %,'20', [33,102,172]/255,'*');
        hold on  
        ylim([Ylowlim,YUpLim])
        xlim([-0.5,inf])
        yticks([0:2:10]);
        ylabel('log_{10} HBV RNA (copies/mL) ','FontSize',15); % ,'Interpreter','latex','FontSize',15)
        xlabel('Time (days)','FontSize',12)
        title(PatientID(ii))

    end
    
 
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

