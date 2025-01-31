% HBVModelSimulator
close all
clear all

% PatientDataIndex = 1 + 7*[0:length(PatientID)];
%% Random sampling for plotting
% PtID0731 = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37]; 
% DoseVector = [100;200;100;100;200;200;200;200;300;300;100;100;100;200;100;100;100;200;200;300;300;300;300;300;300;200;200;300;100]; 
% Cohort100mg = PtID0731(DoseVector == 100);
% Cohort200mg = PtID0731(DoseVector == 200);
% Cohort300mg = PtID0731(DoseVector == 300);
% IDtoPlot0731 = [randsample(Cohort100mg,1);randsample(Cohort200mg,1);randsample(Cohort300mg,2)]
% 
% PtID2158 = [1;2;3;4;6;7;9;11;12;13;15;16;18;19;20;22;23;24;25;26;28;30;31;32;35;37] ;
% DoseVector =  [100;100;600;100;300;500;100;300;300;500;600;600;600;100;300;300;500;100;100;500;300;500;500;500;600;600]; 
% Cohort100mg = PtID2158(DoseVector == 100);
% Cohort300mg = PtID2158(DoseVector == 300);
% Cohort500mg = PtID2158(DoseVector == 500);
% Cohort600mg = PtID2158(DoseVector == 600);
% IDtoPlot2158 = [randsample(Cohort100mg,1);randsample(Cohort300mg,2);randsample(Cohort500mg,1);randsample(Cohort600mg,1)]

 %% print fig or not.               
PrintFitFig = 1; % Print the individual fits or not

%% Choose which trial to simulate
TrialToSimulate = 0731; %  0731; % 

if TrialToSimulate == 0731
    %% Load data
    PatientHBVDNADataStruc = load('H0731_Study101_HBVDNA_RowDataSet_FullPt_Rebound');
    PatientSpecificHBVDNAData = cell2mat(struct2cell(PatientHBVDNADataStruc));
    
    PatientHBVRNADataStruc = load('H0731_Study101_HBVRNA_RowDataSet_FullPt_Rebound');
    PatientSpecificHBVRNAData = cell2mat(struct2cell(PatientHBVRNADataStruc));
    
        PatientHBVALTDataStruc = load('H0731_Study101_ALT_RowDataSet_FullPt_Rebound');
    PatientSpecificHBVALTData = cell2mat(struct2cell(PatientHBVALTDataStruc));
    
    PatientID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];  % [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];
                
    %% Load parameters from monolix fitting for all participants
    ExcludedPartipants = 0;
    ParameterStruc = load('H0731_IndPatientParameters_12Dec24_Rebound'); % load('H0731_IndPatientParameters_13June24'); %
    PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
    PatientNumber = 1;
    
    PatientDataID = PatientID;
    PatientParameterID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];
    % Find patients that have both data and parameters
    [IntersectionDataParameters,ia,ib] = intersect(PatientParameterID,PatientDataID);
    PatientID = PatientParameterID(ia); 
    PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
    PatientNumber = 1;
    DoseVector = [100;200;100;100;200;200;200;200;300;300;100;100;100;200;100;100;100;200;200;300;300;300;300;300;300;200;200;300;100]; 
    DoseTextLocation = 9; % PA.TreatmentEnd/3;
    % Trial specific
    DataTimeIndex = [0;1;7;14;21;28;35;42;56];  
 
    PA.TreatmentEnd = 28;
    DNALLoQ = 0.95; % 20 IU/mL
    RNALLoQ =  2.49; %% RNA LLoQ 250 copies/mL
elseif TrialToSimulate == 2158
    %% Load data
    PatientHBVDNADataStruc = load('H2158_Study101_HBVDNA_RowDataSet_Rebound');
    PatientSpecificHBVDNAData = cell2mat(struct2cell(PatientHBVDNADataStruc));
    PatientHBVRNADataStruc = load('H2158_Study101_HBVRNA_RowDataSet_Rebound');
    PatientSpecificHBVRNAData = cell2mat(struct2cell(PatientHBVRNADataStruc));
    
    PatientHBVALTDataStruc = load('H2158_Study101_ALT_RowDataSet_Rebound');
    PatientSpecificHBVALTData = cell2mat(struct2cell(PatientHBVALTDataStruc));
    
    PatientDataID = [1;2;3;4;6;7;9;11;12;13;15;16;18;19;20;22;23;24;25;26;27;28;30;31;32;33;34;35;37] ;
    DataTimeIndex = [0;1;7;14;15;21;28;56]; % 15;22;29;57]; % [0;1;2;8;15;22;29];
    
    %% Load parameters from monolix fitting % H2158_Study101_IndPatientParameters_SingleRebound_27Oct2023V1
    ParameterStruc = load('H2158_IndPatientParameters_Run206_18June24'); 
    % load('H2158_IndPatientParameters_18June24');  % Run 100
    PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
    PatientParameterID =  [1;2;3;4;6;7;9;11;12;13;15;16;18;19;20;22;23;24;25;26;28;30;31;32;35;37] ;
    % Find patients that have both data and parameters
    [IntersectionDataParameters,ia,ib] = intersect(PatientParameterID,PatientDataID);
    PatientID = PatientParameterID(ia);
    DoseVector = [100;100;600;100;300;500;100;300;300;500;600;600;600;100;300;300;500;100;100;500;300;500;500;500;600;600]; 
    DoseTextLocation = 0.85;
    % Trial specific
    DataTimeIndex = [0;1;7;14;15;21;28;56];
    PA.TreatmentEnd = 14;
    DNALLoQ = 1.30; % 20 IU/mL
    RNALLoQ =  2.39; %% RNA LLoQ 250 copies/mL
end

%% Set up simulation 
% fixed parameters
PA.dT = 0.004;
PA.lambda = 13000000*PA.dT;
%% Simulation time
t0 = 0;
tf =  59;
TotalTime = [t0 tf];

TStarVec = zeros(1,length(PatientID) );
%% Simulate model
for ii =  1:   length(PatientID);
    PatientNumberParameter = ia(ii); % PatientID(ii) ;
    PatientNumberData = ib(ii);
    PatientNumber = PatientParameterMatrix(PatientNumberParameter,1);
    X = PatientParameterMatrix(PatientNumberParameter,2:end);
    PatientDNADataVecPlot = PatientSpecificHBVDNAData(PatientNumberData,2:end) ; % only load treatment data
    PatientRNADataVecPlot = PatientSpecificHBVRNAData(PatientNumberData,2:end) ;  
    PatientALTDataVecPlot = 10.^(PatientSpecificHBVALTData(PatientNumberData,2:end)) ;
    
    %Patient specific parameters 
    PA.epsC = X(1);
    PA.beta = 10^(X(2));
    PA.alpha = X(3);
    PA.delta = X(4);
    PA.pi =  X(5); 
    PA.rhoR = X(6); % Repeated column in parameter set
    PA.rhoV =  X(7); % Repeated column in parameter set
    PA.cR = X(8); %Fixed
    PA.cV = X(9);
    PA.muR = X(10);
    PA.muV = X(11);
    PA.epsN =  X(12); %0.9962; %  X(8);
    % PA.VIC = 10^(X(13)).*5.82; %Conversion factor from IU/mL to copies/mL
    PA.AIC = 10^(X(14));
    PA.AN =  X(15);
    PA.cA = X(16); %Fixed
    PA.k = X(17);

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
    %% Solve the ODE systems
    IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
    %% Solve the Treated ODE system
    [sol] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
      %% Calculate t^*
    Conv = 0.2; % Convex combination coefficient: 0 is closer to 2nd phase, 1 closer to 1st
    EquilibVec = log( PA.rhoR*PA.alpha.*(1-PA.epsC)*IIC./( PA.cR*(PA.rhoR+PA.pi+PA.delta) ) )./log(10);
    A =  RIC - PA.rhoR*(PIC/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ) - (PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-(PA.rhoR+PA.pi+PA.delta) ))) - PA.rhoR*(PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta) - PA.delta) )*(1/(PA.cR-PA.delta)); ;
    B = PA.rhoR*(PA.alpha*(1-PA.epsC)*IIC/((PA.rhoR+PA.pi+PA.delta+PA.muR) - PA.delta) )*(1/(PA.cR-PA.delta));
    TstarVecV1 = log( ((1-Conv)/Conv )*(A/B) )./ ( PA.cR-PA.delta ) + TotalTime(1)  ; 
    TStarVec(ii) =  TstarVecV1; 
    
    %% HBV DNA Fig
    if PrintFitFig ==1
        DNAColor = [221,17,17]./255;
        RNAColor = [0,0,139]./255;
        ALTColor = [10,10,10]./255;  
        Fig = figure(PatientID(ii));
        % Shading for treatment duration
        if TrialToSimulate == 0731;
        Ylowlim = 0; %0.0
        YUpLim = 425;
        elseif TrialToSimulate == 2158;
        Ylowlim = 0; %0.0
        YUpLim = 190; 
        end
                
        TreatmentShadeX  = [0 PA.TreatmentEnd+0.25 PA.TreatmentEnd+0.25  0];
        TreatmentShadeY = [Ylowlim Ylowlim YUpLim YUpLim];
        fill(TreatmentShadeX,TreatmentShadeY,[186,186,186]./255,'FaceAlpha',0.5,'EdgeColor','none')
        hold on 
        % g12 =  plot([TstarVecV1 TstarVecV1],[Ylowlim 160] ,'LineWidth',1.5,'Color',[0,0,0]/255,'LineStyle','--'); %grey  [224,130,20]/255
        hold on
        % Plot ALT  
        g10 =  plot(sol.x(:),sol.y(7,:) ,'LineWidth',2.75,'Color', ALTColor,'LineStyle','-'); %grey 153,112,171
        hold on
        g60 = scatter(DataTimeIndex,PatientALTDataVecPlot,60,ALTColor,'o','filled'); %,'20', [33,102,172]/255,'*');
        hold on
        
        % 0731 limits
         ylim([-5,450])
        % 2158 limits
%         ylim([-5,200])
        % 
        xlim([-0.5,59])
        yticks([0:100:400]);
        xticks([0:20:60]); 
        % Print participant ID and dose 
        % Dose
        if DoseVector(ii) < 600
        formatSpec3 =   ['%dmg QD'] ; 
        DoseString = sprintf(formatSpec3, DoseVector(ii) );
        else  
        DoseString = '300mg BID';
        end   
        formatSpec2 =   ['ID: %d'] ; 
        IDString = sprintf(formatSpec2, PatientID(ii) );
%         text(38.2,375,IDString,'FontSize',14); % text(40,8.9,IDString,'FontSize',14); 
        if TrialToSimulate == 0731; 
            text(38.2,375,IDString,'FontSize',20); % text(40,8.9,IDString,'FontSize',14); 
            formatSpec2 =   ['CAM: Vebicorvir'] ; 
            TrialString = sprintf(formatSpec2);
            % text(35,400,TrialString,'FontSize',14); % text(35,9.5,TrialString,'FontSize',14);
            text(DoseTextLocation,400,DoseString,'FontSize',20); 
        elseif TrialToSimulate == 2158;
            text(38.2,160,IDString,'FontSize',14); % text(40,8.9,IDString,'FontSize',14); 
            formatSpec2 =   ['CAM: ABI-H2158'] ; 
            TrialString = sprintf(formatSpec2);
            text(35,175,TrialString,'FontSize',14); % text(35,9.5,TrialString,'FontSize',14);
            text(DoseTextLocation,175,DoseString,'FontSize',14); 
        end


        

        % title(PatientID(ii)) 
    end


end
 

function [sol] = HBVMultiscaleModel(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2);
sol = ode45(@HBVMultiscaleDynamics,IC,totaltime,opts);
    function dydt = HBVMultiscaleDynamics(t,y);   % [TIC,IIC,PIC,QIC,PA.VIC,PA.RIC,PA.AIC];
        dydt(1) = PA.lambda - PA.beta*y(1)*y(5) - PA.dT*y(1) ; %Differential equation for TIC
        dydt(2) = PA.beta*y(1)*y(5) - PA.delta*y(2);
        dydt(3) = PA.alpha*(1-CAMTreatment(PA,t))*y(2) - (PA.rhoR+PA.pi*(1-PA.epsN)+PA.delta+PA.muR)*y(3);
        dydt(4) = PA.pi*(1-PA.epsN)*y(3) - (PA.rhoV+PA.delta+PA.muV)*y(4); % PA.beta*y(1)*y(5)
        dydt(5) = PA.rhoV*y(4) - PA.cV*y(5);
        dydt(6) = PA.rhoR*y(3) - PA.cR*y(6);
        dydt(7) = PA.s + PA.alphaA*y(2) - PA.cA*y(7);
        dydt = dydt';
    end

end

function [Output] = CAMTreatment(PA,t) % Third derivative of the fitness function G
if t > PA.TreatmentEnd
    Output =  PA.epsC*exp( -PA.k*(t-PA.TreatmentEnd) )./( PA.epsC*( exp( -PA.k*(t-PA.TreatmentEnd) ) -1) +1 ) ;
else
    Output = PA.epsC ;
end
end

