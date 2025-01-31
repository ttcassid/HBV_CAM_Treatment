%Script to simulate parameter sensitivity to data
% This script runs multiple steps of the predictor corrector algorithm but
% not continuation in the classic sense of solving the nonlinear system to
% find the optimal next point.

close all
clear all

%rng default % For reproducibility
rng(4)
% %Colors
% ColorVec = zeros(2,3);
ColorVec(1,:) = [158,1,66]/256;
ColorVec(2,:) = [213,62,79]/256;
ColorVec(3,:) = [244,109,67]/256;
ColorVec(4,:) = [158,1,66]/256; ; %[ 253,174,97]/256;
ColorVec(5,:) = [158,1,66]/256; ; % [ 254,224,139]/256;
ColorVec(6,:) = [53,151,143]./256; % [ 230,245,152]/256;
ColorVec(7,:) = [ 171,221,164]/256;
ColorVec(8,:) = [ 102,194,165]/256;
ColorVec(9,:) = [50,136,189 ]/256;
ColorVec(10,:) = [171,221,164]/256 ;
ColorVec(11,:) = [171,221,164]/256 ;
ColorVec(12,:) = [171,221,164]/256;
ColorVec(13,:) = [171,221,164]/256;
ColorVec(14,:) = [171,221,164]/256;

%% Simulation time
t0 = 0;
tf =  28;
TotalTime = [t0 tf];

%% Load data
PatientHBVDNADataStruc = load('H0731_Study101_HBVDNA_RowDataSet_FullPt'); % load('VL_CD4_DataFile_Matrix_Censored');
PatientSpecificHBVDNAData = cell2mat(struct2cell(PatientHBVDNADataStruc));

PatientHBVRNADataStruc = load('H0731_Study101_HBVRNA_RowDataSet_FullPt'); % load('VL_CD4_DataFile_Matrix_Censored');
PatientSpecificHBVRNAData = cell2mat(struct2cell(PatientHBVRNADataStruc));

PatientID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];  % [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];

%% Load parameters from monolix fitting for all participants

ExcludedPartipants = 0;
ParameterStruc = load('H0731_IndPatientParameters_12Dec24_Rebound'); % load('H0731_IndPatientParameters_13June24'); %
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));

PatientDataID = PatientID;
PatientParameterID = [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];
% Find patients that have both data and parameters
[IntersectionDataParameters,ia,ib] = intersect(PatientParameterID,PatientDataID);
PatientID = PatientParameterID(ia);

%% Set up simulation
disp('Convert V_0 into copies/mL from IU/mL')


%% Pick which data point to perturb
ExperimentalTimeVec =[21]; %
RNAvsDNAPerturb = 2 % 1 -> DNA, 2 -> RNA.

if RNAvsDNAPerturb == 1; % only perturb DNA data 
                load('14-Dec-2024_HBV_dThetaAnalysis_DNA_Results.mat')
                
            else RNAvsDNAPerturb == 2;% only perturb RNA data 
                load('14-Dec-2024_HBV_dThetaAnalysis_RNA_Results.mat') 
end

PatientNumber = 21;
    
if ExcludedPartipants ==1;
    %% Load parameters
    %     PatientNumber = ii;
    X = PatientParameterMatrix(PatientNumber,:);
    % Fit parameters all need to be on the same scale
    X([1,3:6]) = log10(X([1,3:6]));
    %% PatientData
    %Get patient specific data
    PatientNumberData = PatientNumber;
    PatientDNADataVecPlot = PatientSpecificHBVDNAData(PatientNumberData,2:end) ;
    PatientRNADataVecPlot = PatientSpecificHBVRNAData(PatientNumberData,2:end) ;
elseif ExcludedPartipants == 0;
    PatientNumberParameter = ia(PatientNumber); % PatientID(ii) ;
    PatientNumberData = ib(PatientNumber);
    PatientNumber = PatientParameterMatrix(PatientNumberParameter,1);
    X = PatientParameterMatrix(PatientNumberParameter,2:end);
    % Fit parameters all need to be on the same scale
    X([1,3:6]) = log10(X([1,3:6]));
    PatientDNADataVecPlot = PatientSpecificHBVDNAData(PatientNumberData,2:end) ; % only load treatment data
    PatientRNADataVecPlot = PatientSpecificHBVRNAData(PatientNumberData,2:end) ;
end

%% Initialize the data structures, result output files, and data sets
X0 = X([1:6]);
Y0 = X([8:12,14:16]); % Fixed parameters
% Fix the HBV data.
D = [PatientDNADataVecPlot,PatientRNADataVecPlot]; % [PatientSpecificHBVDNAData(PatientNumberData,2:end),PatientSpecificHBVRNAData(PatientNumberData,2:end)];
% Set up perturbed data
PerturbPercentage = 0.9; % what percentage is the pertrubed data point of the original data point.
PerturbStepSize =  [-log10(PerturbPercentage),log10(PerturbPercentage)]; % [-0.05,0.05];
NPoints = 1; % length(PerturbStepSize) ;
PA.TreatmentEnd = 29;
%% Set up predictor 
  
%% Loop
for kk =  1: 1; % length(PatientID);

    PatientNumberData
    Y = X([8:12,14:16]); % Fixed parameters
    X0 = X([1:6]);
    tic
    %% Solve the Treated ODE system
    PA.dT = 0.004;
    PA.lambda = 13000000*PA.dT;
    
    % Fixed model parameters
    PA.cR = Y(1); %Fixed
    PA.cV = Y(2); %Fixed
    PA.muR = Y(3); %Fixed
    PA.muV = Y(4); %Fixed
    PA.epsN =  Y(5); %Fixed
    PA.AIC = 10^(Y(6)); %Fixed
    PA.AN =  Y(7); %Fixed
    PA.cA = Y(8);  %Fixed
    
    %Patient specific parameters
    PA.epsC = min(10^(X(1)),1);
    PA.beta = 10^(X(2));
    PA.alpha = 10^(X(3));
    PA.delta = 10^(X(4));
    PA.pi =  10^(X(5));
    PA.Rho = 10^(X(6));
    PA.rhoR = PA.Rho;
    PA.rhoV = PA.Rho;
    %% Initial condtions
    M = (PA.pi*PA.alpha)/(PA.delta*(PA.muV + PA.delta + PA.rhoV)*(PA.muR + PA.delta + PA.pi + PA.rhoR) );
    TIC = PA.cV./(PA.rhoV.*PA.beta.*M) ;
    IIC = PA.lambda/PA.delta - (PA.cV*PA.dT)/(PA.beta*PA.delta*PA.rhoV*M);
    PIC = (PA.alpha/(PA.muR + PA.delta + PA.pi + PA.rhoR) )*IIC;
    RIC = PA.rhoR.*PIC/PA.cR;
    QIC = PA.lambda.*M-PA.cV.*PA.dT./(PA.beta.*PA.rhoV);
    PA.VIC = PA.lambda.*PA.rhoV.*M./PA.cV-PA.dT/PA.beta;
    %Fixed parameters
    PA.s = 10^(PA.AN)*PA.cA;
    PA.alphaA = PA.cA*( PA.AIC-PA.s/PA.cA)./IIC;
    
    %% Solve the ODE systems
    IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
    
    [sol] =  HBVMultiscaleModel(IC,TotalTime,PA);
    
    PerturbIndex = 5 % change the perturbation time.
    
    for jj = 1:1; %  NPoints % for both perturbation directions
        
        if RNAvsDNAPerturb == 1; % only perturb DNA data
            PerturbedRNAData = PatientRNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]) ; % PatientSpecificHBVRNAData(PatientNumberData,[2:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]);
            PerturbedDNAData = [PatientDNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end]),...
                (1+PerturbStepSize(jj)).*PatientDNADataVecPlot(PerturbIndex)];
            % load('14-Dec-2024_HBV_dThetaAnalysis_DNA_Results.mat')
            
        else RNAvsDNAPerturb == 2;% only perturb RNA data
            PerturbedDNAData =  PatientDNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]) ;
            PerturbedRNAData = [PatientRNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end]),...
                (1+PerturbStepSize(jj)).*PatientRNADataVecPlot(PerturbIndex)];
            % load('14-Dec-2024_HBV_dThetaAnalysis_RNA_Results.mat')
        end 
        %% Perturb data
        FixedDNAData = PatientDNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]);
        FixedRNAData = PatientRNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]);
        D =  [FixedDNAData,FixedRNAData];
        TimeVec = [ExperimentalTimeVec([1:PerturbIndex-1,PerturbIndex+1:end]); ExperimentalTimeVec(PerturbIndex)];
        D1 =  [PerturbedDNAData,PerturbedRNAData] ;
        %% Update parameter
        X = X0'.*(RelativeParameterChangeTime5(:,PatientNumberParameter) ) + X0';
        
        %Patient specific parameters
        PA.epsC = min(10^(X(1)),1);
        PA.beta = 10^(X(2));
        PA.alpha = 10^(X(3));
        PA.delta = 10^(X(4));
        PA.pi =  10^(X(5));
        PA.Rho = 10^(X(6));
        PA.rhoR = PA.Rho;
        PA.rhoV = PA.Rho;
        %% Initial condtions
        M = (PA.pi*PA.alpha)/(PA.delta*(PA.muV + PA.delta + PA.rhoV)*(PA.muR + PA.delta + PA.pi + PA.rhoR) );
        TIC = PA.cV./(PA.rhoV.*PA.beta.*M) ;
        IIC = PA.lambda/PA.delta - (PA.cV*PA.dT)/(PA.beta*PA.delta*PA.rhoV*M);
        PIC = (PA.alpha/(PA.muR + PA.delta + PA.pi + PA.rhoR) )*IIC;
        RIC = PA.rhoR.*PIC/PA.cR;
        QIC = PA.lambda.*M-PA.cV.*PA.dT./(PA.beta.*PA.rhoV);
        PA.VIC = PA.lambda.*PA.rhoV.*M./PA.cV-PA.dT/PA.beta;
        
        %% Solve the ODE systems
        IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
        
        [solUpdate] =  HBVMultiscaleModel(IC,TotalTime,PA);
    end
    
    
end
toc
close all

%% HBV DNA Fig
DNAColor = [221,17,17]./255;
RNAColor = [0,0,139]./255;
DNAColorPerturb = [244,165,130]./255;
RNAColorPerturb = [146,197,222]./255; % [103,169,207]./255;
DNALLoQ = 0.95; % 20 IU/mL
RNALLoQ =  2.49; %% RNA LLoQ 250 copies/mL
DataTimeIndex = [0;1;7;14;21;28];  
    
Fig = figure(1); 
Ylowlim = -1.75; %0.0
YUpLim = 10;
TreatmentShadeX  = [0 PA.TreatmentEnd+0.25 PA.TreatmentEnd+0.25  0];
TreatmentShadeY = [Ylowlim Ylowlim YUpLim YUpLim];
fill(TreatmentShadeX,TreatmentShadeY,[205,205,205]./255,'FaceAlpha',0.45,'EdgeColor','none')
hold on
% g2 =  plot([t0 tf],[RNALLoQ RNALLoQ],'LineWidth',1.75,'Color', RNAColor,'LineStyle','--'); %grey
% hold on
% g2 =  plot([t0 tf],[DNALLoQ DNALLoQ],'LineWidth',1.75,'Color', DNAColor,'LineStyle','--'); %grey
% hold on 
yyaxis left 
g10 =  plot(sol.x(:),log(sol.y(5,:)./5.82)./log(10),'LineWidth',2.75,'Color', DNAColor,'LineStyle','-'); %grey 153,112,171
hold on
g10 =  plot(solUpdate.x(:),log(solUpdate.y(5,:)./5.82)./log(10),'LineWidth',2.75,'Color', DNAColorPerturb,'LineStyle','--'); %grey 153,112,171
hold on
CensoredIndexDNA = (PatientDNADataVecPlot <= DNALLoQ); % Find censored data points
if max(CensoredIndexDNA) > 0
    g60 = scatter(DataTimeIndex(CensoredIndexDNA ==0) ,PatientDNADataVecPlot(CensoredIndexDNA ==0),60,DNAColor,'o','filled'); %,'20', [33,102,172]/255,'*');
    hold on
    g61 = scatter(DataTimeIndex(CensoredIndexDNA > 0),PatientDNADataVecPlot(CensoredIndexDNA > 0),60,DNAColor,'o'); %,'20', [33,102,172]/255,'*');
    hold on
else
    g60 = scatter(DataTimeIndex,PatientDNADataVecPlot,60,DNAColor,'o','filled'); %,'20', [33,102,172]/255,'*');
    hold on
end
ylim([-2.0,10.25])
xlim([-0.5,28])
yticks([0:2.5:10]);
xticks([0:20:60]);
% ylabel('log_{10} HBV DNA (copies/mL) ','FontSize',15);
%% HBV RNA Fig
yyaxis right
g1 =  plot(sol.x(:),log(sol.y(6,:))./log(10),'LineWidth',2.75,'Color', RNAColor,'LineStyle','-'); %grey
hold on
g1 =  plot(solUpdate.x(:),log(solUpdate.y(6,:))./log(10),'LineWidth',2.75,'Color', RNAColorPerturb,'LineStyle','--'); %grey
hold on
CensoredIndexRNA = (PatientRNADataVecPlot <= RNALLoQ); % Find censored data points
if max(CensoredIndexRNA) > 0
    g6 = scatter(DataTimeIndex(CensoredIndexRNA ==0),PatientRNADataVecPlot(CensoredIndexRNA ==0),60,RNAColor,'o','filled'); %,'20', [33,102,172]/255,'*');
    hold on
    g6 = scatter(DataTimeIndex(CensoredIndexRNA > 0),PatientRNADataVecPlot(CensoredIndexRNA >0),60,RNAColor,'o'); %,'20', [33,102,172]/255,'*');
    hold on
else
    g6 = scatter(DataTimeIndex,PatientRNADataVecPlot,60,RNAColor,'o','filled'); %,'20', [33,102,172]/255,'*');
    hold on
    DataTimeIndexPerturb = [0,1,7,14,28,21];
    g6 = scatter(DataTimeIndexPerturb(end),PerturbedRNAData(end),55,RNAColorPerturb,'d','filled'); %,'20', [33,102,172]/255,'*');
end
hold on
ylim([-2.0,10.25])
xlim([-0.5,28])
yticks([0:2.5:10]);
title(PatientSpecificHBVDNAData(PatientNumberData,1)); 
ax = findall(Fig,'Type','axes');
ax.YAxis(1).Color = DNAColor;
ax.YAxis(2).Color = RNAColor;

 

%%

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




