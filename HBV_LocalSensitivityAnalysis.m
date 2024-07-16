% Local sensitivity analysis of the HBV multiscale model
close all
clear all
%% Load parameters

ParameterStruc = load('H0731_IndPatientParameters_13June24');   
PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
PatientID =   [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];  

PatientNumber = 1;

PA.dT = 0.004;
PA.lambda = 13000000*PA.dT;
PA.TreatmentEnd = 28; 
%% Choose parameters to vary 
ParameterNames = {'epsC' 'beta' 'alpha'  'delta' 'pi'  'rhoR' 'rhoV'}; %
NumberOfParameters = length(ParameterNames);

VLNadirLowerParameter = zeros(length(PatientID),NumberOfParameters);
VLNadirUpperParameter = zeros(length(PatientID),NumberOfParameters);
VLNadirBaseline  = zeros(1,length(PatientID));

ReboundTimeLowerParameter = zeros(length(PatientID),NumberOfParameters);
ReboundTimeUpperParameter = zeros(length(PatientID),NumberOfParameters);
TimeToReboundBaseline = zeros(1, length(PatientID));

    %% Simulation time
    t0 = 1;
    tf =  56;
    TotalTime = [t0 tf]; 

for ii = 1: length(PatientID);
    %% Load parameters
    PatientNumber = ii;
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
    PA.epsN =  0; 
    PA.AIC = 10^(X(14));
    PA.AN =  X(15);
    PA.cA = X(16); %Fixed
    
    
    PA.xiP = PA.pi+PA.rhoR+PA.delta+PA.muR;
    PA.xiQ = PA.rhoV+PA.delta+PA.muV;
    %% Set up sampling procedure
    ParameterNames = {'epsC' 'beta' 'alpha'  'delta' 'pi'  'rhoR' 'rhoV' 'cR' 'cV'}; %
    ParameterHomeostasis = [PA.epsC ;PA.beta ; PA.alpha; PA.delta; PA.pi; PA.rhoR; PA.rhoV; PA.cR; PA.cV];
    NumberOfParameters = length(ParameterNames);
    PA.ReboundVL = 0.85.*PA.VIC;
    
    % Create lower and upper paraemter values
    P = 0.10; %Percentage of parameter change
    LowerBoundVector = (1-P).*ParameterHomeostasis;
    LowerBoundVector(2) =  PA.beta.*(1-P);
    UpperBoundVector = (1+P).*ParameterHomeostasis;
    UpperBoundVector(1) =   min(1,PA.epsC.*(1+P)) ; % Cannot have over 100% effectiveness.
    UpperBoundVector(2) =   PA.beta.*(1+P);
    
    %% Baseline simulation
    % Initial condtions
    
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
    %  Solve the ODE systems
    IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];

    % Solve the Treated ODE system for baseline parameters
    [solBaseline] =  HBVMultiscaleModel(IC,TotalTime,PA);
            
         VLNadirBaseline(ii) = min(solBaseline.y(5,:));   
         
        if ~isempty(solBaseline.xe)
            TimeToReboundBaseline(ii) = solBaseline.xe(1,1);
        else
            TimeToReboundBaseline(ii) = inf;
        end
    
    for kk = 1:length(ParameterHomeostasis)
        for jj = 1:kk
            PA.(ParameterNames{jj}) = ParameterHomeostasis(jj);
        end
        %% Simulate model for lower parameter values
        PA.(ParameterNames{kk}) = LowerBoundVector(kk);
         
        % Recalculate initial condtions 
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
        %  Solve the ODE systems
        IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
        % Solve the Treated ODE system
        [solLower] =  HBVMultiscaleModel(IC,TotalTime,PA);

                
         VLNadirLowerParameter(ii,kk) = min(solLower.y(5,:)); 

        if ~isempty(solLower.xe)
            ReboundTimeLowerParameter(ii,kk) = solLower.xe(1,1);
        else
            ReboundTimeLowerParameter(ii,kk) = inf;
        end
        
        
        %% Upper parameter simulation
        PA.(ParameterNames{kk}) = UpperBoundVector(kk);
        
        % Initial condtions
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
        %  Solve the ODE systems
        IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
        % Solve the Treated ODE system
        [solUpper] =  HBVMultiscaleModel(IC,TotalTime,PA);
        
        VLNadirUpperParameter(ii,kk) = min(solUpper.y(5,:)); 
        
        if ~isempty(solUpper.xe)
            ReboundTimeUpperParameter(ii,kk) = solUpper.xe(1,1);
        else
            ReboundTimeUpperParameter(ii,kk) = inf;
        end
        
    end
    
    
    
    
end

%% Aggregate data plotting
AggregateReboundTimeUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateReboundTimeLowerParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLNadirUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLNadirLowerParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLOneWeekUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLOneWeekLowerParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLRNAOneWeekUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateVLRNAOneWeekLowerParameterSensitivity = zeros(1,NumberOfParameters );
AggregateReSensitiseTimeUpperParameterSensitivity = zeros(1,NumberOfParameters );
AggregateReSensitiseTimeLowerParameterSensitivity = zeros(1,NumberOfParameters ); 
for kk = 1:NumberOfParameters % for each parameter
    AggregateVLNadirLowerParameterSensitivity(kk) = nanmedian( [VLNadirLowerParameter(:,kk)./VLNadirBaseline(:)]);
    AggregateVLNadirUpperParameterSensitivity(kk) = nanmedian( [VLNadirUpperParameter(:,kk)./VLNadirBaseline(:)]);
    
    AggregateReboundTimeLowerParameterSensitivity(kk) = nanmedian( [ReboundTimeLowerParameter(:,kk)./TimeToReboundBaseline(:)]);
    AggregateReboundTimeUpperParameterSensitivity(kk) = nanmedian( [ReboundTimeUpperParameter(:,kk)./TimeToReboundBaseline(:)]);
   
end
 
MedianTimeToRebound = median(TimeToReboundBaseline)

% 
Fig1 = figure(1);
xvalues = {'-10%','+10%'};
yvalues = ParameterNames; 
VLNadirHM = heatmap(xvalues,yvalues,100.*[AggregateVLNadirLowerParameterSensitivity',AggregateVLNadirUpperParameterSensitivity'],'ColorbarVisible','off') ; %'Title','Tumour Burden (% of normal)');

% %Time to rebound duration
Fig2 = figure(2);
xvalues = {'-10%','+10%'};
yvalues = ParameterNames;  
ReboundDelayHM = heatmap(xvalues,yvalues,100.*[AggregateReboundTimeLowerParameterSensitivity',AggregateReboundTimeUpperParameterSensitivity'],'ColorbarVisible','off') ; %,'Title','Tumour Doubling Time (% of normal)');

function Loss = Treatment(t,PA); 
if t < PA.TreatmentEnd
    Loss =1-PA.epsC;
else
    Loss = 1;
end
end




function [sol] = HBVMultiscaleModel(totaltime,IC,PA) %ODE model with therapy
opts = odeset('RelTol',1e-6,'AbsTol',1e-6,'MaxStep',1e-2,'Events',@TimeToRebound);
sol = ode45(@HBVMultiscaleDynamics,IC,totaltime,opts);
    function dydt = HBVMultiscaleDynamics(t,y);   % [TIC,IIC,PIC,QIC,PA.VIC,PA.RIC,PA.AIC];
        dydt(1) = PA.lambda - PA.beta*y(1)*y(5) - PA.dT*y(1) ; %Differential equation for TIC
        dydt(2) = PA.beta*y(1)*y(5) - PA.delta*y(2);
        dydt(3) = PA.alpha*Treatment(t,PA)*y(2) - (PA.rhoR+PA.pi*(1-PA.epsN)+PA.delta+PA.muR)*y(3);
        dydt(4) = PA.pi*(1-PA.epsN)*y(3) - (PA.rhoV+PA.delta+PA.muV)*y(4); % PA.beta*y(1)*y(5)
        dydt(5) = PA.rhoV*y(4) - PA.cV*y(5);
        dydt(6) = PA.rhoR*y(3) - PA.cR*y(6);
        dydt(7) = PA.s + PA.alphaA*y(2) - PA.cA*y(7);
        dydt = dydt';
    end
function [value,isterminal,direction] = TimeToRebound(t,y) %event function to calculate period from neutrophil DDE
        
        Rebound = y(5) - PA.ReboundVL ;
        isterminal(1) = 0;   % Continue the integration
        direction(1) = 1;   % increasing direction only
        
        value = [Rebound]; %,Nadir];
        
    end

end

