%% This is the objective function for the fitter for PA.Rb and PA.Da=PA.Db

function [Obj] = HBVModel_ExperimentalPrediction_ObjectiveFunction(X,D,TimeVec,Y);
  
DNAViralLoadData = D(1:length(TimeVec));
RNAViralLoadData = D(length(TimeVec)+1:end);

VLTime = TimeVec;

PA.dT = 0.004;
PA.lambda = 13000000*PA.dT;

%Patient specific parameters (not log transformed)
% PA.epsC = X(1);
% PA.beta = 10^(X(2));
% PA.alpha = X(3);
% PA.delta = X(4);
% PA.pi =  X(5);
% PA.rhoR = X(6);
% PA.rhoV =  X(7);
% PA.VIC = 10^(X(8)).*5.82; %Conversion factor from IU/mL to copies/mL

%Patient specific parameters
PA.epsC = min(10^(X(1)),1);
PA.beta = 10^(X(2));
PA.alpha = 10^(X(3));
PA.delta = 10^(X(4));
PA.pi =  10^(X(5));
PA.rhoR = 10^(X(6));
PA.rhoV = 10^(X(7));
% PA.VIC = 10^(X(8)).*5.82; 

% Fixed model parameters
PA.cR = Y(1); %Fixed
PA.cV = Y(2); %Fixed
PA.muR = Y(3); %Fixed
PA.muV = Y(4); %Fixed
PA.epsN =  Y(5); %Fixed
PA.AIC = 10^(Y(6)); %Fixed
PA.AN =  Y(7); %Fixed
PA.cA = Y(8);  %Fixed

PA.xiP = PA.pi+PA.rhoR+PA.delta+PA.muR;
PA.xiQ = PA.rhoV+PA.delta+PA.muV;

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
%% Simulation time
t0 = 0;
tf =  28;
TotalTime = [t0 tf];

%% Solve the ODE systems
IC = [TIC,IIC,PIC,QIC,PA.VIC,RIC,PA.AIC];
%% Solve the Treated ODE system
[sol] =  HBVMultiscaleModel(IC,TotalTime,PA);
EvalSolHBVDNA = deval(sol,TimeVec,5);
EvalSolHBVRNA = deval(sol,TimeVec,6);
% [log10( EvalSolHBVDNA./5.82 ) , log10(EvalSolHBVRNA ) ]; %Convert copies to IU/mL.

HBVDNAModelObjFunc = log10( EvalSolHBVDNA(~isnan(DNAViralLoadData))./5.82 );
DNAVLDataObjFunc = DNAViralLoadData(~isnan(DNAViralLoadData)); % DNAViralLoadData(~isnan(DNAViralLoadData))-DNAViralLoadData(1) ; % 

HBVRNAModelObjFunc = log10(EvalSolHBVRNA(~isnan(RNAViralLoadData) ) );
RNAVLDataObjFunc = RNAViralLoadData(~isnan(RNAViralLoadData)); %  RNAViralLoadData(~isnan(RNAViralLoadData)) - RNAViralLoadData(1) ; % 

ObjSquare =   sum( HBVDNAModelObjFunc - DNAVLDataObjFunc ).^2  + sum( (HBVRNAModelObjFunc - RNAVLDataObjFunc ).^2 )   ; % enforce two penalties: Ra>Rb and Ra < DaMax ;
Obj = sqrt(ObjSquare);
 
%% Solvers
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



end