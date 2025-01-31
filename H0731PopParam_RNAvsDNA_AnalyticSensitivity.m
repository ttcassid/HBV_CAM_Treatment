%% Script to calculate the analytic sensitivity in circulating hBV RNA and HBV DNA as a function of CAM efficacy 

%Load population parameters
%% Load parameters
%Population estimate  parameters
HBeAgPos = 0
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

    PA.xiP = PA.pi+PA.rhoR+PA.delta+PA.muR;
    PA.xiQ = PA.rhoV+PA.delta+PA.muV;
    
    %Increased efficacy  EpsVec = [0.9,0.95,0.99,0.995];
     EpsVec =  [0.919,0.963,0.988]; 
    
    PA.c = PA.cR;
    f = @(x) (1-exp(-PA.xiP.*x(2))).*(PA.c.*exp(PA.c.*x(2)) );     %  argument within the integral for d R(t)/R_0/d \eps_C d P(t)/P_0/d \eps_C as a function of time and efficacy (x(1)) 
    Upsilon = @(x) (exp(PA.xiQ.*x(2))-1) - PA.xiQ/(PA.xiQ-PA.xiP).*(exp( (PA.xiQ-PA.xiP).*x(2) )-1) ;  
    g = @(x) ( Upsilon(x).*exp(-PA.xiQ.*x(2)) ).*(PA.c*exp(PA.c.*x(2)) );   % argument within the integral for d V(t)/V_0/d \eps_C as a function of time (x(2)) and efficacy (x(1)) 
    RNASensitivityIntegrand = @(t) f([PA.epsC,t]) ; % argument within the integral for d R(t)/R_0/d \eps_C d P(t)/P_0/d \eps_C as a function of time only 
    DNASensitivityIntegrand = @(t) g([PA.epsC,t]);  % argument within the integral for d V(t)/V_0/d \eps_C as a function of time only
    
    % Calculate sensitivity integrals
    IntTime = [0 5];
    RNASensitivityIntegral = integral(@(t) f([PA.epsC,t]) ,IntTime(1),IntTime(2),'arrayvalued', true)*exp(-PA.c*IntTime(2));
    DNASensitivityIntegral = integral(DNASensitivityIntegrand,IntTime(1),IntTime(2),'arrayvalued', true)*exp(-PA.c*IntTime(2));

    RelativeSensitivityDay5 = RNASensitivityIntegral./DNASensitivityIntegral

    IntTime = [0 1];
    RNASensitivityIntegral = integral(@(t) f([PA.epsC,t]) ,IntTime(1),IntTime(2),'arrayvalued', true)*exp(-PA.c*IntTime(2));
    DNASensitivityIntegral = integral(DNASensitivityIntegrand,IntTime(1),IntTime(2),'arrayvalued', true)*exp(-PA.c*IntTime(2));

    RelativeSensitivityDay1 = RNASensitivityIntegral./DNASensitivityIntegral







