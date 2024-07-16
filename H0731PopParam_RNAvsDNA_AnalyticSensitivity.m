%% Script to calculate the analytic sensitivity in circulating hBV RNA and HBV DNA as a function of CAM efficacy 

%Load population parameters
HBeAgPos = 1 % For HBeAg Positive or Negative
if HBeAgPos == 1 
    PA.epsC = 0.963; % [0.919,0.963,0.988]
    PA.beta = 9.9e-11; % 6.8e-7;
    PA.alpha = 365.6; %  0.203 
    PA.delta = 0.022; % 0.062
    PA.pi =  205.6;
    PA.rhoR = 2.37;
    PA.rhoV =  1.87;
    PA.cR = 1; %Fixed
    PA.cV = 1;
    PA.muR = 0;
    PA.muV = 0;
    PA.epsN =  0; %0.9962; %  X(8);
    PA.VIC =  6e8; %Conversion factor from IU/mL to copies/mL
    PA.AIC = 38.4;
    PA.AN =  23.2;
    PA.cA = 0.067; %Fixed 
elseif HBeAgPos == 0

     % HBeAg Positive population parameters
    PA.epsC = 0.963; % [0.919,0.963,0.988]
    PA.beta =  6.8e-7;
    PA.alpha =  0.203;
    PA.delta =  0.062;
    PA.pi =  205.6;
    PA.rhoR = 2.37;
    PA.rhoV =  1.87;
    PA.cR = 1; %Fixed
    PA.cV = 1;
    PA.muR = 0;
    PA.muV = 0;
    PA.epsN =  0; %0.9962; %  X(8);
    PA.VIC =  6e8; %Conversion factor from IU/mL to copies/mL
    PA.AIC = 38.4;
    PA.AN =  23.2;
    PA.cA = 0.067; %Fixed 
end

    PA.xiP = PA.pi+PA.rhoR+PA.delta+PA.muR;
    PA.xiQ = PA.rhoV+PA.delta+PA.muV;
    
    %Increased efficacy using estimates for 100, 200, 300 mg QD.
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







