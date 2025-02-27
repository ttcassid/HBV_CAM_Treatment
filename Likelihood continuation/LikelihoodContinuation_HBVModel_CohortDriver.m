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

%% Load data
% PatientHBVDNADataStruc = load('H0731_Study101_HBVDNA_RowDataSet'); % load('VL_CD4_DataFile_Matrix_Censored');
% PatientSpecificHBVDNAData = cell2mat(struct2cell(PatientHBVDNADataStruc));
% 
% PatientHBVRNADataStruc = load('H0731_Study101_HBVRNA_RowDataSet'); % load('VL_CD4_DataFile_Matrix_Censored');
% PatientSpecificHBVRNAData = cell2mat(struct2cell(PatientHBVRNADataStruc));
% 
% PatientID = [3;5;6;7;8;9;10;11;13;15;16;18;19;20;21;24;26;27;28;29;31;32;33;34;35;36;37];

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

% ExcludedPartipants = 0;
% ParameterStruc = load('Study101_IndPatientParameters_FullParticipant_07Feb24') ; % load('Study101_IndPatientParameters_FullParticipant_25Jan2024');
% PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
% 
% 
% PatientDataID = PatientID;
% PatientParameterID =  [3;5;6;7;8;9;10;11;12;13;15;16;18;19;20;21;23;24;26;27;28;29;31;32;33;34;35;36;37];
% % Find patients that have both data and parameters
% [IntersectionDataParameters,ia,ib] = intersect(PatientParameterID,PatientDataID);
% PatientID = PatientParameterID(ia);

%% For results using finalized parameter set excluding 3 participants
% ExcludedPartipants = 1;
% ParameterStruc =  load('Study101_IndPatientParameters_30Aug2023'); % load('VL_CD4_DataFile_Matrix_Censored');
% PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
% PatientNumber = 1;

%% Set up simulation
disp('Convert V_0 into copies/mL from IU/mL')

% PatientParameterMatrix = cell2mat(struct2cell(ParameterStruc));
% PatientNumber = 10;


PatientNumber = 1;
if ExcludedPartipants ==1;
    %% Load parameters
    %     PatientNumber = ii;
    X = PatientParameterMatrix(PatientNumber,:);
    % Fit parameters all need to be on the same scale
    X([1,3:6]) = log10(X([1,3:6]));
    %% PatientData
    %Get patient specific data
    PatientNumberData = PatientNumber;
    PatientDNADataVecPlot = PatientSpecificHBVDNAData(PatientNumber,2:end) ;
    PatientRNADataVecPlot = PatientSpecificHBVRNAData(PatientNumber,2:end) ;
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
NPoints = length(PerturbStepSize) ;
%% Set up predictor
StepSize = 0.1;
NParam = length(X0);
%Structures necessary to calculate first order derivative
ParameterStencil = eye(NParam,NParam);
XLower = zeros(NParam,NParam);
XUpper = zeros(NParam,NParam);

FValLower = zeros(NParam,length(D));
FValUpper = zeros(NParam,length(D));
ObjMixedSecondDifferenceXD = zeros(NParam,length(D));

% Structures necessary to calculate Hessian matrix
BasisVec = eye(NParam);
k = 1; % Counter for each row.
NEvalsHessian = (NParam)*(NParam-1)/2;
FValPosTemplate = zeros(NEvalsHessian,NParam);
FValNegTemplate = zeros(NEvalsHessian,NParam);
for ii = 1 :NParam-1
    for jj = ii+1:NParam
        FValPosTemplate(k,:) = BasisVec(ii,:)+BasisVec(jj,:);
        FValNegTemplate(k,:) = -BasisVec(ii,:)+BasisVec(jj,:);
        k = k+1;
    end
end
PositiveStencil = [FValPosTemplate;-FValPosTemplate];
NegativeStencil = [FValNegTemplate;-FValNegTemplate];

% Function to identify the Stencil row for the i,j partial derivative;
% x(1)<x(2) necessarily here, only calculating the lower portion of the Hessian due to symmetry
FStencil = @(x) (x(1)-1)*(2*NParam - x(1)-2)./2 + x(2)-1;
% FStencilV1 = @(x) (x(1)-1)*(2*NParam - x(1))./2 + x(2);

%Matrices to store parameter updates and function evaluations
XPositive = zeros(NParam,NParam);
XNegative = zeros(NParam,NParam);
FValPositive = zeros(NParam,1);
FValNegative = zeros(NParam,1);
ObjHessian = zeros(NParam,NParam);

%% Set up output
UpdatedGuess = zeros(NPoints,length(X0));
XInitialSol = X0 ;
XUpdatedSol =   XInitialSol;
PredictedObjFunction = zeros(NPoints,1);
GuessObjFunction = zeros(NPoints,1); 
PredictionfEvals = zeros(NPoints,1);

% Vector to store all relative change (For each data point, individual patients get a row to track relative changes in parameters ( + perturbation in kk-th column, - perturbation in length (Patient Number) + kk for each parameter) 
%RelativeParameterChangeTime1 =  zeros(length(X0),2.*length(PatientID)); 
RelativeParameterChangeTime2 =  zeros(length(X0),2.*length(PatientID));  % zeros(length(PatientID),2.*length(X0)); 
RelativeParameterChangeTime3 =  zeros(length(X0),2.*length(PatientID));
RelativeParameterChangeTime4 =  zeros(length(X0),2.*length(PatientID)); 
RelativeParameterChangeTime5 =  zeros(length(X0),2.*length(PatientID));
RelativeParameterChangeTime6 =  zeros(length(X0),2.*length(PatientID));

ParameterFoldChangeVec = zeros(NPoints,length(X0)); % zeros(NPoints,length(X0)+1);
ParameterAbsChangeVec = zeros(NPoints,length(X0)); % zeros(NPoints,length(X0)+1);

% For plotting
% Fig1 = figure(1);
Fig2 = figure(2);

%% Pick which data point to perturb
ExperimentalTimeVec = [0;1;7;14;21;28];  
RNAvsDNAPerturb = 2 % 1 -> DNA, 2 -> RNA.

MaxFoldChange = 0;% For plotting the tornado plots on same x axis.

%% Loop
for kk =  1: length(PatientID); 
    PatientNumber = kk
    if ExcludedPartipants ==1;
        %% Load parameters
        %     PatientNumber = ii;
        X = PatientParameterMatrix(PatientNumber,:);
        % Fit parameters all need to be on the same scale
        X([1,3:6]) = log10(X([1,3:6]));
        %% PatientData
        %Get patient specific data
        PatientNumberData = PatientNumber;
        PatientDNADataVecPlot = PatientSpecificHBVDNAData(PatientNumber,2:end) ;
        PatientRNADataVecPlot = PatientSpecificHBVRNAData(PatientNumber,2:end) ;
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
    X0 = X([1:6]);
    Y0 = X([8:12,14:16]); % Fixed parameters
    XInitialSol = X0 ;
    XUpdatedSol =   XInitialSol;

    tic
    for nn =  2:length(ExperimentalTimeVec); % for each data point
        PerturbIndex = nn % change the perturbation time.

        for jj = 1: NPoints % for both perturbation directions

            if RNAvsDNAPerturb == 1; % only perturb DNA data
                PerturbedRNAData = PatientRNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]) ; % PatientSpecificHBVRNAData(PatientNumberData,[2:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]);
                PerturbedDNAData = [PatientDNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end]),...
                    (1+PerturbStepSize(jj)).*PatientDNADataVecPlot(PerturbIndex)];
            else RNAvsDNAPerturb == 2;% only perturb RNA data
                PerturbedDNAData =  PatientDNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]) ;
                PerturbedRNAData = [PatientRNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end]),...
                    (1+PerturbStepSize(jj)).*PatientRNADataVecPlot(PerturbIndex)];
            end

            %% Perturb data
            FixedDNAData = PatientDNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]);
            FixedRNAData = PatientRNADataVecPlot([1:PerturbIndex-1,PerturbIndex+1:end,PerturbIndex]);
            D =  [FixedDNAData,FixedRNAData];
            TimeVec = [ExperimentalTimeVec([1:PerturbIndex-1,PerturbIndex+1:end]); ExperimentalTimeVec(PerturbIndex)];
            D1 =  [PerturbedDNAData,PerturbedRNAData] ;
            %     break
            XUpdatedSol = X0;

            %% Calculate numerical gradient
            %Calculate mixed partial d/dD (d/dx G(x,D))
            for ii = 1 :NParam
                XLower(ii,:)    =  XUpdatedSol - StepSize.*ParameterStencil(ii,:);
                XUpper(ii,:)    =  XUpdatedSol + StepSize.*ParameterStencil(ii,:);
                FValLower(ii,:) =  HBVModel_ExperimentalPrediction_Solver(XLower(ii,:),D,TimeVec,Y0) ;
                FValUpper(ii,:) = HBVModel_ExperimentalPrediction_Solver(XUpper(ii,:),D,TimeVec,Y0);
                ObjMixedSecondDifferenceXD(ii,:) = -2.*(FValUpper(ii,:) - FValLower(ii,:))./(2.*StepSize);
            end

            % Create off diagonal elements of Hessian Matrix first
            %Calculate Diagonal elements of Hessian Matrix d^2/dx^2(G(x,D))
            %first, re-initalize to zero from previous runs
            ObjHessian = zeros(NParam,NParam);
            for ii = 1:NParam-1
                for mm = ii+1:NParam
                    HessianIndex = FStencil([ii,mm]);
                    XPositive(1,:) =  XUpdatedSol + StepSize.*PositiveStencil(HessianIndex,:); %.*XInitialSol;
                    XPositive(2,:) =  XUpdatedSol + StepSize.*PositiveStencil(NEvalsHessian+HessianIndex,:); %.*XInitialSol;
                    XNegative(1,:) =  XUpdatedSol + StepSize.*NegativeStencil(HessianIndex,:); %.*XInitialSol;
                    XNegative(2,:) =  XUpdatedSol + StepSize.*NegativeStencil(NEvalsHessian+HessianIndex,:); %.*XInitialSol;
                    for ll = 1:2
                        FValPositive(ll,:) = HBVModel_ExperimentalPrediction_ObjectiveFunction(XPositive(ll,:),D,TimeVec,Y0);
                        FValNegative(ll,:) = HBVModel_ExperimentalPrediction_ObjectiveFunction(XNegative(ll,:),D,TimeVec,Y0);
                    end
                    ObjHessian(ii,mm) = ( FValPositive(1,:)+FValPositive(2,:) - ( FValNegative(1,:)+FValNegative(2,:) ) )./(4*StepSize.^2);
                end
            end
            % As we haven't calculated diagonal elements yet and enforce symmetry
            ObjHessian = ObjHessian + ObjHessian';

            %Calculate Diagonal elements of Hessian Matrix d^2/dx^2(G(x,D))
            for ii = 1:NParam
                ObjHessian(ii,ii) = (HBVModel_ExperimentalPrediction_ObjectiveFunction(XUpper(ii,:),D,TimeVec,Y0) ...
                    - 2.*HBVModel_ExperimentalPrediction_ObjectiveFunction(XUpdatedSol,D,TimeVec,Y0)...
                    +HBVModel_ExperimentalPrediction_ObjectiveFunction(XLower(ii,:),D,TimeVec,Y0) )./(StepSize.^2);
            end

            %Calculate the predictor matrix
            dThetaDData = - ObjHessian\ObjMixedSecondDifferenceXD;

            %
            CorrectorDataVec = zeros(length(D1),1); % To account for the NaN elements in the data vector by replacing them with 0
            CorrectorDataVec(~isnan(D)>0) = (D1(~isnan(D))'-D(~isnan(D))');
            Corrector = dThetaDData*CorrectorDataVec; % (D1(~isnan(D))'-D(~isnan(D))') ; % log(D')./log(10))
            UpdatedGuess(jj,:) = XUpdatedSol + Corrector' ;
        end

        %% Plotting fold change in parameters
        for ii = 1:NPoints
            ParameterFoldChangeVec(ii,:) = [( UpdatedGuess(ii,:) -X0 )./X0 ];
            if isnan(D(length(D)/2*(RNAvsDNAPerturb-1)+nn) ) 
                ParameterFoldChangeVec(ii,:) = NaN;
            end
        end 

% All the positive changes then all the negative changes in column space of these matrices. 
        if nn == 2
           RelativeParameterChangeTime2([1:length(X0)],kk) = ParameterFoldChangeVec(1:NPoints/2,1:end);
           RelativeParameterChangeTime2([1:length(X0)],length(PatientID)+kk) = ParameterFoldChangeVec(1+NPoints/2:end,1:end);
        elseif nn == 3
           RelativeParameterChangeTime3([1:length(X0)],kk) = ParameterFoldChangeVec(1:NPoints/2,1:end);
           RelativeParameterChangeTime3([1:length(X0)],length(PatientID)+kk) = ParameterFoldChangeVec(1+NPoints/2:end,1:end);
        elseif nn == 4
           RelativeParameterChangeTime4([1:length(X0)],kk) = ParameterFoldChangeVec(1:NPoints/2,1:end);
           RelativeParameterChangeTime4([1:length(X0)],length(PatientID)+kk) = ParameterFoldChangeVec(1+NPoints/2:end,1:end);
        elseif nn == 5
           RelativeParameterChangeTime5([1:length(X0)],kk) = ParameterFoldChangeVec(1:NPoints/2,1:end);
           RelativeParameterChangeTime5([1:length(X0)],length(PatientID)+kk) = ParameterFoldChangeVec(1+NPoints/2:end,1:end);
        elseif nn == 6
           RelativeParameterChangeTime6([1:length(X0)],kk) = ParameterFoldChangeVec(1:NPoints/2,1:end);
           RelativeParameterChangeTime6([1:length(X0)],length(PatientID)+kk) = ParameterFoldChangeVec(1+NPoints/2:end,1:end);
        end

    end
end
%% Plotting the median changes

       %  MaxFoldChange = max(1.1*max(RelativeParameterChangeTime2(1,:)),0.1); % update the size of the x axis to be the largest necessary.
        XSpace = [1:NParam];
        TornadoColorVec1 =  [158,1,66]/256;
        TornadoColorVec2 =  [171,221,164]/256;
        ParameterId = {'\eps_C','\beta','\alpha','\delta','\pi','\rho_R','\rho_V'};
% Cycle through each perturbed data point
%% Data point 2
        PopulationParameterFoldChange = zeros(2,NParam);
        for zz = 1:NParam 
        PopulationParameterFoldChange(1,zz) = nanmedian(RelativeParameterChangeTime2(zz,1:length(PatientID)));
        PopulationParameterFoldChange(2,zz) = nanmedian(RelativeParameterChangeTime2(zz,1+length(PatientID):end) );
        end
        figure(2)
        TornadoValues = 100.*[PopulationParameterFoldChange]';
        b =  barh(XSpace,TornadoValues,'stacked'); % barh(XSpace(ii),TornadoValues,'stacked');
        b(1).FaceColor = TornadoColorVec1;
        b(2).FaceColor = TornadoColorVec2;
        yticks([XSpace]');
        yticklabels({});
        % yticklabels(ParameterId);
        % if RNAvsDNAPerturb ==1
        %     title('HBV DNA Perturbation')
        % elseif RNAvsDNAPerturb ==2
        %     title('HBV RNA Perturbation')
        % end

%% Data point 3
        PopulationParameterFoldChange = zeros(2,NParam);
        for zz = 1:NParam 
        PopulationParameterFoldChange(1,zz) = nanmedian(RelativeParameterChangeTime3(zz,1:length(PatientID)));
        PopulationParameterFoldChange(2,zz) = nanmedian(RelativeParameterChangeTime3(zz,1+length(PatientID):end) );
        end
        figure(3)
        TornadoValues = 100.*[PopulationParameterFoldChange]';
        b =  barh(XSpace,TornadoValues,'stacked'); % barh(XSpace(ii),TornadoValues,'stacked');
        b(1).FaceColor = TornadoColorVec1;
        b(2).FaceColor = TornadoColorVec2;
        yticks([XSpace]');
         yticklabels({});
        % yticklabels(ParameterId);
        % if RNAvsDNAPerturb ==1
        %     title('HBV DNA Perturbation')
        % elseif RNAvsDNAPerturb ==2
        %     title('HBV RNA Perturbation')
        % end
%% Data point 4
        PopulationParameterFoldChange = zeros(2,NParam);
        for zz = 1:NParam 
        PopulationParameterFoldChange(1,zz) = nanmedian(RelativeParameterChangeTime4(zz,1:length(PatientID)));
        PopulationParameterFoldChange(2,zz) = nanmedian(RelativeParameterChangeTime4(zz,1+length(PatientID):end) );
        end
        figure(4)
        TornadoValues = 100.*[PopulationParameterFoldChange]';
        b =  barh(XSpace,TornadoValues,'stacked'); % barh(XSpace(ii),TornadoValues,'stacked');
        b(1).FaceColor = TornadoColorVec1;
        b(2).FaceColor = TornadoColorVec2;
        yticks([XSpace]');
        yticklabels({});
        % yticklabels(ParameterId);
        % if RNAvsDNAPerturb ==1
        %     title('HBV DNA Perturbation')
        % elseif RNAvsDNAPerturb ==2
        %     title('HBV RNA Perturbation')
        % end
 %% Data point 5
        PopulationParameterFoldChange = zeros(2,NParam);
        for zz = 1:NParam 
        PopulationParameterFoldChange(1,zz) = nanmedian(RelativeParameterChangeTime5(zz,1:length(PatientID)));
        PopulationParameterFoldChange(2,zz) = nanmedian(RelativeParameterChangeTime5(zz,1+length(PatientID):end) );
        end
        figure(5)
        TornadoValues = 100.*[PopulationParameterFoldChange]';
        b =  barh(XSpace,TornadoValues,'stacked'); % barh(XSpace(ii),TornadoValues,'stacked');
        b(1).FaceColor = TornadoColorVec1;
        b(2).FaceColor = TornadoColorVec2;
        yticks([XSpace]');
        yticklabels({});
       % yticklabels(ParameterId);
        % if RNAvsDNAPerturb ==1
        %     title('HBV DNA Perturbation')
        % elseif RNAvsDNAPerturb ==2
        %     title('HBV RNA Perturbation')
        % end

%% Data point 6
        PopulationParameterFoldChange = zeros(2,NParam);
        for zz = 1:NParam
            PopulationParameterFoldChange(1,zz) = nanmedian(RelativeParameterChangeTime6(zz,1:length(PatientID)));
            PopulationParameterFoldChange(2,zz) = nanmedian(RelativeParameterChangeTime6(zz,1+length(PatientID):end) );
        end
        figure(6)
        TornadoValues = 100.*[PopulationParameterFoldChange]';
        b =  barh(XSpace,TornadoValues,'stacked'); % barh(XSpace(ii),TornadoValues,'stacked');
        b(1).FaceColor = TornadoColorVec1;
        b(2).FaceColor = TornadoColorVec2;
        yticks([XSpace]');
        yticklabels({});
       %  yticklabels(ParameterId );
        % if RNAvsDNAPerturb ==1
        %     title('HBV DNA Perturbation')
        % elseif RNAvsDNAPerturb ==2
        %     title('HBV RNA Perturbation')
        % end
for mm = 2:6 % length(ExperimentalTimeVec);
figure(mm);
xlim([-50 50]);
xticks([-50:20:50])
end

%%  
 t = datetime('today');
 S = char(t);
 if RNAvsDNAPerturb ==1
filename = [S,'_HBV_dThetaAnalysis_DNA_Results']; 
 elseif RNAvsDNAPerturb ==2
filename = [S,'_HBV_dThetaAnalysis_RNA_Results']; 
 end
 
save(filename)

toc


