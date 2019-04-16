clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Description 
%   This code is the MATLAB code for applying the fractional-order
%   differentatioan in multispectral image fusion. (FDIF)
%%  Reference
%   [1] S. Rahmani, M. Strait, D. Merkurjev, M. Moeller, and T. Wittman,
%       "An adaptive IHS Pan-sharpening method," IEEE Geosci. Remote Sens.
%       Lett., vol. 7, no. 4, pp. 746-750, Oct. 2010.
%   [2] A. Azarang and H. Ghassemian, "Application of fractional-order differentiation
%       in multispectral image fusion," Remote Sens. Lett., vol. 9, no. 1,
%       pp. 91-100, Jan. 2018.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the dataset

addpath QuickBird_Data   %% Dataset path
load  PAN;               %% loading the MS image
load   MS;                %% Loading the PAN image

%% Make the PAN and MS data ready for the processing

MSWV_db  = double(MS);
PANWV_db = double(PAN);
MS_ORG   = double(MS);

%% Resizing, Upsampling the MS data to the size of PAN

MSWV_US  = imresize(MSWV_db,  1/4, 'bicubic');
MSWV_US  = imresize(MSWV_US,  4,   'bicubic');
MSWV_DG  = MSWV_US;
PANWV_DS = imresize(PANWV_db, 1/4, 'bicubic');
PANWV_US = imresize(PANWV_DS, 4,   'bicubic');

figure, imshow(uint8(MSWV_DG(:,:,1:3)),'Border','tight')
figure, imshow(uint8(PANWV_DS),'Border','tight')

%% Data Normialization

for i=1:size(MSWV_US,3)
    bandCoeffs(i)      = max(max(MSWV_US(:,:,i)));
    MSWV_US(:,:,i)     = MSWV_US(:,:,i)/bandCoeffs(i);
end

P = PANWV_DS;
panCoeff = max(max(P));
P = P/panCoeff;

%% Primitive detail map

W  =  impGradDes(MSWV_US, P);
I  =  W(1).*MSWV_US(:,:,1)+W(2).*MSWV_US(:,:,2)+W(3).*MSWV_US(:,:,3)+W(4).*MSWV_US(:,:,4);  %% Optimal weights using AIHS
P  =  (P - mean(P(:)))*std(I(:))/std(P(:)) + mean(I(:));                                    %% Histogram matching

%% Proposed framework

refined_map = FractionalDiff(P-I);
Fused_FDIF(:,:,1) = MSWV_US(:,:,1) + (refined_map);
Fused_FDIF(:,:,2) = MSWV_US(:,:,2) + (refined_map);
Fused_FDIF(:,:,3) = MSWV_US(:,:,3) + (refined_map);
Fused_FDIF(:,:,4) = MSWV_US(:,:,4) + (refined_map);

%% Denormalization
for i=1:size(Fused_FDIF, 3)
    Fused_FDIF(:,:,i)  = Fused_FDIF(:,:,i)*bandCoeffs(i);
end


%% Showing the fusion result
figure, imshow(uint8(Fused_FDIF(:,:,1:3)),'Border','tight')
%% Objective assessment of FDIF
addpath Objective_Evaluation
Methods = {'FDIF'};
ERGAS = ERGAS(MS_ORG,Fused_FDIF, 4);
SAM   = SAM(MS_ORG,Fused_FDIF);
RASE  = RASE(MS_ORG,Fused_FDIF);
RMSE  = RMSE(MS_ORG,Fused_FDIF);
UIQI  = uqi(MS_ORG,Fused_FDIF);
CC    = CC(MS_ORG,Fused_FDIF);
T = table(ERGAS, SAM, RASE, RMSE, UIQI, CC, 'RowNames', Methods)

% End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%