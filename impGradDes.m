%% This code is used to estimate the weights of the MultiSpectral (MS) bands based on Adaptive IHS (AIHS) method.
%% References 
% [1]   S. Rahmani, M. Strait, D. Merkurjev, M. Moeller, and T. Wittman, 
%       "An adaptive IHS pan-sharpening method," IEEE Geoscience and Remote Sensing Letters 7, no. 4, 746-750, 2010.
% [2]   A. Azarang, H. E. Manoochehri and N. Kehtarnavaz, 
%       "Convolutional Autoencoder-Based Multispectral Image Fusion," in IEEE Access, vol. 7, pp. 35673-35683, 2019.
%% This code need two inputs: 
%           M         -   Low Resolution MS (LRMS) image to size of PANchromatic (PAN) image
%           P          -   Original PAN image
%% In the outputs you can find the estimated weight for each MS band
%       findalph  -   Spectral weights


function findalph = impGradDes(M, P)

[n, m, d] = size(M);

%% Initializing the optimal weight vector

findalph = ones(d,1);
%% Optimization process
for i=1:d
   for j=1:d
   A(i,j) = sum(sum(M(:,:,i).*M(:,:,j)));
   end
   B(i,1) = sum(sum(P.*M(:,:,i)));
end

tau = 5;
iter = 150000;
gamma1 = 1/200000;
gamma2 = 1;

inv = (eye(d) + 2*tau*gamma1*A)^(-1);

for i = 1:iter
   findalph = inv * (findalph+2*tau*max(-findalph,0)+2*tau*gamma1*B);
end
%% EOF