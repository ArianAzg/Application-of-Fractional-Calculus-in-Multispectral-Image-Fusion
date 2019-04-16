function [refined_map] = FractionalDiff(primitive_map)


%%   Description 
%    This function applies an approximation of fractional-order
%    differentation on the primitive detail map extracted from the PAN image
%    and a linear combination of MS image. 

%%  Reference
%   [1] A. Azarang and H. Ghassemian, "Application of fractional-order differentiation
%       in multispectral image fusion," Remote Sens. Lett., vol. 9, no. 1,
%       pp. 91-100, Jan. 2018.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input

%   -  primitive_map       primitive detail map

TFval = 0.09;         % Trade-off parameter between the spatial and spectral 
%                 injection, higher means higher spectral information. This
%                 parameter is chosen empirically for each satellite. 


M =     (1/8*(1+(-TFval)+(-TFval)*(-TFval+1)/2))*...
        [((-TFval)*(-TFval+1))/2, 0, ((-TFval)*(-TFval+1))/2, 0, ((-TFval)*(-TFval+1))/2;
        0, -TFval, -TFval, -TFval, 0;
        ((-TFval)*(-TFval+1))/2, -TFval, 8, -TFval, ((-TFval)*(-TFval+1))/2;
        0, -TFval, -TFval, -TFval, 0;
        ((-TFval)*(-TFval+1))/2, 0, ((-TFval)*(-TFval+1))/2, 0, ((-TFval)*(-TFval+1))/2;];  %% Superimposed 5*5 Mask 

refined_map = imfilter(double(primitive_map),M,'conv','same');

end