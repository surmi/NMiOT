function [unph,dIn] = TIE_unwrap2DPdct(ph,dx,lambda)

% This function carries out the unwrapping of a two-dimensional wrapped phase. 
% 
% Release Date: 15. August 2017 
% ***************************************************************
% Matlab Syntax: 
% [unph,dIn] = TIE_unwrap2DPdct(ph,dx,lambda); 
% 
% Input parameters:
% ph  				 a two-dimensional wrapped phase (matrix)
% dx                 pixel size 
% lambda             employed wavelength 
% 
% Output parameters:
% unph  		    a two-dimensional real unwrapped phase (matrix)
% dIn               a two-dimensional function that corresponds to the axial derivative of the intesity
% 
% ***************************************************************
% Copyright (c) 2017
% Juan Martinez-Carranza and Tomasz Kozacki
% Institute of Micromechanics and Photonics, 
% Division Engineering Photonics
% Warsaw University of Technology, 
% ul. sw. Andrzeja Boboli 8, 02-525 Warsaw, Poland
%
% e-mail: j.martinez@mchtr.pw.edu.pl
%  
% ***************************************************************
% Copyright (c) 2017, Juan Martinez-Carranza and Tomasz Kozacki
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% + Redistributions of source code must retain the above copyright notice,
%   this list of conditions and the following disclaimer.
% 
% + Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in the
%   documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% **************************************************************************/
% Description:
% 
% The function TIE_unwrap2DPdct carries out the unwrapping operation of a two-dimensional 
% wrapped phase and stores the result in the matrix unph. The unwrapping operation is 
% based on the Transport of Intensity Equation (TIE). The TIE links the phase and the intensity 
% by the Poisson Equation. When solving the Poisson equation 2D, the absolute phase values can
% be retrieved. This TIE-based unwrapping PROCEDURE is described in the paper:
% 
% - J. Martinez-Carranza, K. Falaggis, and T. Kozacki, "Fast and accurate phase-unwrapping 
% algorithm based on the transport of intensity equation," Appl. Opt. 56, 7079–7088 (2017). 
%
% Note that the implementation of the algorithm in this function
% uses the Discrete Cosine Transform instead of the mirror padding to solve
% the boundary condtion.

k = 2*pi/lambda;   
    
php=ph;                               % input wrapped phase
[Ny,Nx] = size(php);
dfx = 1/Nx/dx;  dfy = 1/Ny/dx;    
fx = (0  : (Nx -1))*dfx;  % for dct
fy = (0  : (Ny -1))*dfy;   % for dct
[Fx,Fy] = meshgrid(fx,fy);
F = Fx.^2+Fy.^2;                      % Laplacian in Frequency coordinates


%fist derivative calculation

Uin = (exp(1i*php));
Kernel =( -1i*pi* lambda*F ); 
du = idct2(Kernel.*(dct2(Uin)));
dIn = k*2*real (conj(Uin).*(du)); % Calculation of the axial derivative of the intensity (Eq. (12))
sum(sum(dIn/k))

F(1,1)=Inf; % for dct
sf2c = (1./F);  %for dct

unph = ( idct2( dct2(dIn).*sf2c ) )*1/4/pi^2;  % Solving the Poisson Equation (Eq. (14))



