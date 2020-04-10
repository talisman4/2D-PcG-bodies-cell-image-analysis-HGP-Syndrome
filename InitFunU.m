% This Matlab function is part of the code package released for the article:
%
% "A new high-throughput sequencing-based technology reveals early alteration
% of heterochromatin and deregulation of bivalent genes in Hutchinson-Gilford
% Progeria Syndrome"
% E. Sebestyén1, F. Marullo, F. Lucini, A. Bianchi, C. Petrini, S.  Valsoni,
% I. Olivieri, L. Antonelli, F. Gregoretti, G. Oliva, F. Ferrari and C. Lanzuolo
%
% We kindly request you to acknowledge the authors properly (citation
% or request for permission from the authors) when using this function.
%
% 2019 (C) L. Antonelli, F. Gregoretti, G. Oliva
%
% The package is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% The package is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.

% Inizializing u minimizer function 

function [pfu] = InitFunU(img, iNy, iNx) 

pfu  = abs(img);
maxI = double(max(img(:)));
minI = double(min(img(:)));
pfu  = (pfu-minI)/(maxI - minI); % pfu has values in [0,1] 
