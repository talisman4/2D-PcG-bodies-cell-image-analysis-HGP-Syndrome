% This Matlab function is part of the code package released for the
% article:
%
% "A new high-throughput sequencing-based technology reveals early
% deregulation of bivalent genes in Hutchinson-Gilford Progeria Syndrome"
% E. Sebestyén1, F. Marullo, F. Lucini, A. Bianchi, C. Petrini, S.  Valsoni,
% I. Olivieri, L. Antonelli, F. Gregoretti, G. Oliva, F. Ferrari
%
% We kindly request you to acknowledge the authors properly
% (citation or request for permission from the authors) when using this
% function.
%
% 2019 (C) L. Antonelli, F. Gregoretti, G. Oliva

% Inizializing u minimizer function 

function [pfu] = InitFunU(img, iNy, iNx) 

pfu  = abs(img);
maxI = double(max(img(:)));
minI = double(min(img(:)));
pfu  = (pfu-minI)/(maxI - minI); % pfu has values in [0,1] 
