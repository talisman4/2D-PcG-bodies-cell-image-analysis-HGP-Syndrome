% Inizializing u minimizer function 

function [pfu] = InitFunU(img, iNy, iNx) 

pfu  = abs(img);
maxI = double(max(img(:)));
minI = double(min(img(:)));
pfu  = (pfu-minI)/(maxI - minI); % pfu has values in [0,1] 
