%% -------------------------------------------------------------------------
% This Matlab function is part of the code package released for the article:
%
% "SAMMY-seq reveals early alteration of heterochromatin and deregulation of
% bivalent genes in Hutchinson-Gilford Progeria Syndrome"
% E. Sebestyén, F. Marullo, F. Lucini, C. Petrini, A. Bianchi, S. Valsoni,
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
%
% Isodata Thresholding algorithm: Modified formulation defined using as
% the initial guess for image thresholding the mean of the average 
% of image intensity in the foreground and background region 
% computed by Chan - Vese Model 
%% -------------------------------------------------------------------------

function [thresh] = IsoDataThresh(frame_prime, iNy, iNx, C1, C2, thresh)

%  thresh_old = .5*(C1 + C2);
%  thresh_new = thresh_old; 

  if (C1 > C2)
    C = C1;
  else
    C = C2;
  end
  thresh_new = C;

  c_old = 0.;
  c_new = 0.;
  thresh_old = C;
  % c_new = abs(thresh_new - thresh_old);

  while ( (c_new - c_old) >= 0 )
      c_old = c_new;
      thresh_old = thresh_new;

      n1 = 0;
      n2 = 0;
      Iplus = 0.;
      Iminus = 0.;

      for ii=1:iNy
        for jj=1:iNx
          if ( frame_prime(ii,jj) ~= -255 )  
            if(frame_prime(ii,jj) >= thresh_old)  
              n1 = n1 + 1;;
              Iplus = Iplus + double(frame_prime(ii,jj));
            else
              n2 = n2 + 1;
              Iminus = Iplus + double(frame_prime(ii,jj));
            end %if
          end %if
        end%for
      end%for 
    
      if (n1 ~= 0)
        C1 = Iplus / double(n1);
      else
        C1 = 0.;
      end
      if (n2 ~= 0)
        C2 = Iminus / double(n2);
      else
        C2 = 0.;
      end
      % fprintf('C1 %f n1 %d C2 %f n2 %d \n',C1, n1, C2, n2)
      
      thresh_new  = .5*(C1 + C2);
      c_new = abs(thresh_new - thresh_old);
      % fprintf('c_old %f c_new %f thresh_new %f \n',c_old, c_new, thresh_new);;
  end
  thresh = thresh_new;
end %
