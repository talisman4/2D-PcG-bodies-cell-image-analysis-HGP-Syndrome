% This Matlab function is part of the code package released for the article:
%
% "A new high-throughput sequencing-based technology reveals early alteration
% of heterochromatin and deregulation of bivalent genes in Hutchinson-Gilford
% Progeria Syndrome"
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

function [prob,pdfEst] = PrintStatisticsSingle(real,karaa,lw,xl,yl,tit,bwreal,classes)
    binCtrs = linspace(min(real),max(real),classes);
    binWidth = binCtrs(2) - binCtrs(1);
    
    counts = hist([real'],binCtrs);
    prob = counts ./ (length(real) * binWidth);
    sum_prob = sum(prob);
    prob = prob./sum_prob;

    
    if min(real)<0
        min_real = min(real)-1;
    else
        min_real = 0;
    end%if
    if min(rand)<0
        min_rand = min(rand)-1;
    else
        min_rand = 0;
    end%if
    
    % Real
    %paramEsts = wblfit(real-min_real);
    xgrid = linspace(min(real),max(real),100);
    %pdfEst = wblpdf(xgrid-min_real,paramEsts(1),paramEsts(2));
    pdfEst = pdf(fitdist(real','kernel','width',bwreal),xgrid);
    pdfEst = mat2gray(pdfEst).*max(prob);
    hold('on');
    plot(xgrid,pdfEst,'Color',karaa./1.25,'LineWidth',lw);
    hold('off');
       
    title(tit,'Interpreter','LaTeX');
    %xlabel(xl)%,'Interpreter','LaTeX');
    %ylabel(yl)%,'Interpreter','LaTeX');
    xlabel(xl,'Interpreter','LaTeX','FontSize',14,'FontWeight','bold');
    ylabel(yl,'Interpreter','LaTeX','FontSize',14,'FontWeight','bold');
    %legend(l1,l2,'Location','Best');%NorthEast');
    %set(gca,'yaxislocation','right');
