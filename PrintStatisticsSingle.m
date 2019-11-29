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
