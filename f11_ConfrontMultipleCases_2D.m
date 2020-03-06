%%% Begin Main Function
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
%
% f11_ConfrontMultipleCases_2D uses data produced by mainSeg_pcg2D.m to compute:
% - the number of PcG bodies per nucleus,
% - the area of any Nucleus,
% - the area of any PcG, 
% - the eccemtricity of any Nucleus;
% create the graphs related to the above measures;

function f11_ConfrontMultipleCases_2D(popol, popnick_list, initpath, outpath, legendas, iprint)
%popol useful if you have more populations for a given data set
%popnick_list cell array whose elements are the directoy names containing data produced by mainSeg_pcg2D.m
%initpath path where all the datasets are stored
%outpath path where all the analysis results will be stored
%legendas ell array whose elements are the names used to characterize the populations
%iprint if 1 print out analysis info

if ~exist(outpath, 'dir')
    mkdir(outpath)
end

% path where to save heatmaps
dirgraphs = [outpath '/grafici'];
if ~exist(dirgraphs, 'dir')
    mkdir(dirgraphs)
end

dirCSV = [outpath '/CSV/'];
if ~exist(dirCSV, 'dir')
    mkdir(dirCSV)
end

set(0,'DefaultFigureWindowStyle','docked');

% laura: add as arguments
minPcGxNucleus = 0;
maxPcGxNucleus = intmax;

clear temp;

PcGNumber_hist = cell(0);
PcGAreas_hist = cell(0);
NucleiAreas_hist = cell(0);
NucleiCircularity_hist = cell(0);
PcGNumber_pdf = cell(0);
PcGAreas_pdf = cell(0);
NucleiAreas_pdf = cell(0);
NucleiCircularity_pdf = cell(0);
PcGNumber_realdata = cell(0);
PcGAreas_realdata = cell(0);
NucleiAreas_realdata = cell(0);
NucleiCircularity_realdata = cell(0);
% NUMBER OF POLYCOMBS PER NUCLEUS
Number_per_nucleus = [];
Group_per_nucleus = [];

global n bins karaa lw al psize;

%figure(1), clf(1);
bins = 10;
karaa = double(label2rgb(1:size(popnick_list,2)))./255;
%karaa = zeros(1,size(popnick_list,2),3);
%karaa(1,1,3)=1; %blue
%karaa(1,2,1)=1; %red
%size(karaa)
%lw = 1.25; % linewidth
lw = 4.0; % linewidth
al = 0.001; % alpha for KS test
psize = 4;
w = 'A':'Z';

%title correction
%for kount=1:size(popnick_list,2) % cycle on population
%   legendas2{kount} = strrep(legendas{kount},'_','\_');
%end
fprintf('Analyzing dataset... %s n. of popol %d\n', popol,size(popnick_list,2));
for kount=1:size(popnick_list,2) % cycle on population
    %title correction
    legendas2{kount} = strrep(legendas{kount},'_','\_');
    filenameNofPcG = ['NofPcGxNuclei_' popnick_list{kount} '.csv'];
    filenameNucleiAreas = ['NucleiAreas_' popnick_list{kount} '.csv'];
    filenamePcGAreas = ['PcGAreas_' popnick_list{kount} '.csv'];
    filenameCircularity = ['NucleiCircularity_' popnick_list{kount} '.csv'];
    fprintf('Analyzing popol... %s\n', popnick_list{kount});
    %dare = [initpath popol '/' popnick_list{kount} '/STCmat/'];
    dare = [initpath popol popnick_list{kount} '/STCmat/'];
    STCmat = [popnick_list{kount} '_STC.mat'];
    load([dare STCmat]);
    FirstSTC = STC;
    clear STC;

    % NUMBER OF POLYCOMBS PER NUCLEUS
    if iprint
       fprintf('--> Number of PcG\n');
       fprintf(' -- %s N. of Series %d\n',popnick_list{kount}, size(FirstSTC,2));
    end
    TOTNUMPcG = [];
    Number_per_nucleus = [];
    Group_per_nucleus = [];
    for s=1:size(FirstSTC,2)
        if iprint
           fprintf('    Series %d\n',s);
        end
        voxelvol = FirstSTC(s).Scale.x * FirstSTC(s).Scale.y * FirstSTC(s).Scale.z;
        n = 1;
        while n <= size(FirstSTC(s).NCL,2)
             PcGNum = FirstSTC(s).PcG{n}.NumObjects;
             % code part to discriminate PcG with a certain area
             %for p=1:FirstSTC(s).PcG{n}.NumObjects
             %    parea = size(FirstSTC(s).PcG{n}.PixelIdxList{p},1)*voxelvol;
             %    if parea > 0.2
             %       PcGNum = PcGNum - 1;
             %    end
             %end
             %TOTNUMPcG = [TOTNUMPcG FirstSTC(s).PcG{n}.NumObjects];
             %Number_per_nucleus(end+1) = FirstSTC(s).PcG{n}.NumObjects;
             TOTNUMPcG = [TOTNUMPcG PcGNum];
             Number_per_nucleus(end+1) = PcGNum;
             Group_per_nucleus(end+1,:) = legendas{kount};
             if iprint
                fprintf(' -- -- Nucleus %d N. of PcG %d\n',n, FirstSTC(s).PcG{n}.NumObjects);
             end
             n = n + 1;
        end%while n
    end%for s
    % --- PLOT --- 
    [prob,pdfE] = PrintStatisticsSingle(TOTNUMPcG,karaa(1,kount,:),lw,...
        'percentage of Nuclei','','\bf Number of PcG Bodies per Nucleus',.035,10);
    clear dists randists ks pv;
    hold('on')
    if kount==size(popnick_list,2)
        hold('off');
        legend(legendas);
    end%if
    PcGNumber_hist{kount} = prob;
    PcGNumber_pdf{kount} = pdfE;
    PcGNumber_realdata{kount} = TOTNUMPcG;
    dlmwrite([dirCSV filenameNofPcG],PcGNumber_realdata{kount},'precision','%.6f','-append');

    % PcG VOLUMES
    if iprint
       fprintf('--> PcG VOLUMES \n');
       fprintf(' -- %s N. of Series %d\n',popnick_list{kount}, size(FirstSTC,2));
    end
    TOTPcGVol = [];
    PcGno = 0;
    for s=1:size(FirstSTC,2)
        if iprint
           fprintf('    Series %d\n',s);
        end
        voxelvol = FirstSTC(s).Scale.x * FirstSTC(s).Scale.y * FirstSTC(s).Scale.z;
        n = 1;
        while n <= size(FirstSTC(s).NCL,2)
%%%%%%%%%%%%
             if FirstSTC(s).PcG{n}.NumObjects > 0 & any(any(any(FirstSTC(s).NCL{n}.Nucleus)))
                 PcGno = PcGno + FirstSTC(s).PcG{n}.NumObjects;
                 for p=1:FirstSTC(s).PcG{n}.NumObjects
                     parea = size(FirstSTC(s).PcG{n}.PixelIdxList{p},1)*voxelvol;
                     %if parea > 0 && parea <= 0.2
                     if parea > 0
                        TOTPcGVol = [TOTPcGVol parea];
                        fprintf('Pop %d Series %d Nucleo %d  PcGArea %f\n',kount, s, n, parea);  
                     else 
                        fprintf('Pop %d Series %d Nucleo %d  PcGArea too big %f\n',kount, s, n, parea);  
                     end
                 end%for
             end%if
%%%%%%%%%%%
             n = n + 1;
        end%while n
    end%for s
    % --- PLOT --- 
    [prob,pdfE] = PrintStatisticsSingle(TOTPcGVol,karaa(1,kount,:),lw,...
        'percentage of PcG bodies','','\bf Area of PcG',.035,10);
    clear dists randists ks pv;
    hold('on')
    if kount==size(popnick_list,2)
        hold('off');
        legend(legendas);
    end%if
    PcGAreas_hist{kount} = prob;
    PcGAreas_pdf{kount} = pdfE;
    PcGAreas_realdata{kount} = TOTPcGVol;
    dlmwrite([dirCSV filenamePcGAreas],PcGAreas_realdata{kount},'precision','%.6f','-append');

    % NUCLEI VOLUMES
    if iprint
       fprintf('--> NUCLEI VOLUMES \n');
       fprintf(' -- %s N. of Series %d\n',popnick_list{kount}, size(FirstSTC,2));
    end
    TOTNCLVol = [];
    maxN = 0;
    for s=1:size(FirstSTC,2)
        if iprint
           fprintf('    Series %d\n',s);
        end
        voxelvol = FirstSTC(s).Scale.x * FirstSTC(s).Scale.y * FirstSTC(s).Scale.z;
        %for n=1:size(FirstSTC(s).NCL,2)
        n = 1;
        while n <= size(FirstSTC(s).NCL,2)
%%%%%%%%%%%
             NumPcGxNuclei = FirstSTC(s).PcG{n}.NumObjects;
             if NumPcGxNuclei >= minPcGxNucleus && NumPcGxNuclei <= maxPcGxNucleus
                if  (FirstSTC(s).PcG{n}.NumObjects) > maxN
                    maxN = FirstSTC(s).PcG{n}.NumObjects;
                end%if

                % Nucleus Analisys: Area, Area
                %[bwl_NCL nCC_NCL] = bwconncomp(NCL{n}.Nucleus); % nCC_NCL should be 1...
                CC_NCL  = bwconncomp(FirstSTC(s).NCL{n}.Nucleus); % nCC_NCL should be 1...
                nCC_NCL = CC_NCL.NumObjects;

                if ( nCC_NCL == 1 ) % nucleo intero
                   TOTNCLVol = [TOTNCLVol size(CC_NCL.PixelIdxList{1},1)*voxelvol];
                else               % nucleo scomposto
                   CC_NCLTotVol = 0.;
                   for i=1:nCC_NCL
                       CC_NCLTotVol = CC_NCLTotVol + size(CC_NCL.PixelIdxList{i},1)*voxelvol;
                   end
                   TOTNCLVol = [TOTNCLVol CC_NCLTotVol];
                end%if nCC_NCL
                %fprintf('Pop %d Series %d Nucleo %d Objetcts %d Area %f\n',toti, s, n, nCC_NCL, TOTNCLVol(n));
             end%if NumPcGxNuclei
%%%%%%%%%%%%
             n = n + 1;
        end%while n
    end%for s
    % --- PLOT --- 
    [prob,pdfE] = PrintStatisticsSingle(TOTNCLVol,karaa(1,kount,:),lw,...
        'percentage of Nuclei','','\bf Area of Nuclei',.035,10);
    clear dists randists ks pv;
    hold('on')
    if kount==size(popnick_list,2)
        hold('off');
        legend(legendas);
    end%if
    NucleiAreas_hist{kount} = prob;
    NucleiAreas_pdf{kount} = pdfE;
    NucleiAreas_realdata{kount} = TOTNCLVol;
    dlmwrite([dirCSV filenameNucleiAreas],NucleiAreas_realdata{kount},'precision','%.6f','-append');

    % NUCLEI ECCENTRICITY
    if iprint
       fprintf('--> NUCLEI ECCENTRICITY \n');
       fprintf(' -- %s N. of Series %d\n',popnick_list{kount}, size(FirstSTC,2));
    end
    TOTNCLCircularity = [];
    maxN = 0;
    for s=1:size(FirstSTC,2)
        if iprint
           fprintf('    Series %d\n',s);
        end
        %for n=1:size(FirstSTC(s).NCL,2)
        n = 1;
        while n <= size(FirstSTC(s).NCL,2)
%%%%%%%%%%%
             % Nucleus Analisys: Area, Area
             %[bwl_NCL nCC_NCL] = bwconncomp(NCL{n}.Nucleus); % nCC_NCL should be 1...
             CC_NCL  = bwconncomp(FirstSTC(s).NCL{n}.Nucleus); % nCC_NCL should be 1...
             nCC_NCL = CC_NCL.NumObjects;

             if ( nCC_NCL == 1 ) % nucleo intero
                TOTNCLCircularity = [TOTNCLCircularity FirstSTC(s).NCL{n}.eccentricity];
                %fprintf('Pop %d Series %d Nucleo %d  metric %f\n',kount, s, n,  FirstSTC(s).NCL{n}.eccentricity);  
             else               % nucleo scomposto
                fprintf('Pop %d Series %d Nucleo %d more than 1 Objetct %d\n',kount, s, n, nCC_NCL);  
             end%if nCC_NCL
             %fprintf('Pop %d Series %d Nucleo %d Objetcts %d Area %f\n',kount, s, n, nCC_NCL, TOTNCLVol(n));
%%%%%%%%%%%%
             n = n + 1;
        end%while n
    end%for s
    % --- PLOT --- 
    [prob,pdfE] = PrintStatisticsSingle(TOTNCLCircularity, karaa(1,kount,:),lw,...
        'percentage of Nuclei','','\bf Circularity of Nuclei',.035,10);
    clear dists randists ks pv;
    hold('on')
    if kount==size(popnick_list,2)
        hold('off');
        legend(legendas);
    end%if
    NucleiCircularity_hist{kount} = prob;
    NucleiCircularity_pdf{kount} = pdfE;
    NucleiCircularity_realdata{kount} = TOTNCLCircularity;
clear FirstSTC;
dlmwrite([dirCSV filenameCircularity],NucleiCircularity_realdata{kount},'precision','%.6f','-append')
end%for kount
%pause

figure(9), clf(9);
xgrid = linspace(0,max(TOTNUMPcG),10); grid('on');
bar(xgrid,reshape(cell2mat(PcGNumber_hist),10,size(popnick_list,2)),1);
bar(reshape(cell2mat(PcGNumber_hist),10,size(popnick_list,2)),1);
h = get(gca,'child');
for kount=1:size(popnick_list,2)
    set(h(size(popnick_list,2)-kount+1),'FaceColor',max(0,karaa(1,kount,:)./1.25),'EdgeColor',max(0,karaa(1,kount,:)./1.25),'LineWidth',1);
end%for
legend(legendas2); title('\bf Number of PcG Bodies per Nucleus'); 
%xlim([-.2,1.2]);
xlabel(''), ylabel('percentage of Nuclei');
xgrid = linspace(0,90,10); grid('on');
%axes('Position',[0 0 1 .99],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
%    'Visible','off','Units','normalized','clipping','off');
axes('Position',[0 0 1 .99],'Box','off',...
    'Visible','off','Units','normalized','clipping','off');
set(gcf,'PaperPosition',[0 0 5.5 4]); set(gcf,'PaperSize',[psize psize]);
print(gcf, '-djpeg', '-r100', [outpath popol 'Hists_PcGNumber.jpg']);
copyfile([outpath popol 'Hists_PcGNumber.jpg'], dirgraphs);

figure(10), clf(10);
xgrid = linspace(min(TOTNCLVol),max(TOTNCLVol),10); grid('on');
bar(xgrid,reshape(cell2mat(NucleiAreas_hist),10,size(popnick_list,2)),1);
h = get(gca,'child');
for kount=1:size(popnick_list,2)
    set(h(size(popnick_list,2)-kount+1),'FaceColor',max(0,karaa(1,kount,:)./1.25),'EdgeColor',max(0,karaa(1,kount,:)./1.25),'LineWidth',1);
end%for
legend(legendas2); title('\bf Nuclei Areas');
xlabel(''), ylabel('percentage of Nuclei');
xgrid = linspace(0,90,10); grid('on');
axes('Position',[0 0 1 .99],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
    'Visible','off','Units','normalized','clipping','off');
set(gcf,'PaperPosition',[0 0 5.5 4]); set(gcf,'PaperSize',[psize psize]);
print(gcf, '-djpeg', '-r100', [outpath popol 'Hists_NucleiAreas.jpg']);
copyfile([outpath popol 'Hists_NucleiAreas.jpg'], dirgraphs);

figure(11), clf(11);
%xgrid = linspace(min(TOTPcGVol),max(TOTPcGVol),10); grid('on');
%min(TOTPcGVol)
%max(TOTPcGVol)
xgrid = linspace(0,1,10); grid('on');
bar(xgrid,reshape(cell2mat(PcGAreas_hist),10,size(popnick_list,2)),1);
xlim([-0.5,1.5]);
h = get(gca,'child');
for kount=1:size(popnick_list,2)
    set(h(size(popnick_list,2)-kount+1),'FaceColor',max(0,karaa(1,kount,:)./1.25),'EdgeColor',max(0,karaa(1,kount,:)./1.25),'LineWidth',1);
end%for
legend(legendas2); title('\bf Areas of PcG');
xlabel(''), ylabel('percentage of PcG');
xgrid = linspace(0,90,10); grid('on');
axes('Position',[0 0 1 .99],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
    'Visible','off','Units','normalized','clipping','off');
set(gcf,'PaperPosition',[0 0 5.5 4]); set(gcf,'PaperSize',[psize psize]);
print(gcf, '-djpeg', '-r100', [outpath popol 'Hists_PcGAreas.jpg']);
copyfile([outpath popol 'Hists_PcGAreas.jpg'], dirgraphs);

figure(14), clf(14);
group = []; tha_data = [];
for kount=1:size(popnick_list,2)
    group = [group; repmat({legendas{kount}},size(PcGNumber_realdata{kount},2),1)];
    tha_data = [tha_data PcGNumber_realdata{kount}];
end%for
h = boxplot(tha_data, group, 'Notch', 'off', 'Orientation', 'horizontal', ...
    'OutlierSize', 3, 'Symbol', '.', 'Jitter', .2, 'Positions', size(popnick_list,2):-1:1);
for kount=1:size(h,2)
    set(h(:,size(h,2)-kount+1),'Color',karaa(1,kount,:)./1.25,'LineWidth',lw);
    bbox= findobj(gcf,'tag','Box');
    %for j=1:length(bbox)
    %   patch(get(bbox(j),'XData'),get(bbox(j),'YData'),'y','FaceColor',karaa(1,kount,:)./2,'FaceAlpha',.5);
    %end
end%for
title('\bf \fontname{SansSerif} Number of PcG Bodies per Nucleus '); grid('on');
axes('Position',[0 0 1 .99],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
    'Visible','off','Units','normalized','clipping','off');
set(gcf,'PaperPosition',[0 0 5.5 4]); set(gcf,'PaperSize',[psize psize]);
print(gcf, '-djpeg', '-r100', [outpath popol 'Boxplot_PcGNumber.jpg']);
copyfile([outpath popol 'Boxplot_PcGNumber.jpg'], dirgraphs);

figure(15), clf(15);
group = []; tha_data = [];
for kount=1:size(popnick_list,2)
    group = [group; repmat({legendas{kount}},size(NucleiAreas_realdata{kount},2),1)];
    tha_data = [tha_data NucleiAreas_realdata{kount}];
end%for
h = boxplot(tha_data, group, 'Notch', 'off', 'Orientation', 'horizontal', ...
    'OutlierSize', 3, 'Symbol', '.', 'Jitter', .2, 'Positions', size(popnick_list,2):-1:1);
for kount=1:size(h,2)
    set(h(:,size(h,2)-kount+1),'Color',karaa(1,kount,:)./1.25,'LineWidth',lw);
    bbox= findobj(gcf,'tag','Box');
    %for j=1:length(bbox)
    %   patch(get(bbox(j),'XData'),get(bbox(j),'YData'),'y','FaceColor',karaa(1,kount,:)./2,'FaceAlpha',.5);
    %end
end%for
title('\bf \fontname{SansSerif} Areas of Nuclei '); grid('on');
axes('Position',[0 0 1 .99],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
    'Visible','off','Units','normalized','clipping','off');
set(gcf,'PaperPosition',[0 0 5.5 4]); set(gcf,'PaperSize',[psize psize]);
print(gcf, '-djpeg', '-r100', [outpath popol 'Boxplot_NucleiAreas.jpg']);
copyfile([outpath popol 'Boxplot_NucleiAreas.jpg'], dirgraphs);

figure(16), clf(16);
group = []; tha_data = [];
for kount=1:size(popnick_list,2)
    group = [group; repmat({legendas{kount}},size(PcGAreas_realdata{kount},2),1)];
    tha_data = [tha_data PcGAreas_realdata{kount}];
end%for
h = boxplot(tha_data, group, 'Notch', 'off', 'Orientation', 'horizontal', ...
    'OutlierSize', 3, 'Symbol', '.', 'Jitter', .2, 'Positions', size(popnick_list,2):-1:1);
for kount=1:size(h,2)
    set(h(:,size(h,2)-kount+1),'Color',karaa(1,kount,:)./1.25,'LineWidth',lw);
    bbox= findobj(gcf,'tag','Box');
    %for j=1:length(bbox)
    %   patch(get(bbox(j),'XData'),get(bbox(j),'YData'),'y','FaceColor',karaa(1,kount,:)./2,'FaceAlpha',.5);
    %end
end%for
title('\bf \fontname{SansSerif} Areas of PcG '); grid('on');
axes('Position',[0 0 1 .99],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
    'Visible','off','Units','normalized','clipping','off');
set(gcf,'PaperPosition',[0 0 5.5 4]); set(gcf,'PaperSize',[psize psize]); 
print(gcf, '-djpeg', '-r100', [outpath popol 'Boxplot_PcGAreas.jpg']);
copyfile([outpath popol 'Boxplot_PcGAreas.jpg'], dirgraphs);

figure(17), clf(17);
group = []; tha_data = [];
for kount=1:size(popnick_list,2)
    group = [group; repmat({legendas{kount}},size(NucleiCircularity_realdata{kount},2),1)];
    tha_data = [tha_data NucleiCircularity_realdata{kount}];
end%for
h = boxplot(tha_data, group, 'Notch', 'off', 'Orientation', 'horizontal', ...
    'OutlierSize', 3, 'Symbol', '.', 'Jitter', .2, 'Positions', size(popnick_list,2):-1:1);
for kount=1:size(h,2)
    set(h(:,size(h,2)-kount+1),'Color',karaa(1,kount,:)./1.25,'LineWidth',lw);
    bbox= findobj(gcf,'tag','Box');
    %for j=1:length(bbox)
    %   patch(get(bbox(j),'XData'),get(bbox(j),'YData'),'y','FaceColor',karaa(1,kount,:)./2,'FaceAlpha',.5);
    %end
end%for
title('\bf \fontname{SansSerif} Circularity of Nuclei '); grid('on');
axes('Position',[0 0 1 .99],'Xlim',[0 1],'Ylim',[0 1],'Box','off',...
    'Visible','off','Units','normalized','clipping','off');
set(gcf,'PaperPosition',[0 0 5.5 4]); set(gcf,'PaperSize',[psize psize]);
print(gcf, '-djpeg', '-r100', [outpath popol 'Boxplot_NucleiCircularity.jpg']);
copyfile([outpath popol 'Boxplot_NucleiCircularity.jpg'], dirgraphs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

