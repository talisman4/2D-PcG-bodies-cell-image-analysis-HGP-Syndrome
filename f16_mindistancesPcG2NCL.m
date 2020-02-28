function f16_mindistancesPcG2NCL(popol, popnick_list, outCellTot, initpath, outpath, iprint)

if ~exist(outpath, 'dir')
    mkdir(outpath)
end


fprintf('Analyzing ... %s\n', popol);
for kount=1:size(popnick_list,2) % cycle on populations
    stack_format = popnick_list{kount}
    dirMinDist = [outpath '/' popol '/' stack_format '_MinDistPNG'];
    mkdir(dirMinDist);
    filenameMindistPcG2NCL = [dirMinDist '/' popnick_list{kount} '_MindistancesPcG2NCL.csv'];
    fileIDmindist = fopen(filenameMindistPcG2NCL,'w');
    %dare = [initpath popol '/' popnick_list{kount} '/STCmat/']
    dare = [initpath popol popnick_list{kount} '/STCmat/']
    STCmat = [popnick_list{kount} '_STC.mat']
    load([dare STCmat]);
    FirstSTC = STC;
    clear STC;

    outCell = outCellTot{kount};
    if iprint
       fprintf(' -- %s N. of Series %d\n',popnick_list{kount}, size(FirstSTC,2));
    end
    for s=1:size(FirstSTC,2) % ciclo sulle serie
        filenameformat = strcat(stack_format, sprintf('_Series%03d',s));
        if iprint
           fprintf(' Serie %03d Number of Nuclei %d\n', s,size(FirstSTC(s).NCL,2));
        end
        if outCell(s,1) ~= -1
           n = 1;
           pOutNCL = 1;
           while n <= size(FirstSTC(s).NCL,2)  % cycle on Nuclei
              if  n ~= outCell(s,pOutNCL)
                 fprintf('  Nucleus n. %d Number of PcG %d\n', n, FirstSTC(s).PcG{n}.NumObjects);
                 %defining nucleus contour
                 regIm2stats = regionprops(FirstSTC(s).NCL{n}.Nucleus, 'PixelList');
                 shp_NCLn = boundary(regIm2stats.PixelList);
                 fig_shape = figure(n+100); set(fig_shape, 'Visible', 'off');
                 hold on
                 %plot nuclear contour
                 plot(regIm2stats.PixelList(shp_NCLn,1),regIm2stats.PixelList(shp_NCLn,2),'b','LineWidth',2);
                 %plot nuclear centroid
                 CentroidX = FirstSTC(s).NCL{n}.NuclCentr(2);
                 CentroidY = FirstSTC(s).NCL{n}.NuclCentr(1);
                 plot(CentroidX,CentroidY,'b.','MarkerSize',20);
                 if FirstSTC(s).PcG{n}.NumObjects > 0
                    for k=1:FirstSTC(s).PcG{n}.NumObjects
                       %Compute the distance of the n-th Nuclei (shape) from the k-th PcG centroid
                       distances_shpNCLn2PcGk = sqrt((regIm2stats.PixelList(shp_NCLn,1) - FirstSTC(s).PcG{n}.Centroid{k}(2)).^2 + (regIm2stats.PixelList(shp_NCLn,2) - FirstSTC(s).PcG{n}.Centroid{k}(1)).^2);
                       %Find the closest one
                        [mindist_NCLn2PcGk, index_NCLn2PcGk] = min(distances_shpNCLn2PcGk);
                        %PcG centroid 
                        PcGX = FirstSTC(s).PcG{n}.Centroid{k}(2);
                        PcGY = FirstSTC(s).PcG{n}.Centroid{k}(1);
                        plot(PcGX,PcGY,'gd','MarkerSize',8);
                        %closest point on the Nuclear boundary to PcG
                        closestX = regIm2stats.PixelList(shp_NCLn(index_NCLn2PcGk),1);
                        closestY = regIm2stats.PixelList(shp_NCLn(index_NCLn2PcGk),2);
                        %plot distance
                        plot(closestX,closestY,'r.','MarkerSize',20);
                        line([closestX, PcGX], [closestY, PcGY],'LineWidth', 1.5, 'Color', 'r');
                        %distance between nuclear Centroid and the clostet point (on the nuclear boundary) to PcG
                        dist_NCLnCentr2Closest = sqrt((CentroidX - closestX).^2 + (CentroidY - closestY).^2);
                        ratio_mindistANDdistNCLCenter2Closest = mindist_NCLn2PcGk/dist_NCLnCentr2Closest;
                        %distance between nuclear Centroid and PcG
                        dist_NCLnCentr2PcGk = sqrt((CentroidX - PcGX).^2 + (CentroidY - PcGY).^2);
                        ratio_distPcGkandNCLCenter = dist_NCLnCentr2PcGk/dist_NCLnCentr2Closest;
                        %fprintf(fileIDmindist, 'Series%03d, Nuclei%d, %f, %f, %f, %f, %f\n', s, n, dist_NCLnCentr2Closest, dist_NCLnCentr2PcGk, ratio_distPcGkandNCLCenter, mindist_NCLn2PcGk, ratio_mindistANDdistNCLCenter2Closest);
                        fprintf(fileIDmindist, 'Series%03d, Nuclei%d, %f, %f, %f\n', s, n, dist_NCLnCentr2Closest, mindist_NCLn2PcGk, ratio_mindistANDdistNCLCenter2Closest);
                        %fprintf('nuclei %d PcG %d index_NCLn2PcGk %d (%d,%d) mindist %f\n', n, k, index_NCLn2PcGk, closestX, closestY, mindist_NCLn2PcGk);
                        fprintf('Series%03d, Nuclei%d, PcG %d, %f, %f, %f, %f, %f\n', s, n, k, dist_NCLnCentr2Closest, dist_NCLnCentr2PcGk, ratio_distPcGkandNCLCenter, mindist_NCLn2PcGk, ratio_mindistANDdistNCLCenter2Closest);
                    end%for k
                    print(gcf, '-dpng', '-r100', [dirMinDist '/' filenameformat,  '_mindist_' int2str(n) '.png']);

                 %else
                 %   mindist_NCLn2PcGk = -1.;
                 end%if
              else
                 if iprint
                    fprintf('  discarding Nucleus n. %d \n', n);
                 end
                 pOutNCL = pOutNCL + 1;
              end%if
              close(n+100)
              n = n + 1;
           end%while n
        else
           if iprint
              fprintf(' Discarded Serie %d \n',s);
           end
        end%if
    end%for s
    fclose(fileIDmindist);
end%for kount
