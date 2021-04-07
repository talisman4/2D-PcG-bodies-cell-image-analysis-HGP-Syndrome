%%% Begin Main Function
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
% main program Nuclei Region Segmentation + PcG detection through Isodata thresholding operation
% In general, first channel  is Blu   (ch00/Dapi), 
%             second channel is Green (ch01/PcG),
%             third channel  is Red   (ch02/Lamin)

%%function InitFunU(img, iNy, iNx);

%%clear all;


function mainSeg_pcg2D(population,segmentation_done,only_segmentation, ...
                       dirimages,print_thresh,dir_save,flag_seg,mu,lambda)
% population data input structure: for each population(i), required fields are:
% .stack_basename is the reference string name of the population i
% where names for corresponding images have to be of the format 'stack_basename''stack_format'
% (e.g. stack_format, reflecting filenames, might be '%03d_z%02d_ch%02d.tif' where %03d=serie's number, %02d=image plane number, %02d=channel number)
% .name population string name used for the data output directory
% .xyscale scale x-axis valus (=scale y-axis value)
% can be retrived from Voxel parameter value of DimensionDescription label (for DimID=X or DimID=Y) in the xml image metadata file
% .series vector of string numbers, each corresponding to a serie's number of input images
% .aplane vector of string numbers, each corresponding to a plane of the serie's stack to be analyzed
% .channels number of the channels of each image (at least 2: dapi and pcg)
% .thresh -255 the threshold used for the thresholding operation has to be computed by ISODATA,
%         otherwise uses a fixed threshold
% optional field is:
% .chext vector containing the channel numbers corresponding to dapi/pcg/lamin respectively 
% segmentation_done if 1 segmentation already performed
% only_segmentation if 1 performs only segmentation
% dirimages relative path where input images are stored
% print_thresh if 1 prints the threshold given in input
% dir_save if 1 store output data in the same dir as population data files
% flag_seg nuclei segmentation strategy
% mu curvature weight for active contour function
% lambda regions weight for active contour function


fprintf('----------------------------------------------------------\n');
fprintf('             Nuclei Segmentation and PcG Detection        \n');
fprintf('----------------------------------------------------------\n');

if ~exist('flag_seg','var')
   flag_seg = 0;
else
   switch flag_seg
      case 7
         fprintf('Segmentation on ImDap_seg = (medfilt2(ImPcG) + medfilt2(ImLam))/2\n');
      case 6
         fprintf('Segmentation on ImDap_seg = (medfilt2(ImDap) + ImLam)/2\n');
      case 5
         fprintf('Segmentation on ImDap_seg = (ImPcG + ImLam)/2\n');
      case 4
         fprintf('Segmentation on ImPcG + abs(medfilt2(ImLam,[4,4])\n');
      case 3
         fprintf('Segmentation on (ImDap + medfilt2(ImLam,[4,4])/2)\n');
      case 2
         fprintf('Segmentation on (ImDap + ImLam)/2\n');
      case 1
         fprintf('Segmentation on ImLam\n');
      otherwise
         fprintf('Segmentation on ImDap\n');
     end
end

if ~exist('segmentation_done','var')
   segmentation_done = 0;
end

if ~exist('only_segmentation','var')
   only_segmentation = 0;
end

if ~exist('dir_save','var')
   dir_save = 0; % set to 1 to save output data in the same dir as population data files
end

if ~exist('print_thresh','var')
   print_thresh = 0;
end

if ~exist('lambda','var')
   lambda = 0.0001;
end

if ~exist('mu','var')
   mu = 10;
end

if ~exist('population','var')
   fprintf('population record not defined\n');
   exit;
end

% integers used to indicate channels in file names
dapiext_ = 0;
greenext_ = 1;
redext_ = 2;

% scale factor
sc = 255.;

Scale = struct();
Scale.z = 1;

create_image = 1;

clean_debris = 1; debris_size = 650; % if there are spurios pixels in the segmentation

intensityPcG = 0.;
intensityPcGGreen = 0.;
intensityNCL = 0.;
intensityNCLGreen = 0.;
medintensityPcG = 0.;
medintensityPcGGreen = 0.;
medintensityNCL = 0.;
medintensityNCLGreen = 0.;

for idx = 1:size(population,2) %cycle on population
    clear STC; STC = struct; 
    STC.PcG = cell(0); STC.NCL = cell(0); STC.LumPcG = cell(0); STC.LumDap = cell(0); STC.LumLam = cell(0); 
    STC.Path = []; 
    STC.Scale = struct;

    dirpop = population(idx).name;
    if (~dir_save)
        if ~exist(dirpop, 'dir')
          mkdir(dirpop)
        end
    end
    Scale.x = population(idx).xyscale;
    Scale.y = population(idx).xyscale;
    stack_basename = population(idx).stack_basename;
    stack_format = population(idx).stack_format;
    if isfield(population(idx),'chext') && ~isempty(population(idx).chext)
      if population(idx).channels==3
        fprintf('channel numbers provided dap=%d pcg=%d lam=%d\n', ...
                population(idx).chext(1), population(idx).chext(2), ...
                population(idx).chext(3));
      else
        fprintf('channel numbers provided dap=%d pcg=%d\n', ...
                population(idx).chext(1), population(idx).chext(2));
      end
      dapiext = population(idx).chext(1);
      greenext = population(idx).chext(2);
      redext = population(idx).chext(3);
    else
      fprintf('default channel numbers\n');
      dapiext = dapiext_;
      greenext = greenext_;
      redext = redext_;
    end

    fprintf('----------------------------------------------------------\n');

    K   = 1000;      % Set max number of outer loop iteratioms (default 30)
    GS   = 50;     % Set max number of inner loop iterations (default 50)
    T0  = 0;       % Set to use a preprocessing with k x k weighted average filter for Dapi frame (default do not preproc)
    T1  = 0;       % Set to use a preprocessing with k x k weighted average filter for PcG  frame (default do not preproc)
    T2  = 0;       % Set to use a preprocessing with k x k weighted average filter for Lam  frame (default do not preproc)
    I   = 7;       % Set a k x k weighted average filter before applying IsoData (default use a 4 x 4 filter)

    fprintf('Segmentation Parameter: lambda %f mu %f it_out %d it_in %d I %d \n', lambda, mu, K, GS, I);

    dirSeparatedNuclei = [population(idx).name '/SeparatedNuclei'];
    if dir_save 
       dirSeparatedNuclei = [dirimages dirSeparatedNuclei];
    end 
    if ~exist(dirSeparatedNuclei, 'dir')
       mkdir(dirSeparatedNuclei)
    end
    dirCSV = [population(idx).name '/CSV'];
    if dir_save 
       dirCSV = [dirimages dirCSV];
    end 
    if ~exist(dirCSV, 'dir')
       mkdir(dirCSV)
    end

    if (print_thresh)
       filenameThresh = [ dirCSV '/ThreshPerSeries_' population(idx).name '.csv'];
       fileIDt = fopen(filenameThresh,'w');
    end
    for stack = 1:size(population(idx).series,2) % cycle on Series
        id_stack  = population(idx).series(stack);
        % filenameformat = [stack_basename '%03d_z%02d_ch%02d.tif']; % series,plane,channel
        filenameformatout = [stack_basename '%03d_z%02d'];       % series,plane
        filenameformat = [stack_basename stack_format]; % series,plane,channel
    
        fprintf('Population -- %s N. of Stacks %d\n', population(idx).name, size(population(idx).series,2));

        % set first and last plane of the stack
	pf = population(idx).aplane(stack);
        pl = pf; %2D->last plane is pf
	% retrieve thresh
	thresh = population(idx).thresh;

        fprintf(' --> Stack/Series%03d of 1 plane \n', id_stack);

        dirSeries = strcat('Series', sprintf('%03d',id_stack));
        dirtif =[population(idx).name '/' dirSeries '/tif'];
        if dir_save 
           dirtif = [dirimages dirtif];
        end 
        if ~exist(dirtif, 'dir')
          mkdir(dirtif)
        end
        dirSTC =[population(idx).name '/STCmat'];
        if dir_save 
           dirSTC = [dirimages dirSTC];
        end 
        if ~exist(dirSTC, 'dir')
          mkdir(dirSTC)
        end
        filenameNofNuclei = [ dirCSV '/NucleiPerSeries_' population(idx).name '.csv'];

        % acquire image dimensions
        DapFilename = sprintf( filenameformat, id_stack, pf, dapiext );
        %dirReadImages = [dirimages population(idx).name '/' dirSeries '/'];
        dirReadImages = [dirimages];
        fieldsize = size(imread( [dirReadImages DapFilename] ));
        % initialize vectors defining nuclei regions and PcG bodies
        height = fieldsize(1); width = fieldsize(2); 
        if (~only_segmentation) %Francesco 20/05/2016
            ImPcG_vol = zeros(height,width,1);
            regIm2_vol = false(height,width,1);
        end
   
        if (print_thresh)
           fprintf(fileIDt, 'Series%03d ,', population(idx).series(stack));
	   if (thresh ~= -255)
               fprintf(fileIDt, '%3.8f\n', population(idx).thresh);
	   end
        end
        for plane = pf:pl %cycle on planes

            fprintf('     Images N. %d Series%03d\n', plane, id_stack);

            % reading blu(dapi) image
            DapFilename = sprintf( filenameformat, id_stack, plane, dapiext );
	    ImTemp = imread([dirReadImages DapFilename]);
            if size(ImTemp,3) == 3
              ImDap = double(ImTemp(:,:,3)); 
            else
              ImDap = double(ImTemp);
            end
            %imagesc(ImDap);colormap(gray);title('Dapi Image');
            %imshow(ImDap);
            %pause
            % scaling image pixels between 0 and 255
            maxDap = max(max(ImDap));
            minDap = min(min(ImDap));
            ImDap_sc = (ImDap - minDap)/(maxDap - minDap) * sc;
            %maxDap = max(max(ImDap));
            %minDap = min(min(ImDap));
            %fprintf('dapi -- after the scaling max=%f min=%f\n',maxDap,minDap);
            % saving blu(dapi) image in gray scale
            %m8 = uint8(ImDap);
            %imwrite(m8, DapFilename, 'tiff');
            %movefile(DapFilename,dirgrey)

            % reading green (PcG) image
            PcGFilename = sprintf( filenameformat, id_stack, plane, greenext );
            ImTemp = imread([dirReadImages PcGFilename]);
            if size(ImTemp,3) == 3
              ImPcG = double(ImTemp(:,:,2));
            else
              ImPcG = double(ImTemp);
            end
            maxPcG = max(max(ImPcG));
            minPcG = min(min(ImPcG));
            ImPcG_sc = (ImPcG - minPcG)/(maxPcG - minPcG) * sc;
            %maxGreen = max(max(ImPcG));
            %minGreen = min(min(ImPcG));
            %fprintf('     PcG -- max=%f min=%f\n',maxPcG,minPcG);
            % saving green image in gray scale
            %m8 = uint8(ImPcG);
            %imwrite(m8, PcGFilename, 'tiff');
            %movefile(PcGFilename, dirgrey)
            %figure
            %imagesc(ImPcG);colormap(gray);title('PcG Image');

            % reading red (Lamin) image
            if population(idx).channels==3
               LamFilename = sprintf( filenameformat, id_stack, plane, redext );
               ImTemp = imread([dirReadImages LamFilename]);
               if size(ImTemp,3) == 3
                 ImLam = double(ImTemp(:,:,1));
               else
                 ImLam = double(ImTemp);
               end
               % scaling image pixels between 0 and 255
               %maxLam = max(max(ImLam));
               %minLam = min(min(ImLam));
               %figure
               %imagesc(ImLam);colormap(gray);title('LamPcG Image');
               %ImRed = (ImRed - minRed)/(maxRed - minRed) * sc;
               % saving red image in gray scale
               %m8 = uint8(ImLam);
               %imwrite(m8, LamFilename, 'tiff');
               %movefile(LamFilename, dirgrey)
            else
               ImLam = ImDap;
            end

            if(~segmentation_done)
                % size of images
                [Ny,Nx] = size(ImDap);
    
                % initialization of minimizer
                pfu0 = InitFunU(ImDap, Ny, Nx);
    
                % preprocessing of dapi image by means of mean filter 
                if ( T0 > 0 )
                   w1 = fspecial('average', [T0,T0]);
                   ImDap = imfilter(ImDap,w1); 
                end
    
                % preprocessing of c01 image by means of mean filter 
                if ( T1 > 0 )
                   w1 = fspecial('average', [T1,T1]);
                   ImPcG = imfilter(ImPcG,w1); 
                end
    
                % preprocessing of c02 image by means of mean filter 
                if ( T2 > 0 )
                   w1 = fspecial('average', [T2,T2]);
                   ImLam = imfilter(ImLam,w1); 
                end
    
                % setting vector of segmentation parameters
                SegParameters = [plane; Ny; Nx; lambda; mu; K; GS]; % parameters for active contour function
    
                % region segmentation function 
                switch flag_seg
                   case 7
                      [rows,cols] = find((medfilt2(ImLam) - ImDap)>=0);
                      ImDap_seg = zeros(size(ImDap,1),size(ImDap,2));
                      ImDap_seg([rows,cols]) = (ImDap([rows,cols]) + medfilt2(ImLam([rows,cols])))/2;
                      ImDap_seg = ImDap_seg + ImLam;
                   case 6
                      ImDap_seg = (medfilt2(ImDap,[4,4]) + medfilt2(ImLam))/2;
                   case 5
                      ImDap_seg = (ImPcG + ImLam)/2;
                   case 4
                      ImDap_seg = (ImPcG + medfilt2(ImLam,[4,4]))/2;
                   case 3
                      ImDap_seg = (ImDap + medfilt2(ImLam,[4,4]))/2;
                   case 2
                      ImDap_seg = (ImDap + ImLam)/2;
                   case 1
                      ImDap_seg = ImLam;
                   otherwise
                      ImDap_seg = ImDap;
                end

                % ac_mex output list:
                % regIm highlights foreground objects
                % cntIm shows object contours
                % pfu is the segmented image
                % fMeanIn is the mean fluorescence intensity value of objects
                % fMeanOut is the mean fluorescence intensity value of background
                % den1 is the number of object pixels
                % den2 is the number of background pixels
                [regIm, cntIm, pfu, fMeanIn, fMeanOut, den1, den2] = ac_mex(ImDap_seg, ImPcG, ImLam, pfu0, double(SegParameters)); 
    
                %figure
                %imagesc(regIm);colormap(gray);title('before Fill Holes in Nuclei');
		% invert colors in regIm
                regIm2 = 255 - regIm; % regIm2 has black background and white nuclei regions - regIm is the opposite
                % fill holes in dapi regions of segmented image 
                regIm2 = imfill(regIm2,'holes');
                
                %figure
                %imagesc(regIm2);colormap(gray);title('Fill Holes In Nucleus Regions');
                %pause
    
                % names of output image files: segmented regions and contours
                tempName = sprintf( filenameformatout, id_stack, plane );
                nameReg  = ['rgn_' tempName '.tif'];
                nameCnt  = ['cnt_' tempName '.tif'];
    
                % saving segmented image in tiff format
                imwrite((255-regIm2),nameReg,'tiff'); 
                movefile(nameReg,dirtif)
    
                % defining contours of dapi(blu) image 
                temp_diff = diff(regIm2,1,1);
                dPhiX = [temp_diff; zeros(1,width)];
                temp_diff = diff(regIm2,1,2);
                dPhiY = [temp_diff, zeros(height,1)];
                imageTemp = ImDap_sc;
                size(imageTemp);
                imageTemp( (dPhiX.*dPhiX + dPhiY.*dPhiY) >= eps()) = 255;      
    
                % saving dapi(blu) image contours 
                m8 = uint8(imageTemp);
                imwrite(m8,nameCnt,'tiff');
                movefile(nameCnt,dirtif)
    
                % Filtering PcG image in order to better enhance PcG areas (before applying ISODATA)
                if ( I > 0 )
                   h = 2*I+1;
                   fprintf('     Applying Mean Filter (size %d) for IsoData \n',h);
                   w1 = fspecial('average', [h,h]);
    
                   ImPcGFilt = imfilter(ImPcG,w1);
                   imPrimePcG = ImPcG - ImPcGFilt;
                else
                   imPrimePcG = ImPcG;
                end
    
    %-----------
    
    % Applying ISODATA on PcG image
    
                % excluding background for PcG detection
                imPrimePcG(regIm2 == 0) = -255; 
    
                if (thresh == -255)
                   %fprintf('Initial thresh value == > %f\n',fMeanIn(2));
                   thresh = IsoDataThresh(imPrimePcG, Ny, Nx, fMeanIn, fMeanOut, thresh );
                   fprintf('     Thresh evaluated for PcG Image ==> %f\n',thresh);
                   fprintf(fileIDt, '%3.8f\n', thresh);
                else
		   thresh = population(idx).thresh; % Fixed thresh for ISODATA PcG filtering
                   fprintf('     Applying Fixed thresh on PcG Image == > %f\n',thresh);
                end%if thresh
                regIm = zeros(size(regIm2));
                regIm( regIm2 == 255 & imPrimePcG >= thresh ) = double(255);
                regIm( regIm2 == 255 & imPrimePcG  < thresh ) = double(100);
                % Green intensity of PcG region values in imPrime
                index_intensityPcG = find(regIm == 255);
                regIm_intensityPcG = zeros(size(regIm2));
                if length(index_intensityPcG) > 0
                   regIm_intensityPcG(index_intensityPcG) = imPrimePcG(index_intensityPcG); 
                   temp = sum(sum(regIm_intensityPcG));
                   intensityPcG = intensityPcG + temp;
                   medintensityPcG = medintensityPcG + temp/(length(index_intensityPcG));
                   %fprintf('length(index_intensityPcG) %d \n', length(index_intensityPcG));
                   %pause
                end
                % Green intensity of PcG region values in ImPcG
                regIm_intensityPcG = zeros(size(regIm2));
                if length(index_intensityPcG) > 0
                   regIm_intensityPcG(index_intensityPcG) = ImPcG(index_intensityPcG); 
                   temp = sum(sum(regIm_intensityPcG));
                   intensityPcGGreen = intensityPcGGreen + temp;
                   medintensityPcGGreen = medintensityPcGGreen + temp/(length(index_intensityPcG));
                   %pause
                end
                % Green intensity of NCL region values in imPrime
                index_intensityNCL = find(regIm >= 100);
                regIm_intensityNCL = zeros(size(regIm2));
                if length(index_intensityNCL) > 0
                   regIm_intensityNCL(index_intensityNCL) = imPrimePcG(index_intensityNCL); 
                   temp = sum(sum(regIm_intensityNCL));
                   intensityNCL = intensityNCL + temp;
                   medintensityNCL = medintensityNCL + temp/(length(index_intensityNCL));
                end
                % Green intensity of NCL region values in ImPcG
                regIm_intensityNCL = zeros(size(regIm2));
                if length(index_intensityNCL) > 0
                   regIm_intensityNCL(index_intensityNCL) = ImPcG(index_intensityNCL); 
                   temp = sum(sum(regIm_intensityNCL));
                   intensityNCLGreen = intensityNCLGreen + temp;
                   medintensityNCLGreen = medintensityNCLGreen + temp/(length(index_intensityNCL));
                end
                %fprintf('In imPrime: intensityPcG %f medintensityPcG %f intensityNCL %f medintensityNCL %f \n',intensityPcG, medintensityPcG, intensityNCL, medintensityNCL);
                %fprintf('In ImPcG: intensityPcG %f medintensityPcG %f intensityNCL %f medintensityNCL %f \n',intensityPcGGreen, medintensityPcGGreen, intensityNCLGreen, medintensityNCLGreen);
    
                % FileName of the maps of PcG images (Green)
                namePcG = ['pcg_' tempName '.tif'];
    
                % saving map of Green spots
                m8 = uint8(regIm);
                imwrite(m8,namePcG,'tiff');
                movefile(namePcG,dirtif)
            end%if(~segmentation_done)
            

            if (segmentation_done)
                PcGSegFilename = sprintf( filenameformatout, id_stack, plane );
                namePcG = ['pcg_' PcGSegFilename '.tif'];
	        regIm = imread([dirtif '/' namePcG]);
                % reading green (PcG) image  %% to retrieve original intensity values 20/11/2018
                PcGFilename = sprintf( filenameformat, id_stack, plane, greenext );
                ImTemp = imread([dirReadImages PcGFilename]);
                if size(ImTemp,3) == 3
                  ImPcG = double(ImTemp(:,:,2));
                else
                  ImPcG = double(ImTemp);
                end
            end %Francesco 20/05/2016
            ImSegPcG = regIm>100; % PcG Regions

            if (~only_segmentation) %Francesco 20/05/2016
                %figure, imshow(regIm2)
                p3d = plane + 1;
                regIm2_vol(:,:,p3d) = regIm;
                regIm2_vol(:,1,p3d) = regIm2_vol(:,2,p3d); regIm2_vol(:,end,p3d) = regIm2_vol(:,end-1,p3d); % adjust borders
                regIm2_vol(1,:,p3d) = regIm2_vol(2,:,p3d); regIm2_vol(end,:,p3d) = regIm2_vol(end-1,:,p3d); % adjust borders
                % remove nucleus regions conncted to image borders
                regIm2_vol(:,:,p3d) = imclearborder(regIm2_vol(:,:,p3d),8);
    
                ImDap_vol(:,:,p3d) = ImDap; % volume of fluorescence channel ch00
                ImPcG_vol(:,:,p3d) = ImPcG; % volume of fluorescence channel ch01
                ImLam_vol(:,:,p3d) = ImLam; % volume of fluorescence channel ch02
    
                ImSegPcG_vol(:,:,p3d)=ImSegPcG; 
            end
            fprintf('----------------------------------------------------------\n');
        end%for plane
        NofPlanes = pl - pf +1;
        medintensityPcG = medintensityPcG/NofPlanes;
        medintensityPcGGreen = medintensityPcGGreen/NofPlanes;
        medintensityNCL = medintensityNCL/NofPlanes;
        medintensityNCLGreen = medintensityNCLGreen/NofPlanes;
        fprintf('In imPrime avg (on N. of stack %d): intensityPcG %f medintensityPcG %f intensityNCL %f medintensityNCL %f \n',stack, intensityPcG, medintensityPcG, intensityNCL, medintensityNCL);
        fprintf('In ImPcG avg (on N. of stack %d): intensityPcG %f medintensityPcG %f intensityNCL %f medintensityNCL %f \n',stack, intensityPcGGreen, medintensityPcGGreen, intensityNCLGreen, medintensityNCLGreen);
        

        if(~only_segmentation)
            fieldsize = size(regIm2_vol); % 2 dimensions

            LumPcG = struct;
            LumDap = struct;
            LumLam = struct;
            LumDap.Sum = reshape(sum(sum(ImDap_vol,1),2),size(ImDap_vol,3),1);
            LumDap.Mean = reshape(mean(mean(ImDap_vol,1),2),size(ImDap_vol,3),1);
            LumPcG.Sum = reshape(sum(sum(ImPcG_vol,1),2),size(ImPcG_vol,3),1);
            LumPcG.Mean = reshape(mean(mean(ImPcG_vol,1),2),size(ImPcG_vol,3),1);
            LumLam.Sum = reshape(sum(sum(ImLam_vol,1),2),size(ImLam_vol,3),1);
            LumLam.Mean = reshape(mean(mean(ImLam_vol,1),2),size(ImLam_vol,3),1);
            LumDap.Contrast = [];
            LumPcG.Contrast = [];
            LumLam.Contrast = [];
            for l=1:size(LumDap.Sum,1)
              LumDap.Contrast(1,l) = sum( sum( ( ImDap_vol(:,:,l)-LumDap.Mean(l) ).^2,1 ),2 ) ./ (fieldsize(1)*fieldsize(2));
            end%for l
            for l=1:size(LumPcG.Sum,1)
              LumPcG.Contrast(1,l) = sum( sum( ( ImPcG_vol(:,:,l)-LumPcG.Mean(l) ).^2,1 ),2 ) ./ (fieldsize(1)*fieldsize(2));
            end%for l
            for l=1:size(LumLam.Sum,1)
              LumLam.Contrast(1,l) = sum( sum( ( ImLam_vol(:,:,l)-LumLam.Mean(l) ).^2,1 ),2 ) ./ (fieldsize(1)*fieldsize(2));
            end%for l
       
            % remove small objects from nuclei image
            if clean_debris
              regIm2_vol = bwareaopen(regIm2_vol, debris_size, 6);
            end%if
  
            fprintf('Defining Nuclei image for stack %s\n',dirSeries);
	    % convert nuclei image in binary form
            PSI = boolean(regIm2_vol);
            
	    % reconstructions of nuclei 
            nuclei_CC = bwconncomp(PSI);

	    % label nuclei
            nuclei_L = labelmatrix(nuclei_CC);

	    % compute area for each nucleus object
            nuclei_dim = [];
            for k = 1:nuclei_CC.NumObjects
              nuclei_dim(end+1) = size(nuclei_CC.PixelIdxList{k},1);
            end%fo%r
            %nuclei_dim

            % exclude too small nucleus objects
            mean_dim = mean(nuclei_dim);
            idx_PSI = find(nuclei_dim > .1*mean_dim);
            clear nuclei_dim mean_dim;
            mask_PSI = ismember(nuclei_L, idx_PSI);
            nuclei_L(mask_PSI == 0) = 0;
    
            %% Number of cells (nuclei)
            M = size(idx_PSI,2);
            disp(['Identifying cells: ' int2str(M) ' cells found']);

            % .........................................................................
    
            % PSI2(:,:,:,k), with k from 1 to M, is the nucleus
            % corresponding to the k-th label
            PSI_new = false([fieldsize,M]);
    
            for k=1:M
                PSI_new(:,:,k) = (nuclei_L==idx_PSI(k));
            end
    
            PSI = boolean(PSI_new);
            clear PSI_new regIm2_vol L nuclei_L idx_PSI;
    
            fieldsize = size(PSI);
            disp('Extracting Nuclei...');
     
            % NCL holds pixel data regarding each nucleus
            NCL = cell(0);
    
            made = fieldsize(3);
            fprintf('M= %d\n',made);

            % cycle on nuclei
            for n=1:made
                [i] = find(PSI(:,:,n));
                [x,y] = ind2sub(fieldsize(1:2),i);
                NCL{n} = struct;
                % n-th nucleus
                NCL{n}.Nucleus = PSI( min(x):max(x), min(y):max(y), n );
                NCL{n}.PcG = ImSegPcG_vol( min(x):max(x), min(y):max(y) ).*NCL{n}.Nucleus;
		% discard too small detected PcG objects which are probably just noise.
                NCL{n}.PcG = bwareaopen(NCL{n}.PcG,3,8);

                % if there is more than one nucleus object remove all but the biggest
                temp = bwconncomp(NCL{n}.Nucleus,6); % 6, 18, 26
                if temp.NumObjects>1
                    thamax = 0; thamaxk = 0;

                    for k=1:temp.NumObjects
                        if thamax<size(temp.PixelIdxList{k},1)
                            thamax = size(temp.PixelIdxList{k});
                            thamaxk = k;
                        end%if
                    end%for k

                    for k=1:temp.NumObjects
                        if k~=thamaxk
                            NCL{n}.Nucleus(temp.PixelIdxList{k}) = 0;
                        end%if
                    end%for k

                    NCL{n}.PcG = NCL{n}.PcG.*NCL{n}.Nucleus; % remove spot outside of remaining obj
                end%if

                %if ~isempty(NCL{n}.Nucleus)
                %    void = zeros(size(NCL{n}.Nucleus(:,:,1)));
                %    for q=1:size(NCL{n}.Nucleus,3)
                %        if  sum(sum(NCL{n}.Nucleus(:,:,q)))<=16
                %             NCL{n}.Nucleus(:,:,q) = void;
                %        end%if
                %    end%for q
                %    clear void;
                %end%if isempty

                % compute eccentricity of nucleus NCL{n} to provide an estimate of its "roundness"
                regIm2 = PSI(:,:,n);
                [LregIm2, nCCregIm2] = bwlabel(regIm2);
                regIm2stats = regionprops(LregIm2, 'Area', 'Perimeter', 'PixelList');
 
                %fprintf('Nucleus %d nCCregIm2 %d \n',n, nCCregIm2);

                MIN_AREA = 50;
                MIN_CIRCLE_METRIC = 0.75; %0.75

                CCregIm2area = [];

                for i=1:nCCregIm2
                  CCregIm2area(end+1) = size(regIm2stats(i).PixelList,1);
                end
    
                MIN_AREA = .5*mean(CCregIm2area);
    
                if ( nCCregIm2 == 1 )
                    if(regIm2stats(1).Area >= MIN_AREA)
                        metric = (4*pi*regIm2stats(1).Area)/(regIm2stats(1).Perimeter^2);
                        NCL{n}.eccentricity = metric; % 0 is circle 1 is parabola
                    end
                else
                    NCL{n}.eccentricity = 10.; %
                end
                fprintf('Nucleus %d area %f perimeter %f metric %f \n',n, ...
                         regIm2stats(i).Area, regIm2stats(i).Perimeter, NCL{n}.eccentricity);
            end%for n

            if create_image

                %fig_numbers = figure('Visible', 'off'); imshow(mat2gray(ImDap_vol(:,:,1)));
                fig_numbers = figure('Visible', 'off'); imshow(mat2gray(ImPcG_vol(:,:,1)));
                set(fig_numbers,'Paperposition',[0,0,4,4]);
                %set(fig_numbers, 'Position', [0 0 4 4]);
                for h=1:size(NCL,2)
                    [ti,tj] = find(PSI(:,:,h),h);
                    if (size(ti,1) == 0) || (size(tj,1) == 0)
                      for PSI_plane=1:size(PSI,3)
                        [ti,tj] = find(PSI(:,:,h));
                        if (size(ti,1) ~= 0) && (size(tj,1) ~= 0)
                          break
                        end 
                      end
                    end
                    text(tj(floor(size(tj,1))),ti(floor(size(ti,1))),int2str(h),'FontWeight','bold','Color','g','FontSize',13);
                end%for
                filenameformatnumbers = strcat(stack_basename, sprintf('%03d',id_stack));
                print(fig_numbers, '-dpng', '-r200', [ dirSeparatedNuclei '/' filenameformatnumbers '_numbers.png']);
                clear ti tj;
                close(fig_numbers);
            end%if

            clear temp thamax thamaxk PSI ImPcG_vol ImPcG_vol ImLam_vol ImDap_vol;

            C = cell(0); PcG = cell(0);

            fileID = fopen(filenameNofNuclei,'a');
            fprintf(fileID, 'Series%03d,', id_stack);
            for n=1:size(NCL,2)
                fprintf('2d reconstruction Nucleus...  %d\n',n);
                fig1 = figure(1); set(fig1, 'Visible', 'on'); clf(1);
                %set(fig1,'Position', [0 0 4 4]);
                set(fig1,'Paperposition',[0,0,4,4]);
                set(fig1,'DoubleBuffer','off');
                % PLOT CENTROID
                [Cx,Cy] = ind2sub(size(NCL{n}.Nucleus),find(NCL{n}.Nucleus));
                C{n} = [mean(Cx),mean(Cy)];
                NCL{n}.NuclCentr = C{n};
                % PLOT PROTEINS
                PcG{n} = bwconncomp(NCL{n}.PcG,6);
                %PcG{n} = bwconncomp(NCL{n}.PcG,4);
                fprintf('Nucleo n. %d PcG n. %d\n', n, PcG{n}.NumObjects);
                proimage = zeros(size(NCL{n}.Nucleus));
		proimage(NCL{n}.Nucleus) = 100;
                if PcG{n}.NumObjects > 0
                    karaa = double(label2rgb(1:PcG{n}.NumObjects,'jet',[0,0,0],'shuffle'))./255;
                    for k=1:PcG{n}.NumObjects
                        protemp = zeros(size(NCL{n}.Nucleus));
                        protemp(PcG{n}.PixelIdxList{k}) = 1;
                        [x,y] = ind2sub(size(NCL{n}.Nucleus),PcG{n}.PixelIdxList{k});
                        PcG{n}.PixelSubList{k} = [x,y];
                        PcG{n}.Centroid{k} = [mean(x),mean(y)];
                        PcG{n}.CentDist{k} = sqrt(sum((C{n} - PcG{n}.Centroid{k}).^2));
                        proimage = proimage + protemp.*255;
                    end%for k
                    imtemp = uint8(mat2gray(proimage).*255);
                else 
                    imtemp = uint8(mat2gray(proimage).*100);
                end%if PcG{n}.NumObjects
                imwrite(imtemp, [dirSeparatedNuclei '/' filenameformatnumbers,  '_' int2str(n) '.png']);
                
                fprintf(fileID, '%d,', n);
            end%for n

            close all;
            STC(stack).Path = [dirSTC ]; % dirSTC
            STC(stack).Scale = Scale;
            STC(stack).PcG = PcG;
            STC(stack).NCL = NCL;
            STC(stack).LumPcG = LumPcG;
            STC(stack).LumDap = LumDap;
            STC(stack).LumLam = LumLam;
            clear PcG NCL LumPcG LumDap LumLam;
            fprintf(fileID, '\n');
            fclose(fileID);
        end%if(~only_segmentation)
    end%for stack
    if (print_thresh)
       fclose(fileIDt);
    end
    fprintf('Saving STC ...\n');
    save([dirSTC '/' population(idx).name '_STC'],'STC');
    NofStacks = size(population(idx).series,2);
    medintensityPcG = medintensityPcG/NofStacks;
    medintensityPcGGreen = medintensityPcGGreen/NofStacks;
    medintensityNCL = medintensityNCL/NofStacks;
    medintensityNCLGreen = medintensityNCLGreen/NofStacks;
    fprintf('in ImPrime on pop %s: TOT intensityPcG %f medintensityPcG %f intensityNCL %f medintensityNCL %f \n',population(idx).name, intensityPcG, medintensityPcG, intensityNCL, medintensityNCL);
    fprintf('in ImPcG on pop %s: TOT intensityPcG %f medintensityPcG %f intensityNCL %f medintensityNCL %f \n',population(idx).name, intensityPcGGreen, medintensityPcGGreen, intensityNCLGreen, medintensityNCLGreen);
end%for idx population

fprintf('---------------------------END----------------------------\n');
fprintf('----------------------------------------------------------\n');

% clf;

%%% End of Main Function
