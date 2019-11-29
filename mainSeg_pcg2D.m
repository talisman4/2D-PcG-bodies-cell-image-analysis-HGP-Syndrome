%%% Begin Main Function
%
% calling program of Region Segmentation and Isodata
% immagini sui PcG cell-images
% le immagini sono a colori 
% Generalmente il primo campo e' Blu   (ch00/Dapi), 
%              il secondo     e' Green (ch01/PcG),
%              il terzo       e' Red   (ch02/Lamin)

%%function InitFunU(img, iNy, iNx);

%%clear all;


function mainSeg_pcg2D(population,segmentation_done,only_segmentation,dirimages,stack_format_list,only_thresh,random_case, ...
                     xyscale,zscale,dirold_population_name,dir_save, flag_seg, mu, lambda, lamin)

fprintf('----------------------------------------------------------\n');
fprintf('             Nuclei Segmentation and PcG Detection        \n');
fprintf('----------------------------------------------------------\n');

% Flag che indica la presenza della Lamina come terza immagine
%lamin = 1;

if ~exist('flag_seg','var')
   flag_seg = 0;
else
   switch flag_seg
      case 7
         fprintf('Segmentation on ImDap_seg = (medfilt2(ImPcG) + medfilt2(ImLam))/2\n');
      case 6
         fprintf('Segmentation on ImDap_seg = (medfilt2(ImDap) + ImLam)/2\n');
      case 5
         fprintf('Segmentation on  ImDap_seg = (ImPcG + ImLam)/2\n');
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
only_segmentation = 1;
end

if ~exist('dir_save','var')
dir_save = 0; % flag da attivare per salvare i dati di output sotto la stessa dir
              % della popolazione da esaminare, altrimenti li salva in dirimages
end

if dir_save
   if ~exist('dirold_population_name','var')
       fprintf('dirold_population_name not defined\n');
       exit;
   end
end

if ~exist('only_thresh','var')
only_thresh = 0;
end

if ~exist('random_case','var')
random_case = 0;
end

if ~exist('xyscale','var')
xyscale = 0.12; 
end

if ~exist('zscale','var')
   zscale = 0.21;
end

if ~exist('lambda','var')
   lambda = 0.0001;
end

if ~exist('mu','var')
   mu = 10;
end

if ~exist('population','var')
% parametri della popolazione
fprintf('population record not defined\n');
exit;
end

% Interi che indicano il colore nei nomi dei file
dapiext = 0;
greenext = 1;
redext = 2;

% scale factor
sc = 255.;

Scale = struct();
Scale.x = xyscale;
Scale.y = xyscale;
Scale.z = zscale;

create_video = 0;
create_image = 1;
create_map = 0;

if abs(xyscale - zscale) < 10
    resize_to_cube = 0; % make voxels cubes? (based on xyscale and zscale)
else
    resize_to_cube = 1; % make voxels cubes? (based on xyscale and zscale)
end%if
clean_debris = 1; debris_size = 650; % if there are spurios pixels in the segmentation

old_spherize_labeling = 0;

intensityPcG = 0.;
intensityPcGGreen = 0.;
intensityNCL = 0.;
intensityNCLGreen = 0.;
medintensityPcG = 0.;
medintensityPcGGreen = 0.;
intensityNCL = 0.;
intensityNCLGreen = 0.;
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
    stack_format = stack_format_list{idx};
    %filenameformatout = '%s%03d_z%03d'; % namestack,stack,plane

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

    dirCSV = [population(idx).name '/CSV'];
    if dir_save 
       dirCSV = [dirimages dirCSV];
    end 
    if ~exist(dirCSV, 'dir')
       mkdir(dirCSV)
    end
    if (~only_thresh && ~segmentation_done)
       filenameThresh = [ dirCSV '/ThreshPerSeries_' population(idx).name '.csv'];
       fileIDt = fopen(filenameThresh,'w');
    end
    for stack = 1:size(population(idx).series,2) % cycle on Series
        id_stack  = population(idx).series(stack);
        id_planes = population(idx).planes(stack);
        if (id_planes > 100)
           filenameformat = [stack_format '%03d_z%03d_ch%02d.tif']; % series,plane,channel
           filenameformatout = [stack_format '%03d_z%03d'];       % series,plane
        else
           filenameformat = [stack_format '%03d_z%02d_ch%02d.tif']; % series,plane,channel
           filenameformatout = [stack_format '%03d_z%02d'];       % series,plane
        end
    
        fprintf('Population -- %s N. of Stacks %d\n', population(idx).name, size(population(idx).series,2));

        if (only_thresh == 1)
            pf = population(idx).threshplanes(stack);
            pl = pf;
            thresh = -255; % Initial thresh for ISODATA application before PcG filtering
            fprintf(' Evaluating thresh... %d\n',pf);
        else
            if isfield(population(idx),'first_plane') == 0
              pf = 0;
            else
              pf = population(idx).first_plane(stack);
              fprintf(' First plane... %d\n',pf);
            end
            if isfield(population(idx),'last_plane') == 0
              pl = id_planes-1;
            else
              pl = population(idx).last_plane(stack);
              fprintf(' Last plane... %d\n',pl);
            end
            thresh = population(idx).thresh; % Fixed thresh for ISODATA PcG filtering
            if (only_segmentation)
               fprintf(' Applying fixed thresh... \n');
            end
        end

        fprintf(' --> Stack/Series%03d of %d planes \n', id_stack, population(idx).planes(stack));

        dirSeries = strcat('Series', sprintf('%03d',id_stack));
        dirtif =[population(idx).name '/' dirSeries '/tif'];
        if dir_save 
           dirtif = [dirimages dirtif];
        end 
        if ~exist(dirtif, 'dir')
          mkdir(dirtif)
        end
        % eliminare dirrec3d
        %dirrec3d =[population(idx).name '/rec3d'];
        %if dir_save 
        %   dirrec3d = [dirimages dirrec3d];
        %end 
        %if ~exist(dirrec3d, 'dir')
        %  mkdir(dirrec3d)
        %end
        dirSTC =[population(idx).name '/STCmat'];
        if dir_save 
           dirSTC = [dirimages dirSTC];
        end 
        if ~exist(dirSTC, 'dir')
          mkdir(dirSTC)
        end
        filenameNofNuclei = [ dirCSV '/NucleiPerSeries_' population(idx).name '.csv'];

        % lettura di una immagine per l'acquisizione delle dimensioni
        DapFilename = sprintf( filenameformat, id_stack, pf, dapiext );
        %dirReadImages = [dirimages population(idx).name '/' dirSeries '/'];
        %if dir_save
        dirReadImages = [dirimages dirold_population_name  ];
        %end
        fieldsize = size(imread( [dirReadImages DapFilename] ));
        % inizializzazione dei vettori contenente gli stack
        height = fieldsize(1); width = fieldsize(2); 
        if (~only_segmentation) %Francesco 20/05/2016
            ImPcG_vol = zeros(height,width,id_planes);
            regIm2_vol = false(height,width,id_planes);
        end
   
        if (~only_thresh && ~segmentation_done)
           fprintf(fileIDt, 'Series%03d ,', population(idx).series(stack));
           fprintf(fileIDt, '%d,%3.8f\n', population(idx).threshplanes(stack),population(idx).thresh);
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

            % reading red Lamin) image
            if lamin
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
    
                % preprocessing of dapi frame by means of mean filter 
                if ( T0 > 0 )
                   w1 = fspecial('average', [T0,T0]);
                   ImDap = imfilter(ImDap,w1); 
                end
    
                % preprocessing of c01 frame by means of mean filter 
                if ( T1 > 0 )
                   w1 = fspecial('average', [T1,T1]);
                   ImPcG = imfilter(ImPcG,w1); 
                end
    
                % preprocessing of c02 frame by means of mean filter 
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

                [regIm, cntIm, pfu, fMeanIn, fMeanOut, den1, den2] = ac_mex(ImDap_seg, ImPcG, ImLam, pfu0, double(SegParameters)); 
    
                % if den1 < den2 
                % den1 == pixels of nuclear regions 
                % den2 == pixels of background 
                % else the opposite case
    
                if (den1(1) > den2(1) ) % eliminare ... ?
                   flag_bckg_color = 1;
                else
                   flag_bckg_color = 0;
                end
        
                %figure
                %imagesc(regIm);colormap(gray);title('before Filling Nuclei');
                % filling holes in dapi regions of segmented image 
                regIm2 = 255 - regIm; % regIm2 has black background and white nuclei regions - regIm is the opposite
                regIm2 = imfill(regIm2,'holes');
                
                %figure
                %imagesc(regIm2);colormap(gray);title('Filling Nuclei');
                %pause
    
                % names of output image files: segmented regions and contours
                tempName = sprintf( filenameformatout, id_stack, plane );
                nameReg  = ['rgn_' tempName '.tif'];
                nameCnt  = ['cnt_' tempName '.tif'];
    
                % saving segmented image in tiff and pgm format
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
    
                % Filtering PcG frames in order to apply ISODATA
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
                imPrimePcG2 = imPrimePcG; 
                imPrimePcG(regIm2 == 0) = -255; 
    
                if (thresh == -255)
                   %fprintf('Initial thresh value == > %f\n',fMeanIn(2));
                   thresh = IsoDataThresh(imPrimePcG, Ny, Nx, fMeanIn, fMeanOut, thresh );
                   fprintf('     Thresh evaluated for PcG Image ==> %f\n',thresh);
                else
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
                % eliminazione dei nuclei attaccati al bordo dell'immgine
                regIm2_vol(:,:,p3d) = imclearborder(regIm2_vol(:,:,p3d),8);
                %%regIm2_vol(:,:,p3d) = imfill(regIm2_vol(:,:,p3d),'holes'); inutile
    
                ImDap_vol(:,:,p3d) = ImDap; % volume della fluorescenza ch00
                ImPcG_vol(:,:,p3d) = ImPcG; % volume della fluorescenza ch01
                ImLam_vol(:,:,p3d) = ImLam; % volume della fluorescenza ch02
    
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
            fieldsize = size(regIm2_vol); % 2 dimensioni

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
       
            if clean_debris
              regIm2_vol = bwareaopen(regIm2_vol, debris_size, 6);
            end%if
  
            % inizializzazione di PSI contenente i dati binari delle segmentazioni
            % PSI avra'  0 per background 1 per le regioni nucleari
            fprintf('Defining PSI for stack %s\n',dirSeries);
            PSI = boolean(regIm2_vol);
            
            % SHED - START
    
            % definizione di kat=bwconncomp(PSI), di default 8-connesse
            % kat e' una struttura generata da bwconncomp contenente le componenti connesse di
            % PSI (8-connessi). kat ha 4 campi: 
            % Connectivity:  connetivita'  utilizzata per individuare le componenti connesse (oggetti/nuclei) 
            % ImageSize:     dimensioni di PSI
            % NumObjects:    Numero di componenti connesse (oggetti/nuclei) di PSI
            % PixelIdxList:  1-X-NumObjects cell array di double, dove il k-esimo array e' un vettore 
            %                contenente gli indici lineari (per colonna) del k-esimo oggetto individuato.
            %                size(kat.PixelIdxList{k},1) fornisce il numero di elementi del k-esimo oggetto
    
            le_dim = [];
            % definizione delle componenti connesse di PSI (procedendo per colonne. 
            % kat contiene la multisuperificie di tutti i nuclei
            kat = bwconncomp(PSI);
            L_PSI = labelmatrix(kat);
            % le_dim e' il vettore contenente la 'dimensione' di ciascun oggetto (nucleo) individuato
            for k = 1:kat.NumObjects
              le_dim(end+1) = size(kat.PixelIdxList{k},1);
            end%fo%r
            %le_dim

            % media delle dimensioni di tutti gli oggetti individuati
            mean_dim = mean(le_dim);
    
            % scarto degli oggetti con dimensione inferiore ad un decimo di quella media 
            maxima_pos = [];
            for k = 1:kat.NumObjects
              %se troppo piccolo fai qualcosa (?? laura)
              % calcolo di maxima_pos contenente l'indice di ''posizione media'' degli oggetti individuati
              [x,y] = ind2sub(size(PSI),kat.PixelIdxList{k});
              %z_dim(end+1) = max(z)-min(z);
    
              if size(kat.PixelIdxList{k},1) > .1*mean_dim
                %%maxima_pos(end+1) = kat.PixelIdxList{k}(  round( size(kat.PixelIdxList{k},1)/2 )  );
                %% [x,y,z] = ind2sub(size(PSI),kat.PixelIdxList{k});
                maxima_pos(end+1,:) = [mean(x),mean(y)];
              end%if
            end%for k
            %idx_PSI = find(le_dim > .1*mean_dim & z_dim > 1);
            idx_PSI = find(le_dim > .1*mean_dim);
            clear le_dim mean_dim;
    
            mask_PSI = ismember(L_PSI, idx_PSI);
            L_PSI(mask_PSI == 0) = 0;
    
            %% Number of cells (nuclei)
            M = size(idx_PSI,2);
            disp(['Identifying cells: ' int2str(M) ' cells found']);

            % determina gli array a1, b1 e c1 contenenti gli indici (i,j,z) (voxel) 
            % a partire dalla matrice
            % degli indici lineari della posizione media dei nuclei, maxima_pos, 
            % ridefinisce maxima_pos di dimensione M X 3 (coordinate del voxel
            % in posizione media nello stack widthXheightXslices di ciascuno degli M nuclei)
            %%[a1,b1,c1]=ind2sub(fieldsize,maxima_pos);
            %%clear maxima_pos
            %%maxima_pos(:,1) = a1; maxima_pos(:,2) = b1; maxima_pos(:,3) = c1;
            %%clear a1 b1 c1
    
            c_iso_value = 1.5 * double(mean(PSI(:)) + min(PSI(:)));
            segm_surfaces_transparency = .5;
            % SHED - END
            % .........................................................................
    
            % memorizza ogni singolo nucleo k=1,...,M in PSI_new(:,:,:,k) individuato dall'etichetta L/L_PSI
            % PSI_new(:,:,:,k) ha dimensione fieldsize
            PSI_new = false([fieldsize,M]);
    
            for k=1:M
               %PSI_new(:,:,:,k) = (L_PSI==idx_PSI(k)); %eliminare
                PSI_new(:,:,k) = (L_PSI==idx_PSI(k));
            end
    
            PSI = boolean(PSI_new);
            clear PSI_new regIm2_vol L L_PSI idx_PSI;
    
            fieldsize = size(PSI);
            disp('Extracting Nuclei...');
     
            % inizializzazione di NCL
            NCL = cell(0);
    
            %if size(fieldsize,2)<4
            %    made = 1;
            %else
            %    made = fieldsize(4);
            %end%if
            made = fieldsize(3);
    
            fprintf('made %d\n',made);
            % ciclo sui nuclei individuati
            for n=1:made
                [i] = find(PSI(:,:,n));
                [x,y] = ind2sub(fieldsize(1:2),i);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% planes cut
                %for plane = pf:pl %cycle on planes
                %    plane 
                %    regImNucleus = PSI( min(x):max(x), min(y):max(y), plane, n );
                %    [LregImNucleus, nCCregImNucleus] = bwlabel(regImNucleus);
                %    if ( nCCregImNucleus == 1 )
                %      regImNucluesstats = regionprops(LregImNucleus, 'Area');
                %      currentArea(plane) = regImNucluesstats(1).Area;
                %      fprintf('nucleo %d current plane %d currentArea %d \n',n, plane, currentArea(plane));
                %    else
                %      fprintf('nucleo %d nCCregImNucleus %d \n',n, nCCregImNucleus);
                %      %cut_can_be_performed = 0; 
                %      currentArea(plane) = 0.;
                %    end
                %end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% planes cut
                NCL{n} = struct;
                % definzione dell'ennesimo nucleo
                NCL{n}.Nucleus = PSI( min(x):max(x), min(y):max(y), n );
                NCL{n}.PcG = ImSegPcG_vol( min(x):max(x), min(y):max(y) ).*NCL{n}.Nucleus;
                % rimuove tra gli spot tutti gli oggetti piu' piccoli di 17 pixel con connettivita' < 6
                %NCL{n}.PcG = bwareaopen(NCL{n}.PcG,17,6); % 6, 18, 26
                NCL{n}.PcG = bwareaopen(NCL{n}.PcG,3,8); % 6, 18, 26
                % Pulire i nuclei (lascia solo l'oggetto maggiore ed i relativi spot)
                % qui avviene la separazione dei nuclei 
                temp = bwconncomp(NCL{n}.Nucleus,6); % 6, 18, 26
                if temp.NumObjects>1
                    thamax = 0; thamaxk = 0;
                    % tra i nuclei individuati estrae la dimensione del piu' grande, thamax ed
                    % il suo indice, thamaxk, tra gli M individuati
                    for k=1:temp.NumObjects
                        if thamax<size(temp.PixelIdxList{k},1)
                            thamax = size(temp.PixelIdxList{k});
                            thamaxk = k;
                        end%if
                    end%for k
                    % annullamento delle componenti spurie dovute alla separazione dei nuclei
                    for k=1:temp.NumObjects
                        if k~=thamaxk
                            NCL{n}.Nucleus(temp.PixelIdxList{k}) = 0;
                        end%if
                    end%for k
                    NCL{n}.PcG = NCL{n}.PcG.*NCL{n}.Nucleus; % remove spot outside of remaining obj
                end%if
                % Rimane sempre un quadratino di 16 pixel nelle slice vuote, non so 
                % perche'. Questo pezzo lo rimuove
                if ~isempty(NCL{n}.Nucleus)
                    void = zeros(size(NCL{n}.Nucleus(:,:,1)));
                    for q=1:size(NCL{n}.Nucleus,3)
                        if  sum(sum(NCL{n}.Nucleus(:,:,q)))<=16
                             NCL{n}.Nucleus(:,:,q) = void;
                        end%if
                    end%for q
                    clear void;
                end%if isempty

                % calcolo dell'eccentricità dei nuclei
                regIm2 = PSI(:,:,n);
                [LregIm2, nCCregIm2] = bwlabel(regIm2);
                regIm2stats = regionprops(LregIm2, 'Area', 'Perimeter', 'PixelList');
 
                fprintf('nucleo %d nCCregIm2 %d \n',n, nCCregIm2);

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
                      %figure, imshow(regIm2_grouped{index_regIm2_grouped});
                    end
                else
                    NCL{n}.eccentricity = 10.; %
                end
                fprintf('nucleo %d area %f perimeter %f metric %f \n',n, ...
                         regIm2stats(i).Area, regIm2stats(i).Perimeter, NCL{n}.eccentricity);
            end%for n

            if create_image | create_map % Tutti i nuclei numerati

                %fig_numbers = figure('Visible', 'off'); imshow(mat2gray(ImDap_vol(:,:,1)));
                fig_numbers = figure('Visible', 'off'); imshow(mat2gray(ImPcG_vol(:,:,1)));
                set(fig_numbers,'Paperposition',[0,0,4,4]);
                %set(fig_numbers, 'Position', [0 0 4 4]);
                %sum(sum(PSI))
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
                %filenameformatnumbers = '%s%03d'; % namestack,stack
                %tempName = sprintf( filenameformatnumbers, id_stack )
                filenameformatnumbers = strcat(stack_format, sprintf('%03d',id_stack));
                print(fig_numbers, '-dpng', '-r200', [ dirSeparatedNuclei '/' filenameformatnumbers '_numbers.png']);
                clear ti tj;
                close(fig_numbers);
            end%if

            clear temp thamax thamaxk PSI ImPcG_vol ImPcG_vol ImLam_vol ImDap_vol;

            C = cell(0); PcG = cell(0);
            face_col = [200,100,0]./255;
            face_alp = .2;

            fileID = fopen(filenameNofNuclei,'a');
            fprintf(fileID, 'Series%03d,', id_stack);
            for n=1:size(NCL,2)
                fprintf('Real Case: 2d reconstruction Nucleus...  %d\n',n);
                fig1 = figure(1); set(fig1, 'Visible', 'on'); clf(1);
                %set(fig1,'Position', [0 0 4 4]);
                set(fig1,'Paperposition',[0,0,4,4]);
                set(fig1,'DoubleBuffer','off');
                % PLOT CENTROID
                [Cx,Cy] = ind2sub(size(NCL{n}.Nucleus),find(NCL{n}.Nucleus));
                C{n} = [mean(Cx),mean(Cy)];
                NCL{n}.NuclCentr = C{n};
                %plot3(C{n}(1),C{n}(2),C{n}(3),'go','MarkerSize',2,'LineWidth',2);
                % PLOT PROTEINS
                PcG{n} = bwconncomp(NCL{n}.PcG,6);
                %PcG{n} = bwconncomp(NCL{n}.PcG,4);
                fprintf('Nucleo n. %d PcG n. %d\n', n, PcG{n}.NumObjects);
                proimage = zeros(size(NCL{n}.Nucleus));
                temp = bwconncomp(NCL{n}.Nucleus,6); % 6, 18, 26
                proimage(temp.PixelIdxList{1}) = 100;
                if PcG{n}.NumObjects > 0
                    karaa = double(label2rgb(1:PcG{n}.NumObjects,'jet',[0,0,0],'shuffle'))./255;
                    for k=1:PcG{n}.NumObjects
                        protemp = zeros(size(NCL{n}.Nucleus));
                        protemp(PcG{n}.PixelIdxList{k}) = 1;
                        [x,y] = ind2sub(size(NCL{n}.Nucleus),PcG{n}.PixelIdxList{k});
                        PcG{n}.PixelSubList{k} = [x,y];
                        PcG{n}.Centroid{k} = [mean(x),mean(y)];
                        PcG{n}.CentDist{k} = sqrt(sum((C{n} - PcG{n}.Centroid{k}).^2));
                        %PATCH_3Darray(protemp,karaa(1,k,:));
                        proimage = proimage + protemp.*255;
                    end%for k
                    % PLOT LINES
                    %for k=1:PcG{n}.NumObjects
                    %    line([C{n}(1),PcG{n}.Centroid{k}(1)], ...
                    %    [C{n}(2),PcG{n}.Centroid{k}(2)], ...
                    %    'Color',karaa(1,k,:)/2);
                    %end%for k
                    imtemp = uint8(mat2gray(proimage).*255);
                else 
                    imtemp = uint8(mat2gray(proimage).*100);
                end%if PcG{n}.NumObjects
                imwrite(imtemp, [dirSeparatedNuclei '/' filenameformatnumbers,  '_' int2str(n) '.png']);
                
                % PLOT MEMBRANE
                %N = regionprops(NCL{n}{1},'Area','MajorAxisLength','MinorAxisLength');
                %p = patch(isosurface(permute(NCL{n}.Nucleus,[2,1,3]),.5,'noshare'));
                % size(NCL{n}.Nucleus,3) 
                %Nucleus_verts = isosurface(permute(NCL{n}.Nucleus,[2,1,3]),.5,'noshare');
                %p = patch(Nucleus_verts);
                %set(p,'FaceColor',face_col,'EdgeColor','none','FaceAlpha',face_alp);
                %patch(isocaps(permute(NCL{n}.Nucleus,[2,1,3])),...
                %'FaceColor',face_col,'EdgeColor','none','FaceAlpha',face_alp);
                %set(gca,'projection','perspective');
                %set(gca,'NextPlot','replaceChildren');
                %daspect([1/xyscale 1/xyscale 1/zscale]), view(3); box('on');
                %lightangle(75,45), %lighting('phong');
                %axis('tight'); axis('vis3d'); grid('on');
                %clear k p protemp x y z;
                %clear h a e az el;
                
                %if create_image % create image
                    %tempName = sprintf( filenameformatnumbers, id_stack );
                    %set(gcf,'PaperSize',[7 7]);
                    %set(gcf,'PaperUnits','inches');
                    %set(gcf,'Paperposition',[0,0,4,4]);
                    %print(gcf, '-dpng', '-r100', [dirrec3d '/' filenameformatnumbers, '_' int2str(n) '_3dview.png']); %prova
                    %print(gcf, '-dpng', '-r200', ['results' regexprep(stackfilez{stackfile},'/','_') '_' int2str(n) '_3dview.png']); %prova
                    %saveas(gcf, [dirresults '/' tempName '_' int2str(n) '_3dview.png']);
                    %view(-37.5,15);
                    %print(gcf, '-dpng', '-r200', ['results/' regexprep(stackfilez{stackfile},'/','_') '_' int2str(n) '_side.png']);
                    %view(90,90);
                    %saveas(gcf, [dirresults '/' tempName  '_' int2str(n) '_top.png']);
                    %print(gcf, '-dpng', '-r200', ['results/' regexprep(stackfilez{stackfile},'/','_') '_' int2str(n) '_top.png']);
                    %print(gcf, '-dpng', '-r100', [dirrec3d '/' filenameformatnumbers,  '_' int2str(n) '_top.png']);
                    %close(fig1);
                %end%if create image
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
        end%if(~only_segmentation)
        if (~only_thresh && ~segmentation_done)
           fprintf(fileID, '\n');
           fclose(fileID);
        end
    end%for stack
    if (~only_thresh && ~segmentation_done)
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
