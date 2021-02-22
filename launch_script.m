tic
dirimages = 'images/C001_Ezh2488LamAC647/';
population(1).name = '2019-03-08_C001Ezh2'; % output directory
population(1).xyscale = 0.241;
population(1).stack_basename = 'C001Ezh2_Series';
population(1).stack_format = '%03d_z%02d_ch%02d.tif';
population(1).series =       [7 9 11 13 15 17 19 21 23 25];
population(1).aplane =       [0 0  0  0  0  0  0  0  0  0];
population(1).channels = 3;
population(1).thresh = -255;
population(1).chext = [0 1 2];
segmentation_done = 0; % 1 if segmentation already performed
only_segmentation = 0; % 1 to perform only segmentation
print_thresh = 1;
dir_save = 0;
flag_seg = 2;
mu = 1;
lambda = 0.0001;

mainSeg_pcg2D(population,segmentation_done,only_segmentation, ...
              dirimages,print_thresh,dir_save,flag_seg,mu,lambda);
toc
