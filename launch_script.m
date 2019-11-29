dirimages = 'images/C001_Ezh2488LamAC647/';
dirold_population_name = '';
stack_format_list = {'C001Ezh2_Series'};
xyscale = 0.241;
zscale = 1.;
population(1).name = '2019-03-08_C001Ezh2'; % output directory
population(1).series =       [7 9 11 13 15 17 19 21 23 25];
population(1).planes =       [1 1  1  1  1  1  1  1  1  1];
population(1).threshplanes = [0 0  0  0  0  0  0  0  0  0];
population(1).channels = 3;
population(1).thresh = -255;
segmentation_done = 0;
only_segmentation = 0; % 1 nel caso anche di calcolo del thresh
only_thresh = 0;
random_case = 0;
dir_save = 0;
flag_seg = 2;
mu = 1;
lambda = 0.0001;
lamin = 1;

mainSeg_pcg2D(population,segmentation_done,only_segmentation,dirimages,stack_format_list,only_thresh,...
            random_case,xyscale,zscale,dirold_population_name,dir_save,flag_seg,mu,lambda,lamin);
