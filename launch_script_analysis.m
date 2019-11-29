% nome del dataset (data)
popol = ''; % useful if you have more populations for a given data set
% nome della popolazione sottoforma di nickname (utilizzato anche per il nome della heatmap)
popnick_list = {'2019-03-08_C001Ezh2'}; %,...
                %'2019-03-08_C002Ezh2',...
                %'2019-03-08_C004Ezh2',...
                %'2019-03-08_HGPS164p11Ezh2',...
                %'2019-03-08_HGPS169p11Ezh2',...
                %'2019-03-08_HGPS169p20Ezh2'};

outCellTot = cell(0);

outCellTot{1} = [0; 0; 0; 0; 0;  0;  0;  0;  0;  0;];
%outCellTot{2} = [0; 0; 0; 0; 0;  0;  0;  0;  0;  0;];
%outCellTot{3} = [0; 0; 0; 0; 0;  0;  0;  0;  0;  0;];
%outCellTot{4} = [0; 0; 0; 0; 0;  0;  0;  0;  0;  0;];
%outCellTot{5} = [0; 0; 0; 0; 0;  0;  0;  0;  0;  0;  0;  0;  0;  0;];
%outCellTot{6} = [0; 0; 0; 0; 0;  0;];

initpath = ['']; % final / is needed
outpath = ['2019-03-08_C001Ezh2/']; % final / is needed

titlef = popol;
legendas = {'C001Ezh2'}; %,...
            %'C002Ezh2',...
            %'C004Ezh2',...
            %'HGPS164p11Ezh2',...
            %'HGPS169p11Ezh2',...
            %'HGPS169p20Ezh2'};

iprint = 1;

f11_ConfrontMultipleCases_2D(popol, popnick_list, outCellTot, initpath, outpath, titlef, legendas, iprint);
f16_mindistancesPcG2NCL(popol, popnick_list, outCellTot, initpath, outpath, iprint);
