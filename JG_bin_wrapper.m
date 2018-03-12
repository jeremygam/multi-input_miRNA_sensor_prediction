% Written by Jeremy Gam (jgam@mit.edu)
% This is used to analyze FCS files for miRNA sensor data by binning,
% fitting model to the data to get M and Km, and plotting if necessary

%% Clear items

clear all
close all
clc

addpath('path/to/other/m-files')

%% Define parameters

FCS_dir = 'path/to/FCS/files'; %Dir with only relevant FCS files and scripts
no_dna_control_filename = 'path/to/FCS/files/NoDNA_control.fcs';

bin_color_name = 'Blue'; %Blue=Pacific Blue, Cyan=AmCyan
reporter_color_name = 'Red'; %Red=Texas Red

% Get channel numbers
red_channel_name   = 'Red';
green_channel_name = 'FIT'; %FIT = FITC
blue_channel_name  = 'Blue';

T_logicle = 262144;
M_logicle = 4.5;
r_logicle = -150;
W_logicle = (M_logicle-log10(T_logicle/abs(r_logicle)))/2;
A_logicle = 0;

%% Use parameters

% Open FCS file
[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(no_dna_control_filename);

red_channel_num    = getChannel(fcshdr,'Red');
yellow_channel_num = getChannel(fcshdr,'FIT');
blue_channel_num   = getChannel(fcshdr,'Blue');

bin_color_num      = getChannel(fcshdr,bin_color_name);
reporter_color_num = getChannel(fcshdr,reporter_color_name);


%% Grab the directory and the files within

files = dir(FCS_dir);
filenames = cell(0); 
for i=1:length(files)
    if length(files(i).name) > 3 && strcmp(files(i).name(end-3:end),'.fcs')==1 && files(i).name(1) ~= '.'
        filenames = [filenames files(i).name];
        drawnow
    end
end

%% Read in no dna control FCS file and gate

% Open FCS file
[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(no_dna_control_filename);

% Get channel numbers
red_channel_num   = getChannel(fcshdr,'Red');
yellow_channel_num = getChannel(fcshdr,'FIT');
blue_channel_num  = getChannel(fcshdr,'Blue');

% Gate the cells based on FSC and SSC
fcsdat_gated_by_no_DNA = applyJCGate_JG(no_dna_control_filename,no_dna_control_filename);

close all

%% Use more parameters

%axes_min = 1;
%axes_max = 2.5e+5;
num_bins = 20;

x_data_temp = fcsdat_gated_by_no_DNA(:,bin_color_num);
y_data_temp = fcsdat_gated_by_no_DNA(:,reporter_color_num);

x_data_min = min(x_data_temp);
x_data_max = 50000;
x_data_min = logicleTransform(x_data_min,T_logicle,W_logicle,M_logicle,A_logicle) * M_logicle;
x_data_max = logicleTransform(x_data_max,T_logicle,W_logicle,M_logicle,A_logicle) * M_logicle;

x_data_bin = linspace(x_data_min, x_data_max, num_bins);
%y_data_bin = zeros(1,length(x_data_bin));
%y_bin_ind = discretize(x_data_temp,x_data_bin);

% Define bins
edges = x_data_bin;

%% Define rate constants

%{
Example rate constants and sources:
k_trs = 7; %mRNAs/hour. (Range = 5.8-8.7 mRNAs/hour from Darzacq et al. Nat. Struct. Mol. Biol. 2007. and Schwanhausser et al. Nature 2011).
mRNA_half_life = 240; %min. (Range = 3-5 hours half life Raj et al. PLoS Bio 2006).
mRNA_half_life = 100; %min. (Range = 100 min half life Dar et al. PNAS 2012).
k_deg_m = 0.693 / mRNA_half_life;
k_tln = 1000; %proteins per mRNA per hour. (Range = 750-1300 for high expression from Schwanhausser et al. Nature 2011).
marker_half_life = 2.8; %hours. (Range = 2.8 hour from Halter et al. Cytometry 2007).
marker_half_life = 2.5; %hours. (Range = 2.5 hour from Dar et al. PNAS 2007).
k_deg_marker = 0.693 / marker_half_life;
k_deg_reporter = k_deg_marker;
k_off = 8.8e-5; %sec^-1. (Range = 8.8e-5 sec^-1 from Wee et al. Cell 2012).
K_d = 3.7e-12; %M. (Range = 3.7+/-0.9 pM from Wee et al. Cell 2012).
k_on = k_off/K_d; %1/(M*sec).
k_cat = 6.1e-2; %sec^-1. (Range = 6.1e-2 sec^-1 from Wee et al. Cell 2012).
K_m = 8.4; (nM) Haley and Zamore Nature Struct Molec Biol 2004 
K_m = 8.4e-9*1.54e-12*6.022e23;

mRNA half life = 7 hr Sacchetti et al FEBS letters March 2001
Protein half life = 50 hr Sacchetti et al FEBS letters March 2001
g_R = 0.693/7 = 0.1 1/hr
g_P = = 0.693/50 = 0.014 1/hr

kcat = 7.1e-3 (1/sec) = 25.6 1/hr Haley and Zamore Nature Struct Molec Biol 2004

copies/nucleus = ~1.55 * fluorescence Extrapolated from Cohen et al J
Control Release 2009

k_R = determined by comparing copies/nucleus to mRNA levels from qPCR
k_P = determined by comparing fluorescence to mRNA levels from qPCR
alpha = k_R/g_R
beta = k_cat/g_R = 25.6/0.1 = 256
gamma = g_R * g_P / (k_R * k_P)
delta = k_P/g_P

%}

%old parameters
k_R1 = 7; % mRNA/hour
k_R2 = 7; % mRNA/hour
g_R1 = 0.5; % 1/hour    increase
g_R2 = 0.5; % 1/hour    increase
k_P1 = 8; %protein/mRNA/hour   reduce
k_P2 = 5; %protein/mRNA/hour   reduce
g_P1 = 0.5; % 1/hour   increase
g_P2 = 0.5; % 1/hour   increase
g_R1M = 5; % 1/hour   reduce

K_m = 7.8e3; %[400 4000 40000]; % molec/cell vol    reduce to drag down transition
alpha = k_R1/g_R1;
beta = g_R1M/g_R1;
gamma = g_R2*g_P2/(k_R2*k_P2);
delta = k_P1/g_P1;

%% Read in single-input FCS files in directory and save to struct

number_of_bins = 20;

% Define struct where we will save the data
%binned_data = [];

for i = 1:length(filenames)
    % Define the filepath
    temp_path = strcat(FCS_dir,'/',filenames{i});
    filenames{i}
    
    % Gate the FCS file
    fcsdat_gated_by_no_DNA = applyJCGate_JG(temp_path,no_dna_control_filename);
    
    % Logicle transform fcs data
    fcsdat_transform = logicleTransform(fcsdat_gated_by_no_DNA,T_logicle,W_logicle,M_logicle,A_logicle) * M_logicle;
        
    % Determine which colors are important
    x_raw = fcsdat_gated_by_no_DNA(:,bin_color_num);
    y_raw = fcsdat_gated_by_no_DNA(:,reporter_color_num);
    x_logicle = fcsdat_transform(:,bin_color_num);
    y_logicle = fcsdat_transform(:,reporter_color_num);
    
    % Bin in case it is necessary later: bin raw and transformed data using the x_logicle data set into n bins
    [x_bin,y_bin,x_bin_raw,y_bin_raw] = binning(x_logicle,y_logicle,x_raw,y_raw,number_of_bins);
    
    % Save binned data to struct
    temp_filename = strrep(filenames{i}(1:end-8),'-','_');  
    binned_data.(temp_filename).x_data_raw = x_raw;
    binned_data.(temp_filename).y_data_raw = y_raw;   
    binned_data.(temp_filename).x_data_logicle = x_logicle;
    binned_data.(temp_filename).y_data_logicle = y_logicle; 
    binned_data.(temp_filename).x_data_raw_bin = x_bin_raw;
    binned_data.(temp_filename).y_data_raw_bin = y_bin_raw;      
    binned_data.(temp_filename).x_data_logicle_bin = x_bin;
    binned_data.(temp_filename).y_data_logicle_bin = y_bin;  
    binned_data.(temp_filename).num_cells = length(x_raw);
    
    display(binned_data)
    
end

close all


%% Add information before plotting

% Definitions in the following format: sample_number sample_name description
sample_defs  =  {
    '1',  'sample_001',  'example_description'
    };
             
% Add in descriptions to binned data struct
for i = 1:size(sample_defs,1)
    binned_data.(sample_defs{i,2}).description = sample_defs{i,3};
end

% For plotting logicle axes
linear_axis = [linspace(-100, 90, 20) linspace(100, 900, 9) linspace(1000, 9000, 9) linspace(10000, 90000, 9) linspace(100000, 900000, 9)];
logicle_axis = M_logicle * logicleTransform(linear_axis,T_logicle,W_logicle,M_logicle,A_logicle);
linear_axis_labels = cell(1,length(linear_axis));
linear_axis_labels{1} = '-10^2';
linear_axis_labels{11} = '0';
linear_axis_labels{21} = '10^2';
linear_axis_labels{30} = '10^3';
linear_axis_labels{39} = '10^4';
linear_axis_labels{48} = '10^5';

%% Determine M and Km of the single-input sensors

filenames = sample_defs(:,2);

for i = 1:length(filenames)
    
    disp(filenames{i})
    
    % Logicle transform
    x_transform = binned_data.(filenames{i}).x_data_logicle;
    y_transform = binned_data.(filenames{i}).y_data_logicle;
    M_guess = 10;
    Km_guess = 10;
    
    % Fit Km and M to the file
    [new_param, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(@(x,x_transform) model_and_transform_miRNA(x,x_transform,alpha,beta,gamma,delta,T_logicle,W_logicle,M_logicle,A_logicle) ,[Km_guess M_guess],x_transform,y_transform,[0 0],[10000 10000]);
    conf = nlparci(new_param, residual, 'jacobian', jacobian);

    % Generate curves from fit parameters for reference
    P_2 = [ linspace(-500,0.99,100) logspace(0,7,500)];
    Km_temp = new_param(1);
    M_temp = new_param(2);
    P_1 = delta/2*(sqrt((-alpha*gamma*P_2 + beta*M_temp + Km_temp).^2 + 4*alpha*gamma*P_2*Km_temp) + alpha*gamma*P_2 - beta*M_temp - Km_temp);
    P_2 = logicleTransform(P_2,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;
    P_1 = logicleTransform(P_1,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;
    
    binned_data.(filenames{i}).K_m_final_kP1_fixed = Km_temp;
    binned_data.(filenames{i}).M_final_kP1_fixed = M_temp;
    binned_data.(filenames{i}).x_modeled_kP1_fixed = P_2;
    binned_data.(filenames{i}).y_modeled_kP1_fixed = P_1; 
    
    
    %{
    mir_name = filenames{i}(1:length(filenames{i})-8);
    mir_name = strrep(mir_name,'-','_');
    temp_struct = struct('x_data_logicle',x_transform,'y_data_logicle',y_transform,'x_fit_logicle',P_2,'y_fit_logicle',P_1,'Km',new_param(1),'M',new_param(2),'resnorm',resnorm,'residual',residual,'conf_int',conf,'jacobian',jacobian);
    all_data.(mir_name) = temp_struct;
    %}


end

%% Save binned data struct to mat file

save([FCS_dir '/binned_data.mat'],'binned_data')

%% Load structs

load([FCS_dir '/binned_data.mat'],'binned_data')

%% Do plotting
% Plotting based on Patrick Martineau subplots
% (http://p-martineau.com/perfect-subplot-in-matlab/)

% Define how you want to plot the data. Each row is a separate plot and 
%numbers are from sample_defs
plot_defs = ...
    [1];

% Parameters for figure and panel size
x_axis_name = 'EBFP2';
y_axis_name = 'mKate2';
plotheight=30;
plotwidth=24;
subplotsx=3;
subplotsy=1;   
leftedge=3;
rightedge=0.4;   
topedge=1;
bottomedge=3.3;
spacex=0.2;
spacey=0.2;
fontsize=10;    
sub_pos_JG=subplot_pos_JG(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

% Setting the Matlab figure
f=figure('visible','on')
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

% Color order
color_order =          [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

%legend_combine = cell(1,4);
handle_vector = [];

% Loop to create axes
count = 1;
for i = 1:subplotsy
    for j = 1:subplotsx
        ax = axes('position',sub_pos_JG{i,j},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        
        display(sub_pos_JG{i,j})
        
        % Do a single plot        
        plot_def_temp = plot_defs(count, :);
        display(plot_def_temp)
        count = count + 1;
        
        for k = 1:size(plot_def_temp,2) 
            hold on
            sample_name_temp = sample_defs{plot_def_temp(k),2}
            x_data_temp = binned_data.(sample_name_temp).x_data_logicle_bin;
            y_data_temp = binned_data.(sample_name_temp).y_data_logicle_bin;
            %x_model_temp = binned_data.(sample_name_temp).x_modeled_kP1_fixed;
            %y_model_temp = binned_data.(sample_name_temp).y_modeled_kP1_fixed;
            
            legend_name = binned_data.(sample_name_temp).description;
            %legend_combine{1,k} = legend_name;
            
            handlevector(10*k) = plot(x_data_temp,y_data_temp,'o-','Color',color_order(k,:));
            %handlevector(10*k+1) = plot(x_model_temp,y_model_temp,'-','Color',color_order(k,:));
            axis([0.1 4 0.1 4])
            fig_handle = gca;
            fig_handle.XTick = logicle_axis;
            fig_handle.YTick = logicle_axis;
            fig_handle.XTickLabel = linear_axis_labels;
            fig_handle.YTickLabel = linear_axis_labels;
            grid on
        end
        hold off
        %legend(handlevector([10 20 30 40]),legend_combine,'Location','northwest')

%        if j==subplotsy
%            title(['Title (',num2str(i),')'])
%        end
        
        if i<subplotsx
            set(ax,'xticklabel',[])
        end
        
        if j>1
            set(ax,'yticklabel',[])
        end
        
        if i==subplotsx
            xlabel(x_axis_name)
        end
        
        if j==1
            ylabel(y_axis_name)
        end
        
    end
end