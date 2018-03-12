% Written by Jeremy Gam (jgam@mit.edu)
% Combined script to analyze single-input miRNA sensor library data by
% fitting, plot theoretical curves from the model, calculate residuals

%% Clear items

clear all
close all
clc

addpath('/Users/jeremygam/Dropbox (MIT)/Weiss Lab Docs/Scripts/Constant Analysis Scripts')


%% Define parameters

FCS_dir = '/Users/jeremygam/Desktop/to_analyze/Renamed_LSB_Screen_HeLa/All Plates'; %Dir with FCS files
control_name = '/Users/jeremygam/Desktop/to_analyze/Renamed_LSB_Screen_HeLa/All Plates/HeLa_no_dna_047.fcs'; %Dir of no miRNA target control FCS file
cd(FCS_dir)

T_logicle = 262144;
M_logicle = 4.5;
r_logicle = -10;
W_logicle = (M_logicle-log10(T_logicle/abs(r_logicle)))/2;
A_logicle = 0;


%% Use some parameters

%{
[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(control_name);

red_channel_num    = getChannel(fcshdr,'Red');
yellow_channel_num = getChannel(fcshdr,'FIT');
blue_channel_num   = getChannel(fcshdr,'Blue');

marker_color = blue_channel_num;
report_color = red_channel_num;
%}

marker_color = 8; % Pacific Blue = 7; AmCyan = 8
report_color = 10; % Texas Red = 10


%% Define rate constants

%{
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

k_R1 = 7; % mRNA/hour
k_R2 = 7; % mRNA/hour
g_R1 = 0.5; % 1/hour
g_R2 = 0.5; % 1/hour
k_P1 = 8; %protein/mRNA/hour
k_P2 = 5; %protein/mRNA/hour
g_P1 = 0.5; % 1/hour
g_P2 = 0.5; % 1/hour
g_R1M = 5; % 1/hour

K_m = 7.8e3; %[400 4000 40000]; % molec/cell vol    reduce to drag down transition
alpha = k_R1/g_R1;
beta = g_R1M/g_R1;
gamma = g_R2*g_P2/(k_R2*k_P2);
delta = k_P1/g_P1;


%% Read in control file

% Open FCS file
[fcsdat, fcshdr] = fca_readfcs(control_name);

% Gate the cells based on FSC and SSC
fcsdat_gated_by_no_DNA = applyJCGate_JG(control_name,control_name);


%% Put all raw and logicle data into struct and .mat file

% Grab the directory and the files within
files = dir(FCS_dir);

% Determine all names for fcs files in the directory
filenames = cell(0);
for i=1:length(files)
    %disp(files(i).name)
    if length(files(i).name) > 3 && strcmp(files(i).name(end-3:end),'.fcs')==1
        filenames = [filenames files(i).name];
    end
end

% Define struct where we will save the data
all_data = [];

for i = 1:length(filenames)
    
    % Define the filepath
    temp_path = strcat(FCS_dir,'/',filenames{i});
    disp(temp_path)    
    
    % Open FCS file
    [fcsdat, fcshdr] = fca_readfcs(temp_path);
    
    % Gate the FCS file
    fcsdat_gated_by_no_DNA = applyJCGate_JG(temp_path,control_name);
    
    % Determine which colors are important
    x = fcsdat_gated_by_no_DNA(:,marker_color);
    y = fcsdat_gated_by_no_DNA(:,report_color);
    
    % Logicle transform
    x_transform = logicleTransform(x,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;
    y_transform = logicleTransform(y,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;

    % Save to struct
    mir_name = filenames{i}(1:length(filenames{i})-8);
    mir_name = strrep(mir_name,'-','_');
    temp_struct = struct('x_data_logicle',x_transform,'y_data_logicle',y_transform,'x_data_raw',x,'y_data_raw',y);
    all_data.(mir_name) = temp_struct;

end

% Save the data
save('all_data.mat','all_data')


%% Load raw data into all_data struct

load('all_data.mat')


%% Bin data by EBFP2

% Generate all mirna names for FCS files in the directory
mirs_of_interest = fieldnames(all_data);

% Loop through all miRNAs for analysis
for i = 1:length(mirs_of_interest)
    disp(mirs_of_interest{i})    
    
    
    x_data_temp = all_data.(mirs_of_interest{i}).x_data_logicle;
    y_data_temp = all_data.(mirs_of_interest{i}).y_data_logicle;
    
    x_data_min = min(x_data_temp);
    x_data_max = max(x_data_temp);
    x_data_bin = linspace(x_data_min, x_data_max, 18);
    y_data_bin = zeros(1,length(x_data_bin));
    
    y_bin_ind = discretize(x_data_temp,x_data_bin);
    
    for j = 1:length(x_data_bin)
        y_in_bin = y_data_temp(y_bin_ind == j);
        y_data_bin(j) = median(y_in_bin);
    end
    
    %figure(i)
    %plot(x_data_temp,y_data_temp,'ro',x_data_bin,y_data_bin,'bo')
    
    all_data.(mirs_of_interest{i}).x_data_bin = x_data_bin;
    all_data.(mirs_of_interest{i}).y_data_bin = y_data_bin;
    
end

%% Save data

save('all_data.mat','all_data')


%% Load data from all_data struct

load('all_data.mat')


%% Fit parameters to data

filenames = fieldnames(all_data);

for i = 1:length(filenames)
        
    % Logicle transform
    x_transform = all_data.(filenames{i}).x_data_logicle;
    y_transform = all_data.(filenames{i}).y_data_logicle;
    %x_transform = all_data.(filenames{i}).x_data_bin; %can fit to binned/unbinned
    %y_transform = all_data.(filenames{i}).y_data_bin;
    M_guess = 100;
    Km_guess = 10;
    
    % Fit Km and M to the file
    [new_param, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(@(x,x_transform) model_and_transform_miRNA(x,x_transform,alpha,beta,gamma,delta,T_logicle,W_logicle,M_logicle,A_logicle) ,[Km_guess M_guess],x_transform,y_transform,[0 0],[10000 10000]);
    conf = nlparci(new_param, residual, 'jacobian', jacobian);

    % Generate curves from fit parameters for reference
    P_2 = [ linspace(-500,1,100) logspace(0,7,500)];
    Km_temp = new_param(1);
    M_temp = new_param(2);
    P_1 = delta/2*(sqrt((-alpha*gamma*P_2 + beta*M_temp + Km_temp).^2 + 4*alpha*gamma*P_2*Km_temp) + alpha*gamma*P_2 - beta*M_temp - Km_temp);
    P_2 = logicleTransform(P_2,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;
    P_1 = logicleTransform(P_1,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;
    
    all_data.(filenames{i}).K_m_final_kP1_fixed = Km_temp;
    all_data.(filenames{i}).M_final_kP1_fixed = M_temp;
    all_data.(filenames{i}).x_modeled_kP1_fixed = P_2;
    all_data.(filenames{i}).y_modeled_kP1_fixed = P_1; 
    
    %{
    mir_name = filenames{i}(1:length(filenames{i})-8);
    mir_name = strrep(mir_name,'-','_');
    temp_struct = struct('x_data_logicle',x_transform,'y_data_logicle',y_transform,'x_fit_logicle',P_2,'y_fit_logicle',P_1,'Km',new_param(1),'M',new_param(2),'resnorm',resnorm,'residual',residual,'conf_int',conf,'jacobian',jacobian);
    all_data.(mir_name) = temp_struct;
    %}

end
save('all_data.mat','all_data')


%% Save the best Km and M values to .mat file

save('final_fit.mat','all_data')

%% Load data if necessary

load('final_fit.mat')


%% Export best Km and M values to excel
%{
% Output Km and M to excel

names = fieldnames(all_data);
to_export = cell(length(names)+1,6);

to_export{1,1} = 'miRNA name';
to_export{1,4} = 'error';
to_export{1,5} = 'Km_final';
to_export{1,6} = 'M_final';


for i = 1:length(names)

    to_export{i+1,1} = names{i};
    to_export{i+1,2} = all_data.(names{i}).K_m_best_global_kP1_8;
    to_export{i+1,3} = all_data.(names{i}).M_best_global_kP1_8;
    to_export{i+1,4} = all_data.(names{i}).error_global_kP1_8;
    to_export{i+1,5} = all_data.(names{i}).K_m_final;
    to_export{i+1,6} = all_data.(names{i}).M_final; 

    
end


%excel_out = fin
file_out = [FCS_dir, '/final_fit_kP1_8_logicle.csv'];
cell2csv(file_out, to_export)
%}


%% Find and save residuals for fits from using single iteration of fitting M and Km

%{
mirs_of_interest = {'HEK293FT_hsa_miR_106a_3p'
                    'HEK293FT_hsa_miR_19a_3p'
                    'HEK293FT_hsa_miR_106a_5p'
                    'HEK293FT_hsa_miR_19b_3p'
                    'HEK293FT_hsa_miR_183_5p'
                    'HEK293FT_hsa_miR_323a_5p'};                
%}

mirs_of_interest = fieldnames(all_data);

                
for i = 1:length(mirs_of_interest)

    % Generate curve for model fit    
    Km_temp = all_data.(mirs_of_interest{i}).K_m_final_kP1_fixed;
    M_temp  = all_data.(mirs_of_interest{i}).M_final_kP1_fixed;
    x_fit = [ linspace(-500,1,100) logspace(0,7,500)];
    y_fit = delta/2*(sqrt((-alpha*gamma*x_fit + beta*M_temp + Km_temp).^2 + 4*alpha*gamma*x_fit*Km_temp) + alpha*gamma*x_fit - beta*M_temp - Km_temp);
    x_fit = logicleTransform(x_fit,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;
    y_fit = logicleTransform(y_fit,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;

    x_data_norm = all_data.(mirs_of_interest{i}).x_data_logicle;
    y_data_norm = all_data.(mirs_of_interest{i}).y_data_logicle;
    
    % Cull data points with top 10 fluorescences to avoid logicle
    % transform problems (doesn't converge)
%     [x_data_norm, sort_index] = sort(x_data_norm);
%     y_data_norm = y_data_norm(sort_index);
%     x_data_norm = x_data_norm(1:length(x_data_norm)-10);
%     y_data_norm = y_data_norm(1:length(y_data_norm)-10);
    
    % Generate residuals 
    f = @(x,x_data_norm) model_and_transform_miRNA(x,x_data_norm,alpha,beta,gamma,delta,T_logicle,W_logicle,M_logicle,A_logicle);
    y_model = f([Km_temp M_temp],x_data_norm);
    residual = y_data_norm - y_model;
    
    % Generate distances instead of residuals 
%     f = @(x,x_data_norm) model_and_transform_miRNA(x,x_data_norm,alpha,beta,gamma,delta,T_logicle,W_logicle,M_logicle,A_logicle);
%     y_model = f([Km_temp M_temp],x_fit);
%     dist_matrix = pdist2([x_data_norm y_data_norm],[x_fit' y_fit']);
%     [dist, ind] = min(dist_matrix,[],2);
%     dist = dist.*((y_data_norm>y_fit(ind)')-0.5)*2;
%     x_nearest_model = x_fit(ind);
    
    all_data.(mirs_of_interest{i}).x_residuals = x_data_norm;
    all_data.(mirs_of_interest{i}).y_residuals = residual;

    
end

save('final_fit_kP1_fixed.mat','all_data')


%% Plot residuals for fits from using single iteration of fitting M and Km


mirs_of_interest = {'HEK293FT_hsa_miR_106a_3p'
                    'HEK293FT_hsa_miR_15a_5p'
                    'HEK293FT_hsa_miR_125b_5p'};                

%mirs_of_interest = {'HEK293FT_hsa_miR_19a_3p'};

                
for i = 1:length(mirs_of_interest)

    % Generate curve for model fit    
    Km_temp = all_data.(mirs_of_interest{i}).K_m_final_kP1_fixed;
    M_temp  = all_data.(mirs_of_interest{i}).M_final_kP1_fixed;
    x_fit = [ linspace(-500,1,100) logspace(0,7,500)];
    y_fit = delta/2*(sqrt((-alpha*gamma*x_fit + beta*M_temp + Km_temp).^2 + 4*alpha*gamma*x_fit*Km_temp) + alpha*gamma*x_fit - beta*M_temp - Km_temp);
    x_fit = logicleTransform(x_fit,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;
    y_fit = logicleTransform(y_fit,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;

    x_data_norm = all_data.(mirs_of_interest{i}).x_data_logicle;
    y_data_norm = all_data.(mirs_of_interest{i}).y_data_logicle;
    
    % Cull data points with top 10 fluorescences to avoid logicle
    % transform problems (doesn't converge)
%     [x_data_norm, sort_index] = sort(x_data_norm);
%     y_data_norm = y_data_norm(sort_index);
%     x_data_norm = x_data_norm(1:length(x_data_norm)-10);
%     y_data_norm = y_data_norm(1:length(y_data_norm)-10);
    
    % Generate residuals 
    f = @(x,x_data_norm) model_and_transform_miRNA(x,x_data_norm,alpha,beta,gamma,delta,T_logicle,W_logicle,M_logicle,A_logicle);
    y_model = f([Km_temp M_temp],x_data_norm);
    residual = y_data_norm - y_model;
    
    % Generate distances instead of residuals 
    f = @(x,x_data_norm) model_and_transform_miRNA(x,x_data_norm,alpha,beta,gamma,delta,T_logicle,W_logicle,M_logicle,A_logicle);
    y_model = f([Km_temp M_temp],x_fit);
    dist_matrix = pdist2([x_data_norm y_data_norm],[x_fit' y_fit']);
    [dist, ind] = min(dist_matrix,[],2);
    dist = dist.*((y_data_norm>y_fit(ind)')-0.5)*2;
    x_nearest_model = x_fit(ind);
    
%     x_data_bin = all_data.(mirs_of_interest{i}).x_data_bin;
%     y_data_bin = all_data.(mirs_of_interest{i}).y_data_bin;
    
    
    subplot(length(mirs_of_interest),3,3*i-2)
    plot(x_data_norm,y_data_norm,'.',x_fit,y_fit,'r')
    xlabel('Transf. Marker Fluorescence (logicle transformed)')
    ylabel('Reporter Fluorescence (logicle transformed)')
    legend('miRNA Sensor Data','Modeled miRNA Activity','Location','northwest')
    handle = title(strrep(mirs_of_interest{i},'_','-'));
    positions = get(handle,'Position');
    set(handle,'rotation',90,'position',[-3.5 3 3])
    
    %subplot(1,4,2)
    %plot(x_nearest_model,dist,'o')
    
    %subplot(1,4,3)
    %histogram(dist)
    
    subplot(length(mirs_of_interest),3,3*i-1)
    plot(x_data_norm,residual,'.')
    xlabel('Transf. Marker Fluorescence (logicle transformed)')
    ylabel('Residual Fluorescence')
    
    subplot(length(mirs_of_interest),3,3*i)
    %histogram(residual)
    hist(gca,residual,50)
    set(gca,'view',[90 -90])
    xlabel('Histogram of Residuals')

    
end


%% Plot various theoretical curves for varying M (single miRNA low sensor)

P_2 = logspace(1,8,500);

K_m = 10000000;


for M = [linspace(0, 9, 9) logspace(1, 7, 60)]
    P_1 = delta/2*(sqrt((-alpha*gamma*P_2 + beta*M + K_m).^2 + 4*alpha*gamma*P_2*K_m) + alpha*gamma*P_2 - beta*M - K_m);
    
    % Logicle transform fcs data
    P_1_transform = M_logicle * logicleTransform(P_1,T_logicle,W_logicle,M_logicle,A_logicle);
    P_2_transform = M_logicle * logicleTransform(P_2,T_logicle,W_logicle,M_logicle,A_logicle);
    %P_1_transform = lin2logicle (P_1);
    %P_2_transform = lin2logicle(P_2);
    
    
    %fcsdat_transform = fcsdat;
    figure(40)
    plot(P_2_transform,P_1_transform)
    hold on

    figure(50)
    loglog(P_2,P_1)

    %xlim([1 1e8])
    %ylim([1 1e8])
    hold on
end
figure(40)
hold off
figure(50)
hold off


%% Plot various theoretical curves for varing Km (single miRNA low sensor)

T_logicle = 262144;
M_logicle = 4.5;
r_logicle = -10;
W_logicle = (M_logicle-log10(T_logicle/abs(r_logicle)))/2;
A_logicle = 0;

%old parameters
k_R1 = 7; % mRNA/hour
k_R2 = 7; % mRNA/hour
g_R1 = 0.5; % 1/hour    increase
g_R2 = 0.5; % 1/hour    increase
k_P1 = 10; %protein/mRNA/hour   reduce
k_P2 = 5; %protein/mRNA/hour   reduce
g_P1 = 0.5; % 1/hour   increase
g_P2 = 0.5; % 1/hour   increase
g_R1M = 5; % 1/hour   reduce

K_m = 7.8e3; %[400 4000 40000]; % molec/cell vol    reduce to drag down transition
alpha = k_R1/g_R1;
beta = g_R1M/g_R1;
gamma = g_R2*g_P2/(k_R2*k_P2);
delta = k_P1/g_P1;

%{
mirs_of_interest = {'HEK293FT_hsa_miR_106a_3p'
'HEK293FT_hsa_miR_19a_3p'
'HEK293FT_hsa_miR_106a_5p'
'HEK293FT_hsa_miR_19b_3p'
'HEK293FT_hsa_miR_20a_5p';
%}

%mirs_of_interest = 'HEK_LSB_FF4_B3_B03';
mirs_of_interest = 'HEK293FT_hsa_miR_93_5p';

P_2 = logspace(1,6,500);
M = 316;

for K_m = logspace(-2,6,50)
    P_1 = delta/2*(sqrt((-alpha*gamma*P_2 + beta*M + K_m).^2 + 4*alpha*gamma*P_2*K_m) + alpha*gamma*P_2 - beta*M - K_m);
    
    % Logicle transform fcs data
    P_1_transform = M_logicle * logicleTransform(P_1,T_logicle,W_logicle,M_logicle,A_logicle);
    P_2_transform = M_logicle * logicleTransform(P_2,T_logicle,W_logicle,M_logicle,A_logicle);
    %P_1_transform = lin2logicle (P_1);
    %P_2_transform = lin2logicle(P_2);
    
    
    %fcsdat_transform = fcsdat;
    figure(6)
    
    linear_axis = [linspace(-100, 90, 10) logspace(2, 6, 41)];
    logicle_axis = M_logicle * logicleTransform(linear_axis,T_logicle,W_logicle,M_logicle,A_logicle);
    
    fig_handle = gca;
    fig_handle.XTick = logicle_axis;
    fig_handle.YTick = logicle_axis;
    
    x_temp = all_data.(mirs_of_interest).x_data_raw;
    y_temp = all_data.(mirs_of_interest).y_data_raw;
    
    [x_temp, sort_index] = sort(x_temp);
    y_temp = y_temp(sort_index);
    x_temp = x_temp(1:length(x_temp)-50);
    y_temp = y_temp(1:length(y_temp)-50);
    
    x_temp = M_logicle * logicleTransform(x_temp,T_logicle,W_logicle,M_logicle,A_logicle);
    y_temp = M_logicle * logicleTransform(y_temp,T_logicle,W_logicle,M_logicle,A_logicle);
    
    
    plot(x_temp,y_temp,'r.')
    hold on
        
    plot(P_2_transform,P_1_transform)


    %figure(7)
    %loglog(P_2,P_1)

    %xlim([1 1e8])
    %ylim([1 1e8])
    %hold on
end
figure(6)
hold off
%figure(7)
%hold off

grid on


