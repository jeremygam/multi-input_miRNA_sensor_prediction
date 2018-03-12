function createJCGate_JG(fileName)
% Written by Breanna DiAndreth:
% creates a file with gates for each of the Just Cell gates:
%       gate1: SSC-A vs. FSC-A
%       gate2: FSC-H vs. FSC-W
%       gate3: SSC-H vs. SSC-W
%       gate4: Time vs. Blue
%       gate5: Time vs Yellow
%       gate6: Time vs Red

%can decide if you want to use time gating as well as FSC/SSC.
%modified by Jeremy Gam (jgam@mit.edu) to include time gating
%if user specifies and also to use logicals when gating rather
%than indeces

baseName=fileName(1:end-4);
gateFileName = baseName; %modified by JG. previously gateFileName=['JCGate_' baseName];
alternate_file_name = ['JCGate_' baseName]; %added for backwards compatibility with Bre's previous gate files

all_gates = struct;

if ~( exist([gateFileName '.mat'],'file')==2 || exist([alternate_file_name '.mat'],'file')==2 ) %only make a new gate file if one doesn't exist yet
    
    loop_count = 1;
    
    % read in fcs data in order to gate
    [fcsdat, fcshdr] = fca_readfcs(fileName);
    
    logicals_cumul = ones(size(fcsdat,1), 1);
    
    % generate an array of collected fcs parameters
    par_names = cell(length(fcshdr.par)+1,1);
    for i = 1:length(fcshdr.par)
        par_names{i} = fcshdr.par(i).name;
    end
    par_names{i+1} = 'final';

    % input how you want to gate the first gate
    disp('Valid parameter names: ')
    disp(par_names)
    x_axis_par = input('x-axis parameter for gating: ','s');
    y_axis_par = input('y-axis parameter for gating: ','s');
    gate_mode = input('Gate mode (p=polygon, t=text): ','s');
    plot_scale = 'na';
    if gate_mode == 'p'
        disp('Valid scale names: ')
        disp({'logicle';'loglog';'semilogx';'semilogy';'linear'})
        plot_scale = input('Scale for plot: ','s');
    end
    
    % determine if loop should continue
    not_end = ~strcmp(x_axis_par,'final') & ~strcmp(y_axis_par,'final') & ~strcmp(gate_mode,'final');
    valid_par = ismember(x_axis_par,par_names) & ismember(x_axis_par,par_names) & (strcmp(gate_mode,'p') | strcmp(gate_mode, 't'));
    
    % loop and generate gates until user indicates a 'final'
    while (not_end)
        if (valid_par)
            % read in data
            x_axis_chan = getChannel(fcshdr,x_axis_par);
            y_axis_chan = getChannel(fcshdr,y_axis_par);
            x_axis_data = fcsdat(:,x_axis_chan);
            y_axis_data = fcsdat(:,y_axis_chan);
            
            % remove data that was gated out in previous iterations
            x_axis_data(~logicals_cumul)= NaN;
            y_axis_data(~logicals_cumul)= NaN;
            
            % use polygon to gate
            if gate_mode == 'p'
                figure
                xlabel(x_axis_par)
                ylabel(y_axis_par)
                hold on
                [logicals_temp, gate_temp] = gatePolygon(x_axis_data,y_axis_data,plot_scale);
                logicals_cumul = logicals_cumul .* logicals_temp;
                
                gate_name = ['gate' num2str(loop_count)];
                all_gates.(gate_name).logicals = logicals_cumul;
                all_gates.(gate_name).x_axis_par = x_axis_par;
                all_gates.(gate_name).y_axis_par = y_axis_par;
                all_gates.(gate_name).gate_polygon = gate_temp;
            end
            
            % use text to gate
            if gate_mode == 't'
                gate_string = input('Input your gate matrix in [nx2] format: ','s');
                gate_temp = eval(gate_string);
                logicals_temp = inpolygon(x_axis_data,y_axis_data,gate_temp(:,1),gate_temp(:,2));
                logicals_cumul = logicals_cumul .* logicals_temp;
                
                gate_name = ['gate' num2str(loop_count)];
                all_gates.(gate_name).logicals = logicals_cumul;
                all_gates.(gate_name).x_axis_par = x_axis_par;
                all_gates.(gate_name).y_axis_par = y_axis_par;
                all_gates.(gate_name).gate_polygon = gate_temp;
            end
            
            loop_count = loop_count + 1;
            
        end
        
        % set up next iteration                
        % input how you want to gate the next gate
        disp('Valid parameter names: ')
        disp(par_names)
        x_axis_par = input('x-axis parameter for gating: ','s');
        y_axis_par = input('y-axis parameter for gating: ','s');
        gate_mode = input('Gate mode (p=polygon, t=text): ','s');
        plot_scale = 'na';
        if gate_mode == 'p'
            disp('Valid scale names: ')
            disp({'logicle';'loglog';'semilogx';'semilogy';'linear'})
            plot_scale = input('Scale for plot: ','s');
        end
        
        % determine if loop should continue
        not_end = ~strcmp(x_axis_par,'final') & ~strcmp(y_axis_par,'final') & ~strcmp(gate_mode,'final');
        valid_par = ismember(x_axis_par,par_names) & ismember(x_axis_par,par_names) & (strcmp(gate_mode,'p') | strcmp(gate_mode, 't'));

    end
        

    save(gateFileName,'all_gates')
   
    
else
    %warnmsg=[gateFileName '.mat already exhists'];
    %warning(warnmsg)
end

hold off

end