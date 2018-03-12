function fcsdatGated = applyJCGate_JG(sampleFile, gateFile) %modified by JG for time gating.
%   Function to apply a flow cytometry gate to a fcs file. If
%   the gate file does not exist yet, runs another function 
%   (createJCgate_JG) for the user to generate it.
% 
%   Example:
%       fcsdat_gated_by_no_DNA = applyJCGate_JG(no_dna_control_filename,no_dna_control_filename);
%
%   Written by Breanna DiAndreth, edited by Jeremy Gam
%   bstillo@mit.edu, jgam@mit.edu
%   Last Updated: 2015-12-15;
%

suffix = gateFile(end-2:end);
if strcmp(suffix, 'fcs')
    createJCGate_JG(gateFile); %createJCGate checks to see if there is a mat file and doesn't run if there is one
    baseName=gateFile(1:end-4);
    gateFileName=[baseName '.mat'];% modified by JG. previously was: gateFileName=['JCGate_' baseName '.mat'];

    if exist(['JCGate_' baseName '.mat'],'file')==2 %for backwards compatibility with Bre gate files
        gateFileName = ['JCGate_' baseName '.mat'];
    end
    
else
    gateFileName=gateFile;
end

GF = load(gateFileName);

[fcsdat, fcshdr] = fca_readfcs(sampleFile);

logicals_cumul = ones(size(fcsdat,1),1);

for i = 1:length(fieldnames(GF.all_gates))
    gate_name = ['gate' num2str(i)];
    logicals_cumul = logical(logicals_cumul);
    
    % get data
    x_axis_par = GF.all_gates.(gate_name).x_axis_par;
    y_axis_par = GF.all_gates.(gate_name).y_axis_par;
    x_axis_chan = getChannel(fcshdr,x_axis_par);
    y_axis_chan = getChannel(fcshdr,y_axis_par);
    x_data_temp = fcsdat(:,x_axis_chan);
    y_data_temp = fcsdat(:,y_axis_chan);
    
    % apply previous gates
    x_data_temp(~logicals_cumul)= NaN;
    y_data_temp(~logicals_cumul)= NaN;
    
    [logicals_temp gate_temp] = gatePolygon(x_data_temp,y_data_temp,'na',GF.all_gates.(gate_name).gate_polygon);
    
    logicals_cumul = logicals_cumul .* logicals_temp;
    
end

logicals_cumul = logical(logicals_cumul);
fcsdatGated = fcsdat(logicals_cumul,:);

%{
ch_SSC_A=getChannel(fcshdr,'SSC-A');
ch_FSC_A=getChannel(fcshdr,'FSC-A');
ch_FSC_W=getChannel(fcshdr,'FSC-W');
ch_FSC_H=getChannel(fcshdr,'FSC-H');
ch_SSC_W=getChannel(fcshdr,'SSC-W');
ch_SSC_H=getChannel(fcshdr,'SSC-H');
ch_time=getChannel(fcshdr,'Time');
ch_blue=getChannel(fcshdr,'Cyan');
ch_yellow=getChannel(fcshdr,'FIT');
ch_red=getChannel(fcshdr,'Red');



[P1inds g1]=gatePolygon(fcsdat(:,ch_FSC_A),fcsdat(:,ch_SSC_A),'semilogy',GF.gate1);
[P2inds g2]=gatePolygon(fcsdat(:,ch_FSC_W),fcsdat(:,ch_FSC_H),'linear',GF.gate2);
[P3inds g3]=gatePolygon(fcsdat(:,ch_SSC_W),fcsdat(:,ch_SSC_H),'semilogy',GF.gate3);
JCInds = P1inds & P2inds & P3inds;

if isfield(GF,'gate4') && isfield(GF,'gate5') && isfield(GF,'gate6') % modified by JG. only does this if time is applied in the gate file
    [P4inds g4]=gatePolygon(fcsdat(:,ch_blue),fcsdat(:,ch_red),'loglog',GF.gate4);
    [P5inds g5]=gatePolygon(fcsdat(:,ch_time),fcsdat(:,ch_yellow),'semilogy',GF.gate5);
    [P6inds g6]=gatePolygon(fcsdat(:,ch_time),fcsdat(:,ch_red),'semilogy',GF.gate6);
    %JCInds = P1inds & P2inds & P3inds & P4inds & P5inds & P6inds;
    JCInds = P1inds & P2inds & P3inds & P4inds;

end

fcsdatGated=fcsdat(JCInds,:);

%}

end