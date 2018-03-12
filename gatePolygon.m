function [cells_within_gate, gate_vertices] = gatePolygon(xdat,ydat,axis_scales,position)
% gates a given population.  Uses the polygon vertices specified by
% 'position' or if this is not given, the user is asked to specify the polygon through UI
% Written by Breanna DiAndreth (bstillo@mit.edu), modified by Jeremy Gam (jgam@mit.edu) to use inpolygon

T_logicle = 262144;
M_logicle = 4.5;
r_logicle = -10;
W_logicle = (M_logicle-log10(T_logicle/abs(r_logicle)))/2;
A_logicle = 0;

if nargin<4
%     figure
    if strcmp(axis_scales,'loglog')
        loglog(xdat,ydat,'.','MarkerSize',2)
        set(gca, 'XScale', 'log','YScale', 'log')
    elseif strcmp(axis_scales,'semilogy')
        %dscatter(xdat,ydat,'MARKER','.','MSIZE',2)
        semilogy(xdat,ydat,'.','MarkerSize',2)
        set(gca, 'YScale', 'log')
    elseif strcmp(axis_scales,'semilogx')
        semilogx(xdat,ydat,'.','MarkerSize',2)
        set(gca, 'XScale', 'log')
    elseif strcmp(axis_scales,'logicle')
        xdat = logicleTransform(xdat,T_logicle,W_logicle,M_logicle,A_logicle);
        ydat = logicleTransform(ydat,T_logicle,W_logicle,M_logicle,A_logicle);
        plot(xdat,ydat,'.','MarkerSize',2)
    else
        plot(xdat,ydat,'.','MarkerSize',2)
    end
    %hold on

    h = impoly; %modified by JG to use impoly to determine if event is within gate.
    position = wait(h);

    % make the gate vertices matrix a closed shape
    position = [position; position(1,:)];
    
    % inverse logicle transform if necessary
    if strcmp(axis_scales,'logicle')
        position(:,1) = logicleInverseTransform(position(:,1),T_logicle,W_logicle,M_logicle,A_logicle);
        position(:,2) = logicleInverseTransform(position(:,2),T_logicle,W_logicle,M_logicle,A_logicle);
    end
end

gate_vertices = position;


% determine if in or out of polygon
cells_within_gate = inpolygon(xdat,ydat,position(:,1),position(:,2));


end