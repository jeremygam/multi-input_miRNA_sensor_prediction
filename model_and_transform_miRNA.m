function output = model_and_transform_miRNA(param,x_transform,a,b,c,d,T_logicle,W_logicle,M_logicle,A_logicle)
% This function is the model we use for fitting logicle transformed data
% to. The expected value for y_transform is found by taking the inverse
% logicle transform of x_transform, solving for y using model solution,
% then taking logicle transform again.
% Written by Jeremy Gam (jgam@mit.edu)

    Km  = param(1);
    M = param(2);

    % Since input x data is logicle transformed, take inverse logicle
    x = x_transform./M_logicle;
    x = logicleInverseTransform(x,T_logicle,W_logicle,M_logicle,A_logicle);
    
    y = d/2*(sqrt((-a*c*x + b*M + Km).^2 + 4*a*c*x*Km) + a*c*x - b*M - Km);

    output = logicleTransform(y,T_logicle,W_logicle,M_logicle,A_logicle)*M_logicle;
    
end
