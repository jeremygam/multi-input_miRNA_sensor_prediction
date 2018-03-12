function [x_bin,y_bin,x_bin_raw,y_bin_raw] = binning(x_logicle,y_logicle,x_raw,y_raw,n)
% puts pre-sorted raw and transformed vectors into n equivalent bins

%initialize
y_bin = zeros(n-1,1);
y_bin_raw = zeros(n-1,1);
x_bin = zeros(n-1,1);
x_bin_raw = zeros(n-1,1);
x_min = min(x_logicle);
x_max = max(x_logicle);

%find bin edges of logicle transformed data and the corresponding edges for raw data
x_bin_edges = linspace(x_min, x_max, n+1)';
x_bin_edges_raw = interp1(unique(x_logicle),unique(x_raw),x_bin_edges,'pchip');

%calculate bins
bin_ind = discretize(x_logicle,x_bin_edges);     %invariant under transformation

%sort into bins
for jj = 1:n       %last point in the bin is the right edge so exclude it
    y_in_bin     = y_logicle(bin_ind == jj);
    y_in_bin_raw = y_raw(bin_ind == jj);
    
    %make sure y_in_bin is not empty (otherwise median returns NaN)
    if ~isempty(y_in_bin)
        y_bin(jj) = median(y_in_bin);
    end
    if ~isempty(y_in_bin_raw)
        y_bin_raw(jj) = median(y_in_bin_raw);
    end
    
    %take the center of each bin for the x value
    x_bin(jj) = mean([x_bin_edges(jj);x_bin_edges(jj+1)]);
    x_bin_raw(jj) = mean([x_bin_edges_raw(jj);x_bin_edges_raw(jj+1)]);
end

%drop bins that are empty
x_bin_raw = x_bin_raw(y_bin ~= 0);
y_bin_raw = y_bin_raw(y_bin ~= 0);
x_bin = x_bin(y_bin ~= 0);
y_bin = y_bin(y_bin ~= 0);
%drop bins that are empty
x_bin = x_bin(y_bin_raw ~= 0);
y_bin = y_bin(y_bin_raw ~= 0);
x_bin_raw = x_bin_raw(y_bin_raw ~= 0);
y_bin_raw = y_bin_raw(y_bin_raw ~= 0);

