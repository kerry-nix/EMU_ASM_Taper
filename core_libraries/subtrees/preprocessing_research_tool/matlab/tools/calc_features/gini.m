function result = gini(data)
%     Calculate the Gini coefficient of a numpy array.
%     From https://github.com/oliviaguest/gini

% based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
% from: http://www.statsdirect.com/help/default.htm%nonparametric_methods/gini.htm
% all values are treated equally, arrays must be 1d
data = data(:);
if min(data,[],'all') < 0
    data = data - min(data,[],'all');  % values cannot be negative
end
data = data + 1e-8;  % values cannot be 0
data = sort(data);  % values must be sorted
n = length(data);  % number of array elements
index = [1:n]';  % index per array element
result = (sum((2 * index - n - 1) .* data)) / (n * sum(data));  % Gini coefficient