function ll = line_length(x)
% calculate the line length averaged by number of data points in x
y = x(1:end-1,:);
z = x(2:end,:);

ll = mean(abs(z-y),1);


end