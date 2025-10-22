function [out_values,car_labels] = common_average_reference(values,include,labels)

out_values = values - nanmean(values(:,include),2); 
car_labels = strcat(labels,'-CAR');

end