% Replace ABC with your specific index number digits
index_number = '377';

% Construct the signal variable name based on your index number
signal_variable_name = ['signal' index_number];

% Load the corresponding signal from the mat file
load([signal_variable_name '.mat'], 'xn_test');

% Plot the signal
plot(xn_test);

% Add labels and title for clarity
xlabel('Time');
ylabel('Amplitude');
title('Plot of signal377');




