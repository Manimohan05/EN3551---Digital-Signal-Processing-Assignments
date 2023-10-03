% Load the 'Hallelujah' music clip
load handel;

% Define the signal length
N = 20000;  % First 20,000 samples

% Create signals x, x2, x3, and x4
x = y(1:N);
x2 = x(1:2:N);  % Odd
x3 = x(1:3:N);  % Odd
x4 = x(1:4:N);  % Even
