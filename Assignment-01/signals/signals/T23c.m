% Compute the difference between interpolated_x4 and x
difference_x4 = interpolated_x4(1:length(x)) - x;

% Compute the 2-norm of the difference
norm_difference_x4 = norm(difference_x4);

% Plot the first 50 samples of both signals
figure;
subplot(2, 1, 1); plot(x(1:50)); title('Original Signal x (First 50 Samples)');
subplot(2, 1, 2); plot(interpolated_x4(1:50)); title('Interpolated Signal x4 (First 50 Samples)');
