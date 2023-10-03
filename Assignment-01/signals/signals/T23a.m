% Compute the difference between interpolated_x2 and x
difference_x2 = interpolated_x2(1:length(x)) - x;

% Compute the 2-norm of the difference
norm_difference_x2 = norm(difference_x2);

% Plot the first 50 samples of both signals
figure;
subplot(2, 1, 1); plot(x(1:50)); title('Original Signal x (First 50 Samples)');
subplot(2, 1, 2); plot(interpolated_x2(1:50)); title('Interpolated Signal x2 (First 50 Samples)');
