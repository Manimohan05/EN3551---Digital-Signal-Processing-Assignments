% Define the signal length and sampling rate
N = length(xn_test);  % 1793 samples
fs = 128;  % Hz

% Construct subsets
S1 = xn_test(1:128);
S2 = xn_test(1:256);
S3 = xn_test(1:512);
S4 = xn_test(1:1024);
S5 = xn_test(1:1792);

% Length of each subset
N1 = length(S1);
N2 = length(S2);
N3 = length(S3);
N4 = length(S4);
N5 = length(S5);

% Plot constructed subsets
plot(S1);
xlabel('Time');
ylabel('Amplitude');
title('Subset 1 (S1)');

plot(S2);
xlabel('Time');
ylabel('Amplitude');
title('Subset 2 (S2)');

plot(S3);
xlabel('Time');
ylabel('Amplitude');
title('Subset 3 (S3)');

plot(S4);
xlabel('Time');
ylabel('Amplitude');
title('Subset 4 (S4)');

plot(S5);
xlabel('Time');
ylabel('Amplitude');
title('Subset 5 (S5)');