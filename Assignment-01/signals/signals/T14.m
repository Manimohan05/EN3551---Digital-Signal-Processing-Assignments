% Define parameters
K = 128;  % Number of samples in each subset
L = 14;   % Number of subsets
N = L * K; % Total number of samples

% Initialize an array to store the average DFT
Avg_DFT = zeros(1, K);

% Partition the signal and compute the DFT for each subset
for i = 1:L
    subset = xn_test((i - 1) * K + 1 : i * K);
    subset_DFT = fft(subset,K);
    Avg_DFT = Avg_DFT + subset_DFT;
end

% Calculate the average DFT by dividing by L
Avg_DFT = Avg_DFT / L;

% Create a frequency axis
frequencies = linspace(-fs/2, fs/2, K);

% Plot the magnitude of the average DFT
figure('Position', [100, 100, 800, 300]);
plot(frequencies,abs(Avg_DFT));
title('Magnitude of Average DFT');

xlabel('Frequency (Hz)');
ylabel('Magnitude');