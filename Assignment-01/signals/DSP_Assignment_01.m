%% EN3551: Digital Signal Processing
%  Assignment 1: Detecting Harmonics in Noisy Data and Signal Interpolation using DFT
%  Name        : Randika Perera
%  Index No.   : 200462U




%% Section 3.1: Harmonic Detection

clc;
clear;
close all;


% TASK 1 

fs       = 128;     % Sampling rate
t_start  = 0;       % Starting time of the sampling
t_end    = 14;      % Ending time of the sampling
t_vector = (t_start:1/fs:t_end);
N        = length(t_vector); % Number of Samples. 
% Should be equal to 1793 since 128 x 14 + 1 = 1793

signal_file = load('signal462.mat');   % Load my signal file
signal_data = signal_file.xn_test;


% TASK 2 

% Construct Subsets
S1 = signal_data(1:128);
S2 = signal_data(1:256);
S3 = signal_data(1:512);
S4 = signal_data(1:1024);
S5 = signal_data(1:1792);


% TASK 3

% Apply DFT to each Subset
DFT_S1 = fft(S1);
DFT_S2 = fft(S2);
DFT_S3 = fft(S3);
DFT_S4 = fft(S4);
DFT_S5 = fft(S5);

% Obtain Magnitude Spectra of the Resulting DFT Sequences
magnitude_DFT_S1 = abs(DFT_S1);
magnitude_DFT_S2 = abs(DFT_S2);
magnitude_DFT_S3 = abs(DFT_S3);
magnitude_DFT_S4 = abs(DFT_S4);
magnitude_DFT_S5 = abs(DFT_S5);

% Display the DFT Magnitude Spectra 
figure(1);
subplot(5,1,1);
stem(magnitude_DFT_S1);
title('Magnitude Spectrum of DFT of S1');
xlabel('Index k');
ylabel('Magnitude');
subplot(5,1,2);
stem(magnitude_DFT_S2);
title('Magnitude Spectrum of DFT of S2');
xlabel('Index k');
ylabel('Magnitude');
subplot(5,1,3);
stem(magnitude_DFT_S3);
title('Magnitude Spectrum of DFT of S3');
xlabel('Index k');
ylabel('Magnitude');
subplot(5,1,4);
stem(magnitude_DFT_S4);
title('Magnitude Spectrum of DFT of S4');
xlabel('Index k');
ylabel('Magnitude');
subplot(5,1,5);
stem(magnitude_DFT_S5);
title('Magnitude Spectrum of DFT of S5');
xlabel('Index k');
ylabel('Magnitude');


% TASK 4

K = 128; % Length of each subset
L = 14;  % Number of subsets

% Partition the input signal into L subsets
subsets = reshape(signal_data(1:K*L), K, L);

% Apply DFT to each subset
DFT_subsets = fft(subsets);

% Calculate the arithmetic mean of the L sets of DFT sequences
X_avg = mean(DFT_subsets, 2);

% Find the frequencies corresponding to the peaks of the plot
frequencies = (0:K-1)*(fs/K);

% Find the 8 highest peaks
[~, sorted_indices] = sort(abs(X_avg), 'descend');
highest_peak_indices = sorted_indices(1:8);

% Display the results
figure(2);
stem(frequencies, abs(X_avg));
title('Magnitude Spectrum after DFT Averaging (L=14)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
hold on;
plot(frequencies(highest_peak_indices), abs(X_avg(highest_peak_indices)), 'ro');



% TASK 5 

% Running same code as TASK 4 for different values of K and L until 
% smallest value of L such that the four peaks that correspond to the four
% harmonics not greater than 64 remain clearly visible.

L = 2;          % Number of subsets
K = floor(N/L); % Length of each subset

% Partition the input signal into L subsets
subsets = reshape(signal_data(1:K*L), K, L);

% Apply DFT to each subset
DFT_subsets = fft(subsets);

% Calculate the arithmetic mean of the L sets of DFT sequences
X_avg = mean(DFT_subsets, 2);

% Find the frequencies corresponding to the peaks of the plot
frequencies = (0:K-1)*(fs/K);

% Find the 8 highest peaks
[~, sorted_indices] = sort(abs(X_avg), 'descend');
highest_peak_indices = sorted_indices(1:8);

% Display the results
figure(3);
stem(frequencies, abs(X_avg));
title('Magnitude Spectrum after DFT Averaging (L=2)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
hold on;
plot(frequencies(highest_peak_indices), abs(X_avg(highest_peak_indices)), 'ro');

% It seems like the smallest value for L such that the four peaks that 
% correspond to the four harmonics not greater than 64 remain clearly 
% visible is equal to 2 for the given input signal.


% TASK 6

% Define values of K for the subplots
K_values = [64, 128, 256, 100, 135, 200];

% Create a figure with 6 subplots
figure(4);

for i = 1:length(K_values)
    K = K_values(i);
    
    % Partition the input signal into L subsets
    subsets = reshape(signal_data(1:K*L), K, L);
    
    % Apply DFT to each subset
    DFT_subsets = fft(subsets);
    
    % Calculate the arithmetic mean of the L sets of DFT sequences
    X_avg = mean(DFT_subsets, 2);
    
    % Find the frequencies corresponding to the peaks of the plot
    frequencies = (0:K-1)*(fs/K);
    
    % Find the 8 highest peaks
    [~, sorted_indices] = sort(abs(X_avg), 'descend');
    highest_peak_indices = sorted_indices(1:8);
    
    % Create subplots
    subplot(2, 3, i);
    stem(frequencies, abs(X_avg));
    title(['K = ' num2str(K)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    hold on;
    plot(frequencies(highest_peak_indices), abs(X_avg(highest_peak_indices)), 'ro');
end




%% Section 3.2: Interpolation

clc;
clear;


load handel;
N  = 20000;                     
x  = y(1:N);
x2 = x(1:2:N);
x3 = x(1:3:N);
x4 = x(1:4:N);




% TASK 3a - Interpolate x2 with K=1

K = 1;
DFT_x2 = fft(x2);
N = length(x2);

if mod(N,2)==0  
    N2 = N/2;
    DFT_x2_with_zero_padding = [DFT_x2(1:N2); DFT_x2(N2+1)/2; zeros((K*N)-1,1); DFT_x2(N2+1)/2; DFT_x2((N2+2):N)];
else
    N1 = (N+1)/2;
    DFT_x2_with_zero_padding = [DFT_x2(1:N1); zeros(K*N, 1); DFT_x2((N1+1):N)];
end

x2_Interpolated = (K+1)*ifft(DFT_x2_with_zero_padding);

x2_Interpolated_Extracted = x2_Interpolated(1:((K+1)*(N-1))+1);

Norm_Difference_x2 = norm(x2_Interpolated_Extracted) - norm(x2);
disp(Norm_Difference_x2);

figure(11);
subplot(2, 1, 1);
stem(x2(1:50), 'b', 'Marker', 'o');
title('x2 - Original');
subplot(2, 1, 2)
stem(x2_Interpolated_Extracted(1:50), 'r', 'Marker', 'x');
title('x2 - Interpolated with K=1');




% TASK 3b - Interpolate x3 with K=2

K = 2;
DFT_x3 = fft(x3);
N = length(x3);

if mod(N,2)==0  
    N2 = N/2;
    DFT_x3_with_zero_padding = [DFT_x3(1:N2); DFT_x3(N2+1)/2; zeros((K*N)-1,1); DFT_x3(N2+1)/2; DFT_x3((N2+2):N)];
else
    N1 = (N+1)/2;
    DFT_x3_with_zero_padding = [DFT_x3(1:N1); zeros(K*N, 1); DFT_x3((N1+1):N)];
end

x3_Interpolated = (K+1)*ifft(DFT_x3_with_zero_padding);

x3_Interpolated_Extracted = x3_Interpolated(1:((K+1)*(N-1))+1);

Norm_Difference_x3 = norm(x3_Interpolated_Extracted) - norm(x3);
disp(Norm_Difference_x3);

figure(12);
subplot(2, 1, 1);
stem(x3(1:50), 'b', 'Marker', 'o');
title('x3 - Original');
subplot(2, 1, 2)
stem(x3_Interpolated_Extracted(1:50), 'r', 'Marker', 'x');
title('x3 - Interpolated with K=2');




% TASK 3c - Interpolate x4 with K=3

K = 3;
DFT_x4 = fft(x4);
N = length(x4);

if mod(N,2)==0  
    N2 = N/2;
    DFT_x4_with_zero_padding = [DFT_x4(1:N2); DFT_x4(N2+1)/2; zeros((K*N)-1,1); DFT_x4(N2+1)/2; DFT_x4((N2+2):N)];
else
    N1 = (N+1)/2;
    DFT_x4_with_zero_padding = [DFT_x4(1:N1); zeros(K*N, 1); DFT_x4((N1+1):N)];
end

x4_Interpolated = (K+1)*ifft(DFT_x4_with_zero_padding);

x4_Interpolated_Extracted = x4_Interpolated(1:((K+1)*(N-1))+1);

Norm_Difference_x4 = norm(x4_Interpolated_Extracted) - norm(x4);
disp(Norm_Difference_x4);

figure(13);
subplot(2, 1, 1);
stem(x4(1:50), 'b', 'Marker', 'o');
title('x4 - Original');
subplot(2, 1, 2)
stem(x4_Interpolated_Extracted(1:50), 'r', 'Marker', 'x');
title('x4 - Interpolated with K=3');




%% END OF CODE
