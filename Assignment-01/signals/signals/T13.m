% Compute DFT for each subset
DFT_S1 = fft(S1,N1);
DFT_S2 = fft(S2,N2);
DFT_S3 = fft(S3,N3);
DFT_S4 = fft(S4,N4);
DFT_S5 = fft(S5,N5);

% Compute the magnitude of DFT sequences
Mag_S1 = abs(DFT_S1);
Mag_S2 = abs(DFT_S2);
Mag_S3 = abs(DFT_S3);
Mag_S4 = abs(DFT_S4);
Mag_S5 = abs(DFT_S5);

% Plot the magnitude spectra for each subset
subplot(5, 1, 1); plot(Mag_S1); title('Magnitude Spectrum of S1');
subplot(5, 1, 2); plot(Mag_S2); title('Magnitude Spectrum of S2');
subplot(5, 1, 3); plot(Mag_S3); title('Magnitude Spectrum of S3');
subplot(5, 1, 4); plot(Mag_S4); title('Magnitude Spectrum of S4');
subplot(5, 1, 5); plot(Mag_S5); title('Magnitude Spectrum of S5');

% Adjust the plot layout for better visualization
sgtitle('Magnitude of DFT Sequences for Different Subset Sizes');
