% Function for DFT-based interpolation
dft_interpolate = @(x, K) ifft(upsample(fft(x), K + 1));

% Interpolate x2 with K = 1
K = 1;
interpolated_x2 = dft_interpolate(x2, K);

% Interpolate x3 with K = 2
K = 2;
interpolated_x3 = dft_interpolate(x3, K);

% Interpolate x4 with K = 3
K = 3;
interpolated_x4 = dft_interpolate(x4, K);
