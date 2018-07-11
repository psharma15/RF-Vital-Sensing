% Plotting the frequency spectrum

function plotFreqSpectrum(data, fs)
% Input: 
%       data: Input signal
%       fs  : Sampling frequency in Hz

   %% Time specifications:
   dt = 1/fs;                     % seconds per sample
   stopTime = length(data)/fs;                  % seconds
   t = (0:dt:stopTime-dt)';
   nSample = size(t,1);
   %% Fourier Transform:
   X = fftshift(fft(data));
   %% Frequency specifications:
   dF = fs/nSample;                      % hertz
   f = -fs/2:dF:fs/2-dF;           % hertz
   %% Plot the spectrum:
   figure;
   plot(f,abs(X)/nSample);
   xlabel('Frequency (in hertz)');
   title('Magnitude Response');
end