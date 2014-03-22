function y_dot = robustDiffOneSide(y, dt, N)
% differentiates using smooth noise-robust one sided differentiation formula
%
%   y_dot = robustDiffOneSide(y, dt, N)
%
%% Inputs
% y  - signal/vector to differentiate
% dt - time or distance between points in y
% N  - Number of points to use in differentiation.  This value must be
%      positive integer greater than or equal 2.
%
%% Outputs
% y_dot - differentiated signal/vector
%
%% Description
% robustDiffOneSide differentiates a signal/vector numerically using N
% points before the current point.  Only past information is used to
% calculate the derivative.  In signal processing, this is called causal.
% The larger the value of N, the more high frequency noise is suppressed
% unlike Savitsky-Golay filters and Central Difference methods 
% (see references).  Note that the derivative is not
% estimated at the beginning of y.  This means that (N-1) points at the
% beginning y_dot are NaN.  See the example.
%
%% Example
%   sampleRate = 1000; % Hz
%   dt = 1/sampleRate;
%   t = 0:dt:1;
%   signalFrequency = 40; % Hz
%   noiseMeanFrequency = 350; % Hz
%   noiseMeanAmplitude = 0.2;
%   pureSignal = sin(2*pi*signalFrequency*t);
%   noisySignal = sin(2*pi*signalFrequency*t) + noiseMeanAmplitude*(0.75+0.5*rand(size(t))).*cos(2*pi*noiseMeanFrequency*(0.9+0.1*rand(size(t))).*t+0.3);
%   signal_dot = 2*pi*signalFrequency*cos(2*pi*signalFrequency*t);
%   N = 9;
%   signal_dot_estimate = robustDiffOneSide(noisySignal,dt,9);
%   subplot(211)
%   plot(t,pureSignal,t,noisySignal);
%   xlim([0 4/signalFrequency])
%   xlabel('time [sec]')
%   ylabel('Amplitude')
%   title('y')
%   legend('Pure Signal','Noisy Signal');
%   subplot(212)
%   plot(t, signal_dot, t, signal_dot_estimate9)
%   xlim([0 4/signalFrequency])
%   xlabel('time [sec]')
%   ylabel('Amplitude')
%   title('signal''')
%   legend('Pure Signal Derivative', ['N = ' num2str(N)])
%
%% References
% This function is based on the formulas by Pavel Holoborodko from his
% website: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
% A big thanks is due to Pavel Holoborodko for developing these formulas.
%

%   Jason H. Nicholson jashale@yahoo.com
%   $Revision: 1.0 $  $Date: 2014/03/22 13:24:24 $

%% Input Checking
narginchk(3,3);

% N must be greater than 2
if N>=2
    % Do nothing
else
    error('N must be an integer greater than 2.')
end % end if, check N is greater than 2

yLength = length(y);
% y must be have at least N points to estimate the derivative
if yLength >= N
    % Do nothing
else
    error('y must be greater or equal to N.');
end % check yLength is greater or equal N

%% Calculate Coefficients
c_k = coefficientsOneSide(N);

%% Filter y
y_dot = filter(c_k,2^(N-1)*dt,y(:));

y_dot = reshape([NaN(N,1); 
                 y_dot(N+1:end)],size(y));


end %end function, robustDiffOneSide


%% coefficientsOneSide
% This is the heart of the function.  if you trying to implement this in
% some other programming language, use the function below as a template.
% Both recursive formulas and non recursive formulas exist.
function coefficients = coefficientsOneSide(n)

coefficients = zeros(1,n+1);

%% Recursive Method to Calculate Coefficients
% Note the recurrence is an integer sequence that can be found in The
% On-line Encyclopedia of Integer Sequences under A008315. Link:
% http://oeis.org/A008315
 
coefficients(1) = 1;
for iRow = 2:n
    previousCoefficients = coefficients;
    for iColumn = 2:((iRow+1)/2)
        coefficients(iColumn) = previousCoefficients(iColumn-1) + previousCoefficients(iColumn);
    end
end

%% Non-Recursive Method for Calculating Coefficients
% k = 1:(n+1)/2;
% coefficients(1) = 1;
% coefficients(2:ceil((n+1)/2)) = gamma(n)./(gamma(k+1).*gamma(n-k)).*(n-2*k)./(n-k);

%% Post Processing
% Reflect coefficients about the center of the vector and multiply by -1
coefficients(ceil((n+1)/2+1):end) = -coefficients(floor((n+1)/2):-1:1);
end % end function, coefficientsOnside