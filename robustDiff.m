function y_dot = robustDiff(y, dt, N)
%ROBUSTDIFF differentiates using smooth noise-robust differentiation formula
%
%   y_dot = robustDiff(y, dt, N)
%
%% Inputs
% y  - signal/vector to differentiate
% dt - time or distance between points in y
% N  - Number of points to use in differentiation.  This value must be
%      positive odd integer greater than or equal 5.
%
%% Outputs
% y_dot - differentiated signal/vector
%
%% Description
% robustDiff differentiates a signal/vector numerically using N
% points.  Both future information and past information are used to
% calculate the derivative.  In signal processing, this is called non-causal.
% The larger the value of N, the more high frequency noise is suppressed
% unlike Savitsky-Golay filters and Central Difference methods (see references).  
% Note that the derivative is not estimated at the edges of y.  This means that
% (N-1)/2 points at the beginning and end of y_dot are NaN.  See the example.
%
%% Example
%   dt = 0.001; % sampling rate of 1000Hz 
%   t = 0:dt:3; % sec
%   noiseFrequency = 400; % Hz
%   noise = 10*rand(size(t)); % Noise is 10% of signal
%   frequency = 1; %Hz
%   y = 100*sin(2*pi*frequency*t) + noise;
%   N = 21; % Number of points to use to estimate derivative
%   y_dot_estimate = robustDiff(y, dt, N);
%   y_dot_actual = 100*2*pi*frequency*cos(2*pi*frequency*t);
%   subplot(211);
%   plot(t, y);
%   title('y vs. t');
%   subplot(212);
%   plot(t, y_dot_actual, 'DisplayName', 'y''_{actual} of sin(t)'); 
%   hold('all')
%   plot(t, y_dot_estimate, 'DisplayName', 'y''_{estimate} of sin(t) + noise');
%   legend('show');
%   hold('off');
%   disp(['Beginning (N-1)/2 points and ending (N-1)/2 points of ' ...
%         'y_dot_estimate are NaN']);
%   y_dot_estimate(1:(N-1)/2)
%   y_dot_estimate(end-(N-1)/2+1:end)
%
%% References
% This function is based on the formulas by Pavel Holoborodko from his
% website: http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
% A big thanks is due to Pavel Holoborodko for developing these formulas.
%

%   Jason H. Nicholson jashale@yahoo.com
%   $Revision: 1.1 $  $Date: 2014/03/02 16:15:00 $

%% Input Checking
narginchk(3,3);

% N must be odd and greater than 5
if rem(N,2)==1 && N>=5
    % Do nothing
else
    error('N must be postive, odd integer greater than 5.')
end % end if, check N is odd and greater than 5

yLength = length(y);
% y must be have at least N points to estimate the derivative
if yLength >= N
    % Do nothing
else
    error('y must be greater or equal to N.');
end % check yLength is greater or equal N

%% Calculate Coefficients
% Equation for coefficients
% $$\displaystyle {c_k = \frac{1}{2^{2m+1}}\left[{2m\choose m-k+1}-{2m\choose m-k-1}\right]},\quad \displaystyle{m=\frac{N-3}{2}},\quad M=\frac{N-1}{2}$
%
% See reference for more information.
m = (N-3)/2;
M = (N-1)/2;
k = M:-1:1;

% Note that dividing by 2^(2*m+1) should be a bitshift but I don't know how
% to do this in matlab for type double
c_k = (binomialCoefficient(2*m, m-k+1) - binomialCoefficient(2*m, m-k-1))/2^(2*m+1);

%% Calculate coefficients for filter function
b = [c_k 0 -c_k(end:-1:1)];

%% Filter y
% The filter command only computes with past and present information so a
% shift of elements will be required after this step.
y_dot_intermediate = filter(b,1,y(:));


% divide by dt.  Discard first N-1 elements of y_dot_intermediate
y_dot_intermediate = y_dot_intermediate(N:end)/dt;

% shift elements, replace beginning and ending elements, and reshape to the
% same size as y
y_dot = reshape([NaN(M,1); 
                 y_dot_intermediate; 
                 NaN(M,1)],size(y));


end %end function, robustDiff



%% binomialCoefficient
function coefficients = binomialCoefficient(n, k)
% calculates the binomial coeffiecents given k which is a vector and n
% which is a scalar

% preallocated coefficents
coefficients = zeros(size(k));

% find values of k greater than or equal 0 and less than or equal n
index = k >= 0 & k <= n;

% Only calculate coefficients for logical true elements in index
coefficients(index) = arrayfun(@(x) nchoosek(n, x), k(index));
end
