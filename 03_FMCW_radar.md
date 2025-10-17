# FMCW Radar System

- `AIM`: To develop a MATLAB program that computes the range, azimuth angle, and elevation angle of multiple targets, using implementations of mathematical functions.
## Version-1 (Simplest)
```matlab
clear; clc;
% Multiple target positions (each column is a target: [x; y; z])
positions = [
    100, 150, 200;
    50,  60,  70;
    20,  25,  30
];
num_targets = size(positions, 2);
% Pre-allocate arrays for results
ranges = zeros(1, num_targets);
azimuths_deg = zeros(1, num_targets);
elevations_deg = zeros(1, num_targets);
% Constants
deg_per_rad = 180 / pi;
% Loop over targets
for k = 1:num_targets
    x = positions(1, k);
    y = positions(2, k);
    z = positions(3, k);
    % Compute range
    ranges(k) = sqrt(x^2 + y^2 + z^2);
    % Compute azimuth angle (atan2 returns radians)
    azimuth = atan2(y, x);
    % Compute elevation angle (asin returns radians)
    elevation = asin(z / ranges(k));
    % Convert radians to degrees
    azimuths_deg(k) = azimuth * deg_per_rad;
    elevations_deg(k) = elevation * deg_per_rad;
end
% Display results
for k = 1:num_targets
    fprintf('Target %d: Range = %.2f m, Azimuth = %.2f°, Elevation = %.2f°\n', ...
        k, ranges(k), azimuths_deg(k), elevations_deg(k));
end
```

## Version-2 (Maam's implmentation)
```matlab
clear; clc;
% Multiple target positions (each column is a target: [x; y; z])
positions = [
100, 150, 200;
50, 60, 70;
20, 25, 30
];
num_targets = size(positions, 2);
% Pre-allocate arrays for results
ranges = zeros(1, num_targets);
azimuths_deg = zeros(1, num_targets);
elevations_deg = zeros(1, num_targets);
% Constants
PI = 3.141592653589793;
deg_per_rad = 180 / PI;

% --- Custom math functions ---
function s = my_sqrt(val)
if val &lt;= 0
s = 0;
return;
end
guess = val / 2;
for i = 1:10
guess = 0.5 * (guess + val / guess);
end
s = guess;
end
function a = my_atan(t)
PI = 3.141592653589793;
sign_t = 1;
if t &lt; 0
t = -t;
sign_t = -1;
end
if t &lt;= 1
t2 = t*t;
a = t - (t2*t)/3 + (t2*t2*t)/5 - (t2*t2*t2*t)/7;
a = sign_t * a;
else
a = sign_t * (PI/2 - my_atan(1/t));
end
end
function angle = my_atan2(y,x)
PI = 3.141592653589793;
if x &gt; 0
angle = my_atan(y/x);
elseif x &lt; 0 &amp;&amp; y &gt;= 0
angle = my_atan(y/x) + PI;
elseif x &lt; 0 &amp;&amp; y &lt; 0
angle = my_atan(y/x) - PI;
elseif x == 0 &amp;&amp; y &gt; 0
angle = PI/2;
elseif x == 0 &amp;&amp; y &lt; 0
angle = -PI/2;
else
angle = 0;
end
end
function angle = my_asin(x)
x3 = x*x*x;
x5 = x3*x*x;
angle = x + x3/6 + (3*x5)/40;
end
% --- Loop over targets ---
for k = 1:num_targets
x = positions(1,k);
y = positions(2,k);
z = positions(3,k);

ranges(k) = my_sqrt(x*x + y*y + z*z);
azimuth = my_atan2(y, x);
elevation = my_asin(z / ranges(k));
azimuths_deg(k) = azimuth * deg_per_rad;
elevations_deg(k) = elevation * deg_per_rad;
end
% Display results
for k = 1:num_targets
fprintf(&#39;Target %d: Range = %.2f m, Azimuth = %.2f°, Elevation = %.2f°\n&#39;, ...
k, ranges(k), azimuths_deg(k), elevations_deg(k));
end
```
## Version-3 (Mine)
```matlab
% MATLAB program to compute range, azimuth angle, and elevation angle
% for multiple targets given their Cartesian coordinates (x, y, z).
%
% Features:
% - Accepts user input for target coordinates or uses example data.
% - Computes:
%   - Range R = sqrt(x^2 + y^2 + z^2)
%   - Azimuth angle theta = atan2(y, x) in degrees
%   - Elevation angle phi = asin(z / R) in degrees
% - Supports multiple targets via vectorized operations.
% - Includes error checking for input validity.
% - Displays results in a formatted table.
%
% Notes:
% - The formula R = c * Tr / 2 is not used as Cartesian coordinates are provided.
%   (If time delay Tr is needed, it can be added as an alternative input.)
% - All angles are output in degrees for readability.
%
% Example usage:
% Run the script and enter coordinates when prompted, or use default example data.
% Clear workspace and command window
clear;
clc;
% Main script
fprintf('=== Target Parameter Calculator ===\n');
fprintf('This program computes range, azimuth, and elevation angles for multiple targets.\n\n');
% Ask user if they want to input data or use example
choice = input('Enter 1 to input target coordinates, 2 for example data: ');
if choice == 1
    % User input for number of targets
    num_targets = input('Enter the number of targets: ');
    while num_targets < 1 || ~isnumeric(num_targets) || mod(num_targets, 1) ~= 0
        fprintf('Error: Number of targets must be a positive integer.\n');
        num_targets = input('Enter the number of targets: ');
    end
    
    % Initialize coordinate arrays
    x = zeros(num_targets, 1);
    y = zeros(num_targets, 1);
    z = zeros(num_targets, 1);
    
    % Input coordinates for each target
    for i = 1:num_targets
        fprintf('Enter coordinates for target %d:\n', i);
        x(i) = input('x-coordinate (m): ');
        y(i) = input('y-coordinate (m): ');
        z(i) = input('z-coordinate (m): ');
        while ~isnumeric(x(i)) || ~isnumeric(y(i)) || ~isnumeric(z(i))
            fprintf('Error: Coordinates must be numeric.\n');
            x(i) = input('x-coordinate (m): ');
            y(i) = input('y-coordinate (m): ');
            z(i) = input('z-coordinate (m): ');
        end
    end
else
    % Use example data
    fprintf('Using example data:\n');
    x = [1; 4; 7]; % Example x-coordinates (m)
    y = [2; 5; 8]; % Example y-coordinates (m)
    z = [3; 6; 9]; % Example z-coordinates (m)
    num_targets = length(x);
    fprintf('Example targets:\n');
    for i = 1:num_targets
        fprintf('Target %d: (x=%.2f, y=%.2f, z=%.2f)\n', i, x(i), y(i), z(i));
    end
end
% Call function to compute parameters
[R, theta, phi] = compute_target_params(x, y, z);
% Display results in a formatted table
fprintf('\n=== Results ===\n');
fprintf('Target | Range (m) | Azimuth (deg) | Elevation (deg)\n');
fprintf('-------|-----------|---------------|----------------\n');
for i = 1:num_targets
    fprintf('%6d | %9.2f | %13.2f | %14.2f\n', i, R(i), theta(i), phi(i));
end
% Optional: Plot targets in 3D (uncomment to enable)
figure;
scatter3(x, y, z, 'filled');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Target Positions in 3D Space');
grid on;
% Function to compute range, azimuth, and elevation
function [R, theta, phi] = compute_target_params(x, y, z)
    % Input validation
    if length(x) ~= length(y) || length(y) ~= length(z)
        error('Error: x, y, z must have the same length.');
    end
    if ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(z)
        error('Error: Coordinates must be numeric.');
    end
    
    % Compute range R = sqrt(x^2 + y^2 + z^2)
    R = sqrt(x.^2 + y.^2 + z.^2);
    
    % Check for zero range to avoid division by zero
    if any(R == 0)
        warning('Zero range detected. Setting azimuth and elevation to 0 for those targets.');
        theta = zeros(size(R));
        phi = zeros(size(R));
        return;
    end
    
    % Compute azimuth angle theta = atan2(y, x) in degrees
    theta = atan2(y, x) * (180 / pi);
    
    % Compute elevation angle phi = asin(z / R) in degrees
    phi = asin(z ./ R) * (180 / pi);
end
```
