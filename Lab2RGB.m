function [R,G,B] = Lab2RGB(L, a, b)



if nargin == 1
  b = L(:,:,3);
  a = L(:,:,2);
  L = L(:,:,1);
end

% Thresholds
T1 = 0.008856;
T2 = 0.206893;

[M, N] = size(L);
s = M * N;
L = reshape(L, 1, s);
a = reshape(a, 1, s);
b = reshape(b, 1, s);

% Compute Y
fY = ((L + 16) / 116) .^ 3;
YT = fY > T1;
fY = (~YT) .* (L / 903.3) + YT .* fY;
Y = fY;

% Alter fY slightly for further calculations
fY = YT .* (fY .^ (1/3)) + (~YT) .* (7.787 .* fY + 16/116);

% Compute X
fX = a / 500 + fY;
XT = fX > T2;
X = (XT .* (fX .^ 3) + (~XT) .* ((fX - 16/116) / 7.787));

% Compute Z
fZ = fY - b / 200;
ZT = fZ > T2;
Z = (ZT .* (fZ .^ 3) + (~ZT) .* ((fZ - 16/116) / 7.787));

% Normalize for D65 white point
X = X * 0.950456;
Z = Z * 1.088754;

% XYZ to RGB
MAT = [ 3.240479 -1.537150 -0.498535;
       -0.969256  1.875992  0.041556;
        0.055648 -0.204043  1.057311];

RGB = max(min(MAT * [X; Y; Z], 1), 0);

R = reshape(RGB(1,:), M, N);
G = reshape(RGB(2,:), M, N);
B = reshape(RGB(3,:), M, N); 

if nargout < 2
  R = uint8(round(cat(3,R,G,B) * 255));
end
