clc
clear 
close all
addpath('fdata');
[FILENAME,PATHNAME]=uigetfile('*.jpg','Select the Underwater Image');
FilePath=strcat(PATHNAME,FILENAME);
disp('The Image File Location is');
disp(FilePath);
DataArray=imread(FilePath);
figure,imshow(DataArray);
title('Original image');
input_img=DataArray;
% Convert the input image to double precision for processing to avoid loss
% of pixels 
input_img=im2double(input_img);
enhancedImage = Color_Transform(input_img);
figure,imshow(enhancedImage);
title('Color Transform image');


% Define the size of the local window for calculating statistics
window_size = 15;
% Compute the local mean and standard deviation for each channel
local_mean = imfilter(input_img, ones(window_size) / window_size^2, 'replicate');
local_mean_sq = imfilter(input_img.^2, ones(window_size) / window_size^2, 'replicate');
local_std = sqrt(local_mean_sq - local_mean.^2);
% Calculate the detail image by subtracting the local mean
detail_image = input_img - local_mean;
figure,imshow(detail_image);
title('Detailed image');
% Compute the attenuation image
attenuation_image = 1 - exp(-local_std);
figure,imshow(attenuation_image);
title('Attenuation image');
size(attenuation_image)
% Perform color transfer using the attenuation image
color_transferred_img = zeros(size(input_img));
for channel = 1:3
    color_transferred_img(:,:,channel) = double((enhancedImage(:,:,channel))) + double(attenuation_image(:,:,channel)) + double(detail_image(:,:,channel));
end 
figure,imshow(uint8(color_transferred_img));
title('Color Corrected Image');

% Read the input underwater image
input_image = uint8(color_transferred_img);

% Convert the input image to CIELAB color space

% Convert RGB to Lab
R=color_transferred_img(:,:,1);
G=color_transferred_img(:,:,2);
B=color_transferred_img(:,:,3);
if max(max(R)) > 1.0 || max(max(G)) > 1.0 || max(max(B)) > 1.0
  R = double(R) / 255;
  G = double(G) / 255;
  B = double(B) / 255;
end
% Set a threshold
T = 0.008856;
[M, N] = size(R);
s = M * N;
RGB = [reshape(R,1,s); reshape(G,1,s); reshape(B,1,s)];
% RGB to XYZ
MAT = [0.412453 0.357580 0.180423;
       0.212671 0.715160 0.072169;
       0.019334 0.119193 0.950227];
XYZ = MAT * RGB;
% Normalize for D65 white point
X = XYZ(1,:) / 0.950456;
Y = XYZ(2,:);
Z = XYZ(3,:) / 1.088754;
XT = X > T;
YT = Y > T;
ZT = Z > T;
Y3 = Y.^(1/3); 
fX = XT .* X.^(1/3) + (~XT) .* (7.787 .* X + 16/116);
fY = YT .* Y3 + (~YT) .* (7.787 .* Y + 16/116);
fZ = ZT .* Z.^(1/3) + (~ZT) .* (7.787 .* Z + 16/116);
L = reshape(YT .* (116 * Y3 - 16.0) + (~YT) .* (903.3 * Y), M, N);
a = reshape(500 * (fX - fY), M, N);
b = reshape(200 * (fY - fZ), M, N);
% Extract the L, a, and b channels
L_channel = L;
a_channel = a;
b_channel = b;
% Define local block size and compute integral and squared integral maps
% Local contrast enhancement of the L channel

block_size = 5; % Adjust block size as needed
integral_map = cumsum(cumsum(L_channel, 1), 2);
squared_integral_map = cumsum(cumsum(L_channel.^2, 1), 2);
% Compute local mean and variance using integral maps
local_mean = integral_map(block_size:end, block_size:end) - integral_map(1:end-block_size+1, block_size:end) ...
    - integral_map(block_size:end, 1:end-block_size+1) + integral_map(1:end-block_size+1, 1:end-block_size+1);
local_variance = squared_integral_map(block_size:end, block_size:end) - squared_integral_map(1:end-block_size+1, block_size:end) ...
    - squared_integral_map(block_size:end, 1:end-block_size+1) + squared_integral_map(1:end-block_size+1, 1:end-block_size+1);
% Compute enhancement control factor alpha
alpha = local_variance ./ (local_variance + mean(local_variance(:))); % Adjust as needed
% Expand local_mean to match the size of L_channel
local_mean_expanded = padarray(local_mean, [block_size-1, block_size-1], 'pre');
alpha_expanded = padarray(alpha, [block_size-1, block_size-1], 'pre');
size(local_mean_expanded)
size(L_channel)
% Apply local contrast enhancement
enhanced_L_channel = local_mean_expanded + alpha_expanded .* (L_channel - local_mean_expanded);

% Color balance of a and b channels
% Compute mean intensity values of a and b channels
mean_a = mean(a_channel(:));
mean_b = mean(b_channel(:));
% Perform compensation for color balance
balanced_a_channel = a_channel + (mean_b - mean_a) / (mean_a + mean_b) * a_channel;
balanced_b_channel = b_channel + (mean_a - mean_b) / (mean_a + mean_b) * b_channel;

% Combine the enhanced L channel and balanced a and b channels
enhanced_lab_image = cat(3, enhanced_L_channel, balanced_a_channel, balanced_b_channel);

% Convert the enhanced Lab image back to RGB color space
enhanced_rgb_image = LACC(DataArray); 
Result = (LACE(enhanced_rgb_image));
% Display the enhanced image
figure,
imshow(uint8(Result));
title('Enhanced Underwater Image');

% Performance Analysis
uiqm = UIQM(uint8(Result));
disp('UIQM')
disp(uiqm)
uciqe = UCIQE(uint8(Result));
disp('UCIQE')
disp(uciqe)
