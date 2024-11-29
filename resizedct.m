clear;
clc;
close all;

% �?�?c ảnh đầu vào
input_image = imread('lenag.jpeg');
if size(input_image, 3) == 3
    input_image = rgb2gray(input_image); % Chuyển sang grayscale
end

[M, N] = size(input_image);

% Kích thước khối (N1, N2)
block_size = [8, 8];
N1 = block_size(1);
N2 = block_size(2);

% Cắt ảnh sao cho kích thước chia hết cho kích thước khối
M_trimmed = floor(M / N1) * N1;
N_trimmed = floor(N / N2) * N2;
input_image_trimmed = input_image(1:M_trimmed, 1:N_trimmed);
imwrite(uint8(input_image_trimmed), 'input_image.jpeg');

% Chia ảnh thành các khối không chồng lấn
image_blocks = mat2cell(input_image_trimmed, N1 * ones(1, M_trimmed / N1), N2 * ones(1, N_trimmed / N2));

% Tính DCT 2D cho từng khối
dct_blocks = cellfun(@dct2, image_blocks, 'UniformOutput', false);

% �?p dụng bộ l�?c thông thấp trong mi�?n DCT
L = 8; % �?ộ dài bộ l�?c
h = fir1(L-1, 0.8); % Bộ l�?c thông thấp với tần số cắt 0.8
h_2d = h' * h;

% Tính DCT của bộ l�?c
h_dct = dct2(h_2d, N1, N2);
filtered_blocks = cellfun(@(block) block .* h_dct, dct_blocks, 'UniformOutput', false);

% Giảm kích thước khối
downsampled_blocks = cellfun(@(block) ...
    (block(1:ceil(N1/2), 1:ceil(N2/2)) - block(N1-ceil(N1/2)+1:end, N2-ceil(N2/2)+1:end))/sqrt(2), ...
    filtered_blocks, 'UniformOutput', false);

% Tính IDCT 
reconstructed_blocks = cellfun(@idct2, downsampled_blocks, 'UniformOutput', false);
output_image = cell2mat(reconstructed_blocks);

% Chuẩn hóa ảnh đầu ra
output_image = uint8(mat2gray(output_image) * 255);
figure;

subplot(1, 2, 1);
imshow(input_image_trimmed, []);
title('original image');

subplot(1, 2, 2);
imshow(output_image, []);
title('processed image');
imwrite(uint8(output_image), 'output_image.jpg');