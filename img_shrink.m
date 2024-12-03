clear;
clc;
close all;

% Anh dau vao
input_image = imread('lena.jpg');
if size(input_image, 3) == 3
    input_image = rgb2gray(input_image); % chuyen sang grayscale
end

[M, N] = size(input_image);

% Kich thuoc khoi (N1, N2)
block_size = [8, 8];
N1 = block_size(1);
N2 = block_size(2);

M_trimmed = floor(M / N1) * N1;
N_trimmed = floor(N / N2) * N2;
input_image_trimmed = input_image(1:M_trimmed, 1:N_trimmed);
imwrite(uint8(input_image_trimmed), 'input_image.jpeg');

% Chia anh thanh cac khoi
image_blocks = mat2cell(input_image_trimmed, N1 * ones(1, M_trimmed / N1), N2 * ones(1, N_trimmed / N2));

% Tinh DCT cho tung khoi
dct_blocks = cellfun(@dct2, image_blocks, 'UniformOutput', false);

% Ap dung bo loc thong thap trong mien DCT
L = 8; % Do dai bo loc
h = fir1(L-1, 0.8); % Bo loc thong thap voi tan so cat 0.8
h_2d = h' * h;

% TInh DCT cua bo loc
h_dct = dct2(h_2d, N1, N2);
filtered_blocks = cellfun(@(block) block .* h_dct, dct_blocks, 'UniformOutput', false);

% Giam kich thuoc
downsampled_blocks = cellfun(@(block) ...
    (block(1:ceil(N1/2), 1:ceil(N2/2)) - block(N1-ceil(N1/2)+1:end, N2-ceil(N2/2)+1:end))/sqrt(2), ...
    filtered_blocks, 'UniformOutput', false);

% Tinh IDCT 
reconstructed_blocks = cellfun(@idct2, downsampled_blocks, 'UniformOutput', false);
output_image = cell2mat(reconstructed_blocks);

% Chuan hoa anh dau ra
output_image = uint8(mat2gray(output_image) * 255);
figure;

subplot(1, 2, 1);
imshow(input_image_trimmed, []);
title('original image');

subplot(1, 2, 2);
imshow(output_image, []);
title('processed image');
imwrite(uint8(output_image), 'output_image.jpg');
