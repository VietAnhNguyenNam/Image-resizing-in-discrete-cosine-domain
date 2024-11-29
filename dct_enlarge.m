img = imread("small.jpg");
img = rgb2gray(img);
subplot(1,2,1);
imshow(img)
axis on;
title("original image: 668x1000 pixels");

blockSize = 8;
scaleFactor = 2;
gaussianSigma = 2.5;
output = enlarge(img, blockSize, scaleFactor, gaussianSigma);
% output = uint8(rescale(output, 0, 255));
subplot(1,2,2);
imshow(output, [])
axis on;
title("processed image: 1336x2000 pixels");


function res = blockProcessing(block, N, scaleFactor, Hrm)
    [rows, cols] = size(block);

    if rows ~= N || cols ~= N
        Hrm = imresize(Hrm, [rows cols], 'bicubic');
    end
    
    X = dct2(block);

    Y = zeros(size(X)*scaleFactor);
    Xu = zeros(size(X));
    
    Xu(1:ceil(rows/scaleFactor),1:ceil(cols/scaleFactor)) = X(1:ceil(rows/scaleFactor),1:ceil(cols/scaleFactor)); 
    Xu = (Xu-rot90(Xu, 2))/2^0.5;
    
    Y(1:N,1:N) = Xu;

    res = idct2((Hrm).*Y);
end

function res = enlarge(image, blockSize, scaleFactor, gaussianSigma)
    %low pass filter
    n = blockSize*scaleFactor;
    i=(0:n-1)';
    j=0:n-1;
    tmp = cos((2*i*j+i)*pi/(2*n));
    L = n;
    [x, y] = meshgrid(-L/2:L/2-1, -L/2:L/2-1);
    h2d = exp(-(x.^2 + y.^2) / (2 * gaussianSigma^2));
    h2d(1:L/2,1:L/2) = h2d(L/2+1:end,L/2+1:end);
    h2d(L/2+1:end,:) = 0; h2d(:,L/2+1:end) = 0;
    h2d = h2d / sum(h2d(:));
    B = 2 * tmp;
    Hrm = B*h2d*B';

    % block processing
    res = zeros(size(image)*scaleFactor);
    [rowno, colno] = size(image);
    for i=1:blockSize:rowno-blockSize+1
        for j=1:blockSize:colno-blockSize+1
            blk = blockProcessing(image(i:i+blockSize-1,j:j+blockSize-1), blockSize, scaleFactor, Hrm);
            row = (i-1)*scaleFactor+1; col = (j-1)*scaleFactor+1;
            res(row:row+blockSize*scaleFactor-1,col:col+blockSize*scaleFactor-1) = blk;
        end
    end
end