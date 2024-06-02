d1 = 64;
d2 = d1;
X = zeros(d1);
row1 = 27;
row2 = 37;
X(row1:row2,10:54) = 1;
X(10:54,row1:row2) = 1;
imagesc(X);
colorbar;

X2 = checkerboard(d1/8);
imagesc(X2);
colorbar;

imshow([X, X2]);

imshowpair(X, X2, 'montage')