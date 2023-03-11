x = double(rgb2gray(imread('QR_code.jpg')));
H = fspecial('average',5);
y = awgn(imfilter(x,H,'circular'),20,'measured');
snr_cnc1 = zeros(9,1);
snr_cnc2 = zeros(9,1);
snr_convex = zeros(9,1);
lambdas = 10:2.5:30;
%imshow(int8(min(max(x,0),255)))
%imshow(int8(min(max(y,0),255)))
tic
for i=1:9
    [xhat1, vhat1, res_norm_hist1] = srls_GMC_imaging(y, 'deblurring', lambdas(i), H=H, gamma=0.95, acceleration = 'aa2', splitting='FB');
    snr_cnc1(i) = snr(x,xhat1-x);
end
toc
tic
for i=1:9
    [xhat2, vhat2, res_norm_hist2] = srls_GMC_imaging(y, 'deblurring', lambdas(i), H=H, gamma=0.95, acceleration = 'aa2', splitting='FBF');
    snr_cnc2(i) = snr(x,xhat2-x);
end
toc
tic
for i=1:9
    [xhat3, vhat3, res_norm_hist3] = srls_GMC_imaging(y, 'deblurring', lambdas(i), H=H, gamma=0, acceleration = 'aa2', splitting='FB');
    snr_convex(i) = snr(x,xhat3-x);
end
toc
hold on
plot(lambdas,snr_cnc1)
plot(lambdas,snr_cnc2)
plot(lambdas,snr_convex)
legend('CNC1','CNC2','Location','east')
xlabel('lambda')
ylabel('SNR')
xlim([10 30])
hold off

x = checkerboard(32)*255;
mask = binornd(1,0.5,size(x));
y = awgn(x,10,'measured');
snr_cnc = zeros(10,1);
snr_convex = zeros(10,1);
lambdas = 10:10:100;
tic
for i=1:10
    [xhat1, vhat1, res_norm_hist1] = srls_GMC_imaging(y.*mask, 'inpainting', lambdas(i), gamma=0.95, acceleration = 'aa2',mask=mask, splitting='FBF');
    snr_cnc(i) = snr(x,xhat1-x);
end
toc
tic
for i=1:10
    [xhat2, vhat2, res_norm_hist2] = srls_GMC_imaging(y.*mask, 'inpainting', lambdas(i), gamma=0.95, acceleration = 'aa2',mask=mask, splitting='FB');
    snr_convex(i) = snr(x,xhat2-x);
end
toc
hold on
plot(lambdas,snr_cnc)
plot(lambdas,snr_convex)
legend('CNC','Convex','Location','east')
xlabel('lambda')
ylabel('SNR')
xlim([10,100])
hold off

x = checkerboard(32)*255;
mask = binornd(1,0.5,size(x));
y = awgn(x,10,'measured');
snr_cnc = zeros(10,1);
snr_convex = zeros(10,1);
lambdas = 1000:1000:10000;
tic
for i=1:10
    [xhat1, vhat1, res_norm_hist1] = srls_GMC_imaging(y.*mask, 'matrix completion', lambdas(i), gamma=0.95, acceleration = 'aa2',mask=mask, splitting='FBF');
    snr_cnc(i) = snr(x,xhat1-x);
end
toc
for i=1:10
    [xhat2, vhat2, res_norm_hist2] = srls_GMC_imaging(y.*mask, 'matrix completion', lambdas(i), gamma=0.95, acceleration = 'aa2',mask=mask, splitting='FB');
    snr_convex(i) = snr(x,xhat2-x);
end
hold on
plot(lambdas,snr_cnc)
plot(lambdas,snr_convex)
legend('CNC','Convex','Location','east')
xlabel('lambda')
ylabel('SNR')
xlim([1000,10000])
hold off