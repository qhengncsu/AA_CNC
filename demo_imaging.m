x = double(rgb2gray(imread('QR_code.jpg')));
H = fspecial('average',5);
y = awgn(imfilter(x,H,'circular'),20,'measured');
snr_cnc = zeros(12,1);
snr_convex = zeros(12,1);
lambdas = 2.5:2.5:30;
[xhat1, vhat1, res_norm_hist1] = srls_GMC_imaging(y, 'deblurring', 5, H=H, gamma=0.98, acceleration = 'aa2');
[xhat2, vhat2, res_norm_hist2] = srls_GMC_imaging(y, 'deblurring', 5, H=H, gamma=0.98, acceleration = 'original');
for i=1:12
    [xhat1, vhat1, res_norm_hist1] = srls_GMC_imaging(y, 'deblurring', lambdas(i), H=H, gamma=0.98, acceleration = 'aa2');
    snr_cnc(i) = snr(x,xhat1-x);
end
for i=1:12
    [xhat2, vhat2, res_norm_hist2] = srls_GMC_imaging(y, 'deblurring', lambdas(i), H=H, gamma=0, acceleration = 'aa2');
    snr_convex(i) = snr(x,xhat2-x);
end
hold on
plot(lambdas,snr_cnc)
plot(lambdas,snr_convex)
legend('CNC','Convex','Location','east')
xlabel('lambda')
ylabel('SNR')
xlim([2.5 30])
hold off

x = double(rgb2gray(imread('QR_code.jpg')));
mask = binornd(1,0.8,size(x));
y = awgn(x,20,'measured');
snr_cnc = zeros(16,1);
snr_convex = zeros(16,1);
lambdas = 2.5:2.5:40;
for i=1:16
    [xhat1, vhat1, res_norm_hist1] = srls_GMC_imaging(y.*mask, 'inpainting', lambdas(i), gamma=0.98, acceleration = 'aa2',mask=mask);
    snr_cnc(i) = snr(x,xhat1-x);
end
for i=1:16
    [xhat2, vhat2, res_norm_hist2] = srls_GMC_imaging(y.*mask, 'inpainting', lambdas(i), gamma=0, acceleration = 'aa2',mask=mask);
    snr_convex(i) = snr(x,xhat2-x);
end
hold on
plot(lambdas,snr_cnc)
plot(lambdas,snr_convex)
legend('CNC','Convex','Location','east')
xlabel('lambda')
ylabel('SNR')
xlim([2.5 40])
hold off