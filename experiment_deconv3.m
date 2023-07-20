rng('default')
x = double(rgb2gray(imread('data/QR_code.jpg')));
H = fspecial('average',3);
y = awgn(imfilter(x,H,'circular'),15,'measured');
snr_cnc1 = zeros(8,1);
snr_cnc2 = zeros(8,1);
snr_cnc3 = zeros(8,1);
snr_cnc4 = zeros(8,1);
snr_convex1 = zeros(8,1);
snr_convex2 = zeros(8,1);
snr_convex3 = zeros(8,1);
snr_convex4 = zeros(8,1);
lambdas = 80:-10:10;
tic
[xhat1, vhat1, res_norm_hist1] = srls_GMC_matrix(y, 'deblurring', lambdas(1), 'H', H, 'gamma', 0.8, 'acceleration', 'original', 'splitting', 'FB');
snr_cnc1(1) = snr(x,xhat1-x);
for i=2:8
    [xhat1, vhat1, res_norm_hist1] = srls_GMC_matrix(y, 'deblurring', lambdas(i), 'H', H, 'gamma', 0.8, 'acceleration', 'original', 'splitting', 'FB', 'xv0', [xhat1(:);vhat1(:)]);
    snr_cnc1(i) = snr(x,xhat1-x);
end
time1 = toc
tic
[xhat2, vhat2, res_norm_hist2] = srls_GMC_matrix(y, 'deblurring', lambdas(1), 'H', H, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FB');
snr_cnc2(1) = snr(x,xhat2-x);
for i=2:8
    [xhat2, vhat2, res_norm_hist2] = srls_GMC_matrix(y, 'deblurring', lambdas(i), 'H', H, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FB', 'xv0' , [xhat2(:);vhat2(:)]);
    snr_cnc2(i) = snr(x,xhat2-x);
end
time2 = toc
tic
[xhat3, vhat3, res_norm_hist3] = srls_GMC_matrix(y, 'deblurring', lambdas(1), 'H', H, 'gamma', 0.8, 'acceleration', 'original', 'splitting', 'FBF');
snr_cnc3(1) = snr(x,xhat3-x);
for i=2:8
    [xhat3, vhat3, res_norm_hist3] = srls_GMC_matrix(y, 'deblurring', lambdas(i), 'H', H, 'gamma', 0.8, 'acceleration', 'original', 'splitting', 'FBF', 'xv0', [xhat3(:);vhat3(:)]);
    snr_cnc3(i) = snr(x,xhat3-x);
end
time3 = toc
tic
[xhat4, vhat4, res_norm_hist4] = srls_GMC_matrix(y, 'deblurring', lambdas(1), 'H', H, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FBF');
snr_cnc4(1) = snr(x,xhat4-x);
for i=2:8
    [xhat4, vhat4, res_norm_hist4] = srls_GMC_matrix(y, 'deblurring', lambdas(i), 'H', H, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FBF', 'xv0', [xhat4(:);vhat4(:)]);
    snr_cnc4(i) = snr(x,xhat4-x);
end
time4 = toc
tic
[xhat5, vhat5, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(1), 'H', H, 'gamma', 0, 'acceleration', 'original', 'splitting', 'FB');
snr_convex1(1) = snr(x,xhat5-x);
for i=2:8
    [xhat5, vhat5, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(i), 'H', H, 'gamma', 0, 'acceleration', 'original', 'splitting', 'FB', 'xv0', [xhat5(:);vhat5(:)]);
    snr_convex1(i) = snr(x,xhat5-x);
end
time5 = toc
tic
[xhat6, vhat6, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(1), 'H', H, 'gamma', 0, 'acceleration', 'aa2', 'splitting', 'FB');
snr_convex2(1) = snr(x,xhat6-x);
for i=2:8
    [xhat6, vhat6, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(i), 'H', H, 'gamma', 0, 'acceleration', 'aa2', 'splitting', 'FB', 'xv0', [xhat6(:);vhat6(:)]);
    snr_convex2(i) = snr(x,xhat6-x);
end
time6 = toc
tic
[xhat7, vhat7, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(1), 'H', H, 'gamma', 0, 'acceleration', 'original', 'splitting', 'FBF');
snr_convex3(1) = snr(x,xhat7-x);
for i=2:8
    [xhat7, vhat7, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(i), 'H', H, 'gamma', 0, 'acceleration', 'original', 'splitting', 'FBF', 'xv0', [xhat7(:);vhat7(:)]);
    snr_convex3(i) = snr(x,xhat7-x);
end
time7 = toc
tic
[xhat8, vhat8, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(1), 'H', H, 'gamma', 0, 'acceleration', 'aa2', 'splitting', 'FBF');
snr_convex4(1) = snr(x,xhat8-x);
for i=2:8
    [xhat8, vhat8, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(i), 'H', H, 'gamma', 0, 'acceleration', 'aa2', 'splitting', 'FBF', 'xv0', [xhat8(:);vhat8(:)]);
    snr_convex4(i) = snr(x,xhat8-x);
end
time8 = toc

[xhat1, ~, ~ ] = srls_GMC_matrix(y, 'deblurring', 30, 'H', H, 'gamma', 0, 'acceleration', 'aa2', 'splitting', 'FBF');
[xhat2, ~, ~ ] = srls_GMC_matrix(y, 'deblurring', 60, 'H', H, 'gamma', 0.8, 'acceleration', 'aa2', 'splitting', 'FBF'); 
save("results/deconv3.mat","res_norm_hist1","res_norm_hist2","res_norm_hist3","res_norm_hist4","x","y","xhat1","xhat2",...
    "snr_convex1","snr_convex2","snr_convex3","snr_convex4","snr_cnc1","snr_cnc2","snr_cnc3","snr_cnc4",...
    "time1","time2","time3","time4","time5","time6","time7","time8","lambdas")

