rng('default')
x = double(rgb2gray(imread('QR_code.jpg')));
H = fspecial('average',3);
y = awgn(imfilter(x,H,'circular'),15,'measured');
snr_cnc1 = zeros(8,1);
snr_cnc2 = zeros(8,1);
snr_cnc3 = zeros(8,1);
snr_cnc4 = zeros(8,1);
snr_convex = zeros(8,1);
lambdas = 80:-10:10;
tic
[xhat1, vhat1, res_norm_hist1] = srls_GMC_matrix(y, 'deblurring', lambdas(1), H=H, gamma=0.8, acceleration = 'original', splitting='FB');
snr_cnc1(1) = snr(x,xhat1-x);
for i=2:8
    [xhat1, vhat1, res_norm_hist1] = srls_GMC_matrix(y, 'deblurring', lambdas(i), H=H, gamma=0.8, acceleration = 'original', splitting='FB', xv0 = [xhat1(:);vhat1(:)]);
    snr_cnc1(i) = snr(x,xhat1-x);
end
time1 = toc
tic
[xhat2, vhat2, res_norm_hist2] = srls_GMC_matrix(y, 'deblurring', lambdas(1), H=H, gamma=0.8, acceleration = 'aa2', splitting='FB');
snr_cnc2(1) = snr(x,xhat2-x);
for i=2:8
    [xhat2, vhat2, res_norm_hist2] = srls_GMC_matrix(y, 'deblurring', lambdas(i), H=H, gamma=0.8, acceleration = 'aa2', splitting='FB', xv0 = [xhat2(:);vhat2(:)]);
    snr_cnc2(i) = snr(x,xhat2-x);
end
time2 = toc
tic
[xhat3, vhat3, res_norm_hist3] = srls_GMC_matrix(y, 'deblurring', lambdas(1), H=H, gamma=0.8, acceleration = 'original', splitting='FBF');
snr_cnc3(1) = snr(x,xhat3-x);
for i=2:8
    [xhat3, vhat3, res_norm_hist3] = srls_GMC_matrix(y, 'deblurring', lambdas(i), H=H, gamma=0.8, acceleration = 'original', splitting='FBF', xv0 = [xhat3(:);vhat3(:)]);
    snr_cnc3(i) = snr(x,xhat3-x);
end
time3 = toc
tic
[xhat4, vhat4, res_norm_hist4] = srls_GMC_matrix(y, 'deblurring', lambdas(1), H=H, gamma=0.8, acceleration = 'aa2', splitting='FBF');
snr_cnc4(1) = snr(x,xhat4-x);
for i=2:8
    [xhat4, vhat4, res_norm_hist4] = srls_GMC_matrix(y, 'deblurring', lambdas(i), H=H, gamma=0.8, acceleration = 'aa2', splitting='FBF', xv0 = [xhat4(:);vhat4(:)]);
    snr_cnc4(i) = snr(x,xhat4-x);
end
time4 = toc
tic
[xhat5, vhat5, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(1), H=H, gamma=0, acceleration = 'aa2', splitting='FB');
snr_convex(1) = snr(x,xhat5-x);
for i=2:8
    [xhat5, vhat5, ~] = srls_GMC_matrix(y, 'deblurring', lambdas(i), H=H, gamma=0, acceleration = 'aa2', splitting='FB', xv0 = [xhat5(:);vhat5(:)]);
    snr_convex(i) = snr(x,xhat5-x);
end
time5 = toc

[xhat, ~, ~ ] = srls_GMC_matrix(y, 'deblurring', 60, H=H, gamma=0.8, acceleration = 'aa2', splitting='FB'); 

t=tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 960 640])
nexttile
imshow(uint8(x))
title('original')
nexttile
imshow(uint8(y))
title('blurred and noisy')
nexttile
imshow(uint8(xhat))
title('restored by CNC')
nexttile
plot(log(res_norm_hist1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
title('FB, lambda=10')
l = legend('show','Location','northeast')
nexttile
plot(log(res_norm_hist3), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist4), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
title('FBF, lambda=10')
l = legend('show','Location','northeast')
nexttile
plot(lambdas,snr_convex,'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Convex, FB, aa2 (%0.0f s)',time5))
hold on
plot(lambdas,snr_cnc1,'c-', 'LineWidth', 1.5, 'DisplayName', sprintf('FB, original (%0.0f s)',time1))
hold on
plot(lambdas,snr_cnc2,'m-', 'LineWidth', 1.5, 'DisplayName', sprintf('FB, aa2 (%0.0f s)',time2))
hold on
plot(lambdas,snr_cnc3,'g-', 'LineWidth', 1.5, 'DisplayName', sprintf('FBF, original (%0.0f s)',time3))
hold on
plot(lambdas,snr_cnc4,'-','color','#77AC30','LineWidth', 1.5, 'DisplayName', sprintf('FBF, aa2 (%0.0f s)',time4))
l = legend('show','Location','southeast')
xlabel('lambda');
ylabel('SNR');
xlim([10,80])
title("SNR vs lambda")