x = double(checkerboard(32));
y = awgn(x,15,'measured');
snr_cnc1 = zeros(10,1);
snr_cnc2 = zeros(10,1);
snr_cnc3 = zeros(10,1);
snr_cnc4 = zeros(10,1);
snr_convex = zeros(10,1);
lambdas = 1:1:10;
mask = binornd(1, 0.2, size(x,1), size(x,2));
i = 1;
[xhat1, vhat1, res_norm_hist1] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(i), mask=mask, gamma=0.8, acceleration = 'original', splitting='DY',printevery=1);
[xhat2, vhat1, res_norm_hist1] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(i), mask=mask, gamma=0.8, acceleration = 'original', splitting='DY',xv0 = [y(:)+randn(256*256,1);y(:)+randn(256*256,1)],printevery=1);

tic
for i=1:10
    [xhat1, vhat1, res_norm_hist1] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(i), mask=mask, gamma=0.8, acceleration = 'original', splitting='FB');
    snr_cnc1(i) = snr(x,xhat1-x);
end
time1 = toc
tic
for i=1:10
    [xhat2, vhat2, res_norm_hist2] = srls_GMC_matrix(y.*mask, 'matrix completion', 1, mask=mask, gamma=0.8, acceleration = 'original', splitting='FB');
    snr_cnc2(i) = snr(x,xhat2-x);
end
time2 = toc
tic
for i=1:10
    [xhat3, vhat3, res_norm_hist3] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(i), mask=mask, gamma=0.8, acceleration = 'original', splitting='FBF');
    snr_cnc3(i) = snr(x,xhat3-x);
end
time3 = toc
tic
for i=1:10
    [xhat4, vhat4, res_norm_hist4] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(i), mask=mask, gamma=0.8, acceleration = 'original', splitting='DY', lower = 0, upper = 1);
    snr_cnc4(i) = snr(x,xhat4-x);
end
time4 = toc

tic
for i=1:10
    [xhat5, vhat5, res_norm_hist5] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(i), mask=mask, gamma=0, acceleration = 'aa2', splitting='FB');
    snr_convex(i) = snr(x,xhat5-x);
end
time5 = toc


tic
[xhat1, ~, res_norm_hist1] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(1), mask=mask, gamma=0.8, acceleration = 'original', splitting='FB');
time_FB = toc
tic
[xhat2, ~, res_norm_hist2] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(1), mask=mask, gamma=0.8, acceleration = 'aa2', splitting='FB');
time_FB_aa2 = toc
tic
[xhat3, ~, res_norm_hist3] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(1), mask=mask, gamma=0.8, acceleration = 'original', splitting='FBF');
time_FBF = toc
tic
[xhat4, ~, res_norm_hist4] = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(1), mask=mask, gamma=0.8, acceleration = 'aa2', splitting='FBF');
time_FBF_aa2 = toc

xhat4 = srls_GMC_matrix(y.*mask, 'matrix completion', lambdas(5), mask=mask, gamma=0.8, acceleration = 'aa2', splitting='FBF');



t=tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 960 640])
nexttile
imshow(uint8(x))
title('original')
nexttile
imshow(uint8(y.*mask + 255.*(1-mask)))
title('noisy, 80% pixels missing')
nexttile
imshow(uint8(xhat4))
title('restored')
nexttile
plot(log(res_norm_hist1), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('original (%0.0f s)',time_FB))
hold on
plot(log(res_norm_hist2), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('aa2 (%0.0f s)',time_FB_aa2) )
xlabel('iteration');
ylabel('norm of redidual (log scale)');
title('FB, lambda=100')
l = legend('show','Location','northeast')
nexttile
plot(log(res_norm_hist3), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('original (%0.0f s)',time_FBF))
hold on
plot(log(res_norm_hist4), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('aa2 (%0.0f s)',time_FBF_aa2) )
xlabel('iteration');
ylabel('norm of redidual (log scale)');
title('FBF, lambda=100')
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
l = legend('show','Location','southwest')
xlabel('lambda');
ylabel('SNR');
xlim([50 1050])
title("SNR vs lambda")

times_convex = csvread('times_convex.csv');
times_convex_original = csvread('times_convex_original.csv');