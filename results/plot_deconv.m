load deconv3;
t=tiledlayout(2,3, 'Padding', 'none', 'TileSpacing', 'compact'); 
set(gcf,'Position',[10 10 960 640])
nexttile
imshow(uint8(x))
title('original')
nexttile
imshow(uint8(y))
title('blurred and noisy')
nexttile
imshow(uint8(xhat1))
title('restored by convex, SNR=24.09')
nexttile
imshow(uint8(xhat2))
title('restored by CNC, SNR=25.13')
nexttile
plot(log(res_norm_hist1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
title('FBS, $\lambda=10$, $\gamma=0.8$','Interpreter','latex')
l = legend('show','Location','northeast')
nexttile
plot(log(res_norm_hist3), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Original')
hold on
plot(log(res_norm_hist4), 'r-', 'LineWidth', 1.5, 'DisplayName', 'AA')
xlabel('iteration');
ylabel('residual norm (log scale)');
title('FBFS, $\lambda=10$, $\gamma=0.8$','Interpreter','latex')
l = legend('show','Location','northeast')
l = legend('show','Location','southeast')