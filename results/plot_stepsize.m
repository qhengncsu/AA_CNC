gammas = 0:0.01:0.99;
mus_fb = zeros(100,1);
mus_fbf = zeros(100,1);
for i=1:100
    gamma = gammas(i);
    gamma_matrix = [1-gamma,gamma;-gamma,gamma];
    mus_fb(i) = 2*min(1,(1-gamma)/gamma);
    mus_fbf(i) = 1/norm(gamma_matrix);
end
plot(gammas,mus_fb, 'b-', 'LineWidth', 1.5, 'DisplayName', 'FBS')
hold on
plot(gammas,mus_fbf, 'r-', 'LineWidth', 1.5, 'DisplayName', 'FBFS')
l = legend('show','Location','northeast')
xlabel('$$\gamma$$','Interpreter','latex');
ylabel('step size');
xlim([0,0.99])
ylim([0,2.2])
title("Maximum step size: FBS vs FBFS")