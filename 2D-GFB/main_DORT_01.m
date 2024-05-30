% clc
clear, close all
eps0 = 8.8541878128e-12; mu0  = 4.0*pi*1.0e-7; c0 = 1.0/sqrt(eps0*mu0);
freq = c0; landa = c0/freq;

N = 21; L = 2;

r = linspace(-5.0, 5.0, 21)*landa; r = [r; 0.0*ones(size(r))];
% r = linspace(-5.0, 5.0, N)*landa; r = [r; 0.0*ones(size(r))];
% r = (-N/2+1:N/2)*landa/2; r = [r; 0.5*ones(size(r))];

% x = [+0.5*landa, -0.25*landa; -7*landa, -8*landa]; L = 2;
x = [-1.0*landa, 2.0*landa; -6*landa, -10*landa]; % L = 1;

% x = [+2.0*landa, 0.0*landa, -1.0*landa; -8*landa, -6*landa, -5*landa]; L = 3;
% x = [-1.0*landa; -6*landa]; L = 1;
% x = [+0.0*landa; -5*landa]; L = 1;
figure; plot(r(1, :), r(2, :), 'o'); hold on; plot(x(1, :), x(2, :), 'xr')
G = zeros(N, L);
for j0 = 1:L
    gamma = zeros(N, 1);
    for i0 = 1:N
        gamma(i0) = G_func(r(:, i0), x(:, j0), landa);
%         gamma(i0) = awgn(gamma(i0), 10, 'measured');
    end
    G(:, j0) = gamma;
end

% R = [1.0, 0; 0, 1.0];
R = eye(L, L);
% R = [1.0];

H = G*R*(G.');
[U, S, V] = svd(H); 
figure; stem(1:N, diag(S)/max(diag(S)), 'k', 'LineWidth', 3); xlabel('number', 'FontSize',20)
ylabel('Singular Value', 'FontSize',20); grid on; axis([1 N 0 1]); ax = gca; % current axes
ax.FontSize = 20;
% plot(x,y,'--gs',...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5])

% [V, D] = eig((H')*H);
% [d, ind] = sort(diag(D), 'descend'); D = D(ind,ind); V = V(:,ind);

si = 1;
V2 = U;
% v1 = sqrt(S(si, si))*V2(:, si);
v1 = S(si, si)*S(si, si)*V2(:, si);

% si = 1:2;
% V2 = U;
% v1 = V2(:, si);
% for i0 = 1:length(si)
%     v1(:, i0) = v1(:, i0)*sqrt(S(i0, i0));
% end

% si = 2;
% V2 = U;
% v1 = sqrt(S(si, si))*V2(:, si);
% bachpropagation
step = landa/10.0;
xo1 = (-5*landa:step:5*landa);          nxo = length(xo1);
yo1 = (-12*landa:step:0*landa);         nyo = length(yo1);
[xo2, yo2] = meshgrid(xo1, yo1); xo = xo2(:); yo = yo2(:);
% plot(yo, xo, 'x'); axis('xy')

Dort_res = zeros(length(xo), 1);
for i0 = 1:length(xo)
    gamma = zeros(N, 1);
    x_ = [xo(i0); yo(i0)];
    for j0 = 1:N
        gamma(j0) = G_func(r(:, j0), x_, landa);
    end
%     Dort_res(i0) = (gamma.')*v1;
    Dort_res(i0) = sum((gamma')*v1)^2;
%     Dort_res(i0) = sum((v1')*(1./gamma));
end
norm_m = reshape(Dort_res, [nyo, nxo]);
figure; imagesc(xo1, yo1, abs(norm_m)/max(abs(norm_m(:)))); axis xy; 
hold on; 
% plot(-0.5, -0.5, 'or');  
plot(r(1, :), r(2, :), 'mo', 'LineWidth', 3, 'MarkerSize', 15); hold on; plot(x(1, :), x(2, :), 'xr', 'LineWidth', 3, 'MarkerSize', 15)
xlabel('x-axis (m)', 'FontSize',24)
ylabel('y-axis (m)', 'FontSize',24)
ax = gca; % current axes
ax.FontSize = 30; grid on
legend('Antenna', 'Target')

res_ = abs(norm_m)/max(abs(norm_m(:)));
[M, I] = max(res_);
[~, I2] = max(M); plot(xo1((I2)), yo1(I(I2)), 'sk', 'LineWidth', 3, 'MarkerSize', 15)
legend('Antenna', 'Estimated Object', 'Actual Object')

% err = sqrt((xo1(I2)-x(1,1))^2 +  (yo1(I(I2))-x(2,1))^2)
err = sqrt((xo1(I2)-x(1,2))^2 +  (yo1(I(I2))-x(2,2))^2)

% zlabel('P_M_U_S_I_C', 'FontSize',20)
% figure; mesh(xo2, yo2, norm_m); axis xy
% xlabel('x-axis (m)', 'FontSize',20)
% ylabel('y-axis (m)', 'FontSize',20)
% zlabel('P_M_U_S_I_C', 'FontSize',20)

