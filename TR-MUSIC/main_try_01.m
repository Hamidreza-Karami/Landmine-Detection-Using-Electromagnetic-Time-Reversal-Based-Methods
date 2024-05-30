
clc, clear, close all

file_temp1 = 'w_obj_r';
file_temp2 = 'w_obj_i';

% import frequency
sr = import_real_part('w_obj_r1_1.txt');
freq = sr(:,1);

na = 3; nf = length(freq);
s_w_obj = zeros(na, na, nf);
for i0 = 1:3
    for j0 = 1:3
        
        file_name = [file_temp1, num2str(i0), '_', num2str(j0), '.txt'];
        sr = import_real_part(file_name);
%         plot(sr(:,1), sr(:,2)); hold on
        file_name = [file_temp2, num2str(i0), '_', num2str(j0), '.txt'];
        si = import_imag_part(file_name);
%         plot(si(:,1), si(:,2)); hold off
        s_w_obj(i0, j0, 1:nf) = sr(:,2) + 1j*si(:,2);
        
    end
end

file_temp1 = 'wo_obj_r';
file_temp2 = 'wo_obj_i';

sr = import_real_part('wo_obj_r1_1.txt');
freq_wo = sr(:,1);
% nf2 = length(freq_wo);
s_wo_obj = zeros(na, na, nf);
for i0 = 1:3
    for j0 = 1:3
        
        file_name = [file_temp1, num2str(i0), '_', num2str(j0), '.txt'];
        sr = import_real_part(file_name);
        sr2 = interp1(freq_wo, sr, freq);
%         figure; plot(sr(:,1), sr(:,2)); hold on; plot(sr2(:,1), sr2(:,2), '--r');
        file_name = [file_temp2, num2str(i0), '_', num2str(j0), '.txt'];
        si = import_imag_part(file_name);
        si2 = interp1(freq_wo, si, freq);
%         figure; plot(si(:,1), si(:,2)); hold on; plot(si2(:,1), si2(:,2), '--r');
        s_wo_obj(i0, j0, 1:nf) = sr2(:,2) + 1j*si2(:,2);
        
    end
end

f_ind = 501; freq(f_ind)
H = s_w_obj(1:na, 1:na, f_ind) - s_wo_obj(1:na, 1:na, f_ind);
[U, S, V] = svd(H); figure; stem(1:3, diag(S)/max(diag(S)), 'k', 'LineWidth', 3)
xlabel('number','FontSize',24)
ylabel('Singular value','FontSize',24)
%%
L = 1;
uo = U(:, L+1:end);

e_r = import_E('E1_y_r.txt'); e1y = e_r(:, 6) + 1j*e_r(:, 7);
e_r = import_E('E2_y_r.txt'); e2y = e_r(:, 6) + 1j*e_r(:, 7);
e_r = import_E('E3_y_r.txt'); e3y = e_r(:, 6) + 1j*e_r(:, 7);
xs = e_r(:, 1); zs = e_r(:, 3); 

ni0 = size(e_r, 1);
norm_music = zeros(ni0, 1);
for i0 = 1:ni0
    gamma = zeros(3, 1);
    gamma(1) = e1y(i0, 1);
    gamma(2) = e2y(i0, 1);
    gamma(3) = e3y(i0, 1);
%     gamma(4) = x4(i0, 3);
%     gamma(5) = x5(i0, 3);
    norm_music(i0) = 1.0/(norm((gamma')*uo)^2 + 0*eps + 0.0);
end

xo = linspace(-200, 200, 240);
yo = linspace(50, 500, 240);
[xq, yq] = meshgrid(xo, yo);
vq = griddata(xs, zs, norm_music(:), xq, yq);

figure; imagesc(xo, yo, abs(vq)/max(vq(:))); hold on; plot(-0, 300, 'or', 'MarkerSize', 8, 'LineWidth', 3); 
% plot(0, 150, 'or');
axis equal; axis([-175 175 50 500]); colorbar
ax = gca;
ax.FontSize = 24;
xlabel('x-axis (m)','FontSize',24)
ylabel('z-axis (m)','FontSize',24)
legend('Actual Location')

norm_m = abs(vq)/max(vq(:));
[~, ind] = max(norm_m(:));
plot(xq(ind), yq(ind), 'xk', 'MarkerSize', 12, 'LineWidth', 3)
err = sqrt((xq(ind)+0)^2 + (yq(ind)-300)^2)*1

%%







