%% ========================================================================
% 2.4 作业完整代码：满足输入输出要求（已修正 vertcat 错误）
% 包含算例0（2.3桁架）、算例1（三对角）、算例2（非正定）、算例4（Poisson）
% 输出：results_2_4.txt + 云图
% ========================================================================

clear; clc; close all;

% 输出文件
fid = fopen('results_2_4.txt', 'w');
fprintf(fid, '2.4 作业：平衡方程组求解与误差分析\n');
fprintf(fid, '====================================\n\n');

%% ========================= 算例0：2.3 桁架结构（复用2.3模块）=========================
fprintf(fid, '========== 算例0：2.3 一维两单元杆结构 ==========\n');
K_ff0 = [300, -200; -200, 200];
rhs0 = [0; 10];
[L0, D0, status0] = ldlt_factor(K_ff0);
if status0 == 0
    d_F0 = ldlt_solve(L0, D0, rhs0);
    fprintf(fid, 'LDL^T 分解成功，节点2位移 u2 = %.6f (理论0.1)，节点3位移 u3 = %.6f (理论0.15)\n', d_F0(1), d_F0(2));
    fprintf(fid, '解向量 = [%.6f; %.6f]\n', d_F0(1), d_F0(2));
else
    fprintf(fid, 'LDL^T 分解失败\n');
end

fprintf(fid, '\n========== 算例0：2.3 二维两杆桁架结构 ==========\n');
K_ff0_2 = [0.35355339, 0.35355339; 0.35355339, 1.35355339];
rhs0_2 = [10; 0];
[L0_2, D0_2, status0_2] = ldlt_factor(K_ff0_2);
if status0_2 == 0
    d_F0_2 = ldlt_solve(L0_2, D0_2, rhs0_2);
    fprintf(fid, '节点3位移 u3 = %.6f (理论38.284271)，v3 = %.6f (理论-10.000000)\n', d_F0_2(1), d_F0_2(2));
    fprintf(fid, '解向量 = [%.6f; %.6f]\n', d_F0_2(1), d_F0_2(2));
end

%% ========================= 算例1：三对角对称正定矩阵 =========================
fprintf(fid, '\n========== 算例1：三对角对称正定矩阵 ==========\n');
n_list = [10, 100, 500, 1000];
fprintf(fid, '%-6s %-12s %-12s %-12s %-12s %-20s\n', 'n', '求解时间(s)', '相对残差', '相对误差', '非零元数', '解向量(前3个)');
for n = n_list
    e = ones(n,1);
    K = spdiags([-e, 2*e, -e], -1:1, n, n);
    a_exact = ones(n,1);
    R = K * a_exact;
    
    tic;
    a_num = K \ R;
    t = toc;
    
    rel_res = norm(R - K*a_num) / norm(R);
    rel_err = norm(a_num - a_exact) / norm(a_exact);
    nz = nnz(K);
    if n >= 3
        sol_str = sprintf('[%.4f, %.4f, %.4f,...]', a_num(1), a_num(2), a_num(3));
    else
        sol_str = sprintf('%s', mat2str(a_num'));
    end
    fprintf(fid, '%-6d %-12.4e %-12.2e %-12.2e %-12d %-20s\n', n, t, rel_res, rel_err, nz, sol_str);
end
fprintf(fid, '求解器: MATLAB 内置稀疏求解器 (Intel MKL PARDISO / SuiteSparse)\n');
fprintf(fid, '存储格式: CSC (Compressed Sparse Column)\n\n');

%% ========================= 算例2：非正定矩阵检测 =========================
fprintf(fid, '========== 算例2：非正定矩阵检测 ==========\n');
K_bad = [1, 2; 2, 1];
R_bad = [1; 1];
fprintf(fid, '测试矩阵 K = [[1,2];[2,1]], 右端项 R = [1;1]\n');
[~, D_bad, status_bad] = ldlt_factor(K_bad);
if status_bad == 1
    fprintf(fid, '✓ LDL^T 分解检测到非正主元 (D(2)=%.2e <=0)，正确报错。\n', D_bad(2));
else
    fprintf(fid, '未检测到非正主元，但矩阵非正定，结果可能错误。\n');
end
fprintf(fid, '说明：缺少边界约束时，总体刚度矩阵可能出现零主元或负主元。\n\n');

%% ========================= 算例4：二维 Poisson 方程有限元求解 =========================
fprintf(fid, '========== 算例4：二维 Poisson 方程 (线性三角形元) ==========\n');
fprintf(fid, '求解器: MATLAB 内置稀疏求解器 (CSC 格式)\n');
nx_list = [50, 100, 200];
for idx = 1:length(nx_list)
    nx = nx_list(idx);
    ny = nx;
    [coords, elements] = generate_tri_mesh(0, 1, 0, 1, nx, ny);
    nNodes = size(coords,1);
    nElem = size(elements,1);
    
    t_assemble = tic;
    K = sparse(nNodes, nNodes);
    R = zeros(nNodes,1);
    for e = 1:nElem
        nodes = elements(e,:);
        x = coords(nodes,1);
        y = coords(nodes,2);
        [Ke, Re] = tri3_poisson_ke_rhs(x, y, @poisson_f);
        K(nodes, nodes) = K(nodes, nodes) + Ke;
        R(nodes) = R(nodes) + Re;
    end
    t_assemble = toc(t_assemble);
    
    fixed_dof = find_boundary_nodes(coords);
    free_dof = setdiff(1:nNodes, fixed_dof);
    K_ff = K(free_dof, free_dof);
    R_f = R(free_dof);
    
    t_solve = tic;
    u_free = K_ff \ R_f;
    t_solve = toc(t_solve);
    
    u = zeros(nNodes,1);
    u(free_dof) = u_free;
    u_exact = sin(pi * coords(:,1)) .* sin(pi * coords(:,2));
    
    err_max = max(abs(u - u_exact));
    err_L2 = sqrt(sum((u - u_exact).^2) / sum(u_exact.^2));
    rel_res = norm(R - K*u) / norm(R);
    
    fprintf(fid, '\n--- 网格 %dx%d ---\n', nx, ny);
    fprintf(fid, '节点数: %d, 单元数: %d\n', nNodes, nElem);
    fprintf(fid, '非零元个数: %d\n', nnz(K_ff));
    fprintf(fid, '装配时间: %.3f s, 求解时间: %.3f s\n', t_assemble, t_solve);
    fprintf(fid, '相对残差: %.2e\n', rel_res);
    fprintf(fid, '最大节点误差: %.4e, 相对 L2 误差: %.4e\n', err_max, err_L2);
    
    fprintf(fid, '解向量 (前5个节点):\n');
    for i = 1:min(5, nNodes)
        fprintf(fid, '  u(%d)=%.6f (精确:%.6f)\n', i, u(i), u_exact(i));
    end
    
    if nx == 50
        figure;
        trisurf(elements, coords(:,1), coords(:,2), u);
        title(sprintf('Poisson 数值解 (nx=ny=%d)', nx));
        xlabel('x'); ylabel('y'); zlabel('u');
        shading interp; colormap jet; colorbar;
        saveas(gcf, sprintf('poisson_nx%d.png', nx));
        
        figure;
        trisurf(elements, coords(:,1), coords(:,2), u - u_exact);
        title(sprintf('误差分布 (nx=ny=%d)', nx));
        xlabel('x'); ylabel('y'); zlabel('Error');
        shading interp; colormap jet; colorbar;
        saveas(gcf, sprintf('error_nx%d.png', nx));
    end
end

fprintf(fid, '\n程序运行完毕。\n');
fclose(fid);
fprintf(1, '结果已写入 results_2_4.txt，并生成云图。\n');

%% ======================= 辅助函数 ==========================
function [L, D, status] = ldlt_factor(K)
    n = size(K,1);
    L = eye(n);
    D = zeros(n,1);
    A = K;
    for j = 1:n
        d = A(j,j);
        for p = 1:j-1
            d = d - L(j,p)^2 * D(p);
        end
        if d <= 1e-12
            status = 1;
            return;
        end
        D(j) = d;
        for i = j+1:n
            l = A(i,j);
            for p = 1:j-1
                l = l - L(i,p) * D(p) * L(j,p);
            end
            L(i,j) = l / D(j);
        end
    end
    status = 0;
end

function x = ldlt_solve(L, D, b)
    n = length(b);
    y = zeros(n,1);
    for i = 1:n
        y(i) = b(i);
        for j = 1:i-1
            y(i) = y(i) - L(i,j) * y(j);
        end
    end
    z = y ./ D;
    x = zeros(n,1);
    for i = n:-1:1
        x(i) = z(i);
        for j = i+1:n
            x(i) = x(i) - L(j,i) * x(j);
        end
    end
end

function [coords, elements] = generate_tri_mesh(xmin, xmax, ymin, ymax, nx, ny)
    x = linspace(xmin, xmax, nx+1);
    y = linspace(ymin, ymax, ny+1);
    [X, Y] = meshgrid(x, y);
    coords = [X(:), Y(:)];
    nNodeX = nx+1;
    nNodeY = ny+1;
    elements = [];
    for j = 1:ny
        for i = 1:nx
            n1 = (j-1)*nNodeX + i;
            n2 = (j-1)*nNodeX + i+1;
            n3 = j*nNodeX + i;
            n4 = j*nNodeX + i+1;
            elements = [elements; n1, n2, n3; n2, n4, n3];
        end
    end
end

function [Ke, Re] = tri3_poisson_ke_rhs(x, y, f_func)
    % 修正 vertcat 错误：显式构建 3×3 矩阵
    % x, y 是长度为 3 的列向量
    A_mat = [1, 1, 1;
             x(1), x(2), x(3);
             y(1), y(2), y(3)];
    area = abs(det(A_mat)) / 2;
    
    B = [y(2)-y(3), y(3)-y(1), y(1)-y(2);
         x(3)-x(2), x(1)-x(3), x(2)-x(1)] / (2*area);
    Ke = area * (B' * B);
    centroid_x = mean(x);
    centroid_y = mean(y);
    f_center = f_func(centroid_x, centroid_y);
    Re = (area / 3) * f_center * ones(3,1);
end

function fval = poisson_f(x, y)
    fval = 2 * pi^2 * sin(pi*x) .* sin(pi*y);
end

function fixed_dof = find_boundary_nodes(coords)
    tol = 1e-10;
    x = coords(:,1);
    y = coords(:,2);
    on_boundary = (abs(x-0) < tol) | (abs(x-1) < tol) | (abs(y-0) < tol) | (abs(y-1) < tol);
    fixed_dof = find(on_boundary);
end