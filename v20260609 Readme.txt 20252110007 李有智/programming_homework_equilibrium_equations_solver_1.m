%% ========================================================================
% 2.4 作业完整代码（无2.3桁架部分）
% 包含：算例1（三对角矩阵）、算例2（非正定检测）、任务2（病态矩阵误差分析）、算例4（Poisson有限元）
% 输出：results_2_4.txt + 云图
% ========================================================================

clear; clc; close all;

% 打开输出文件
fid = fopen('results_2_4.txt', 'w');
fprintf(fid, '2.4 有限元平衡方程组求解与误差分析结果\n');
fprintf(fid, '==========================================\n\n');

%% ========================= 算例1 & 算例2 & 病态矩阵分析 =========================
run_ldlt_tests(fid);

%% ========================= 算例4：二维 Poisson 方程有限元求解 =========================
fprintf(fid, '\n========== 算例4：二维 Poisson 方程有限元求解 (线性三角形元) ==========\n');
run_poisson_fem(fid);

fprintf(fid, '\n程序运行完毕。\n');
fclose(fid);
fprintf(1, '所有结果已写入 results_2_4.txt\n');

%% ======================== 以下为函数定义 ================================

%% ---------- 稠密 LDL^T 分解与求解函数 ----------
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

function [r_norm, rel_r_norm] = residual_norm(K, a, R)
    r = R - K * a;
    r_norm = norm(r);
    rel_r_norm = r_norm / norm(R);
end

%% ---------- 算例1 & 算例2 & 病态矩阵分析的主函数 ----------
function run_ldlt_tests(fid)
    % 算例1：三对角对称正定矩阵
    fprintf(fid, '========== 算例1：三对角对称正定矩阵 ==========\n');
    n_list = [10, 100, 500];
    fprintf(fid, '   n       求解时间(s)    相对残差      相对误差\n');
    for n = n_list
        K = diag(2*ones(n,1)) + diag(-1*ones(n-1,1),1) + diag(-1*ones(n-1,1),-1);
        a_exact = ones(n,1);
        R = K * a_exact;
        
        tic;
        [L, D, status] = ldlt_factor(K);
        if status ~= 0
            fprintf(fid, '%4d     LDL^T 分解失败（非正定）\n', n);
            continue;
        end
        a_num = ldlt_solve(L, D, R);
        t = toc;
        [~, rel_res] = residual_norm(K, a_num, R);
        rel_err = norm(a_num - a_exact) / norm(a_exact);
        fprintf(fid, '%4d     %12.4e      %8.2e      %8.2e\n', n, t, rel_res, rel_err);
    end
    
    % 算例2：非正定矩阵检测
    fprintf(fid, '\n========== 算例2：非正定矩阵检测 ==========\n');
    K_bad = [1, 2; 2, 1];
    [~, ~, status] = ldlt_factor(K_bad);
    if status == 1
        fprintf(fid, '✓ 成功检测到非正定主元，分解中止。\n');
    else
        fprintf(fid, '✗ 未检测到非正定，算法有误。\n');
    end
    fprintf(fid, '说明：缺少边界约束时，总体刚度矩阵可能出现零主元或负主元。\n');
    
    % 任务2：病态矩阵误差分析
    fprintf(fid, '\n========== 病态矩阵误差分析 ==========\n');
    K_ill = [1.0, 1.0; 1.0, 1.0001];
    a_exact = [1; 1];
    R_ill = K_ill * a_exact;
    condK = cond(K_ill);
    fprintf(fid, '条件数 cond(K) = %.2f\n', condK);
    
    % 双精度求解
    [L, D, ~] = ldlt_factor(K_ill);
    a_double = ldlt_solve(L, D, R_ill);
    [~, rel_res_d] = residual_norm(K_ill, a_double, R_ill);
    err_d = norm(a_double - a_exact) / norm(a_exact);
    fprintf(fid, '双精度: 解=[%.8f, %.8f], 相对残差=%.2e, 相对误差=%.2e\n', ...
            a_double(1), a_double(2), rel_res_d, err_d);
    
    % 模拟4位有效数字舍入（四舍五入到4位小数）
    K4 = round(K_ill, 4);
    R4 = round(R_ill, 4);
    [L4, D4, ~] = ldlt_factor(K4);
    a_4 = ldlt_solve(L4, D4, R4);
    [~, rel_res_4] = residual_norm(K4, a_4, R4);
    err_4 = norm(a_4 - a_exact) / norm(a_exact);
    fprintf(fid, '4位精度: 解=[%.8f, %.8f], 相对残差=%.2e, 相对误差=%.2e\n', ...
            a_4(1), a_4(2), rel_res_4, err_4);
    fprintf(fid, '结论：病态问题中残差小不一定解准确，条件数放大了误差。\n');
end

%% ---------- 算例4：二维 Poisson 方程有限元求解 ----------
function run_poisson_fem(fid)
    nx_list = [50, 100, 200];   % 网格规模
    for idx = 1:length(nx_list)
        nx = nx_list(idx);
        ny = nx;
        [coords, elements] = generate_tri_mesh(0, 1, 0, 1, nx, ny);
        nNodes = size(coords, 1);
        nElem = size(elements, 1);
        
        % 装配
        t_assemble = tic;
        K = sparse(nNodes, nNodes);
        R = zeros(nNodes, 1);
        for e = 1:nElem
            nodes = elements(e, :);
            x = coords(nodes, 1);
            y = coords(nodes, 2);
            [Ke, Re] = tri3_poisson_ke_rhs(x, y, @poisson_f);
            K(nodes, nodes) = K(nodes, nodes) + Ke;
            R(nodes) = R(nodes) + Re;
        end
        t_assemble = toc(t_assemble);
        
        % Dirichlet 边界条件 u=0
        fixed_dof = find_boundary_nodes(coords);
        free_dof = setdiff(1:nNodes, fixed_dof);
        K_ff = K(free_dof, free_dof);
        R_f = R(free_dof);
        
        % 求解
        t_solve = tic;
        u_free = K_ff \ R_f;
        t_solve = toc(t_solve);
        
        % 重构位移
        u = zeros(nNodes, 1);
        u(free_dof) = u_free;
        
        % 精确解
        u_exact = sin(pi * coords(:,1)) .* sin(pi * coords(:,2));
        err_max = max(abs(u - u_exact));
        err_L2 = sqrt( sum((u - u_exact).^2) / sum(u_exact.^2) );
        rel_res = norm(R - K*u) / norm(R);
        
        % 输出
        fprintf(fid, '\n--- 网格 %dx%d ---\n', nx, ny);
        fprintf(fid, '节点数: %d, 单元数: %d\n', nNodes, nElem);
        fprintf(fid, '非零元个数: %d\n', nnz(K_ff));
        fprintf(fid, '装配时间: %.3f s, 求解时间: %.3f s\n', t_assemble, t_solve);
        fprintf(fid, '相对残差: %.2e\n', rel_res);
        fprintf(fid, '最大节点误差: %.4e, 相对 L2 误差: %.4e\n', err_max, err_L2);
        
        % 绘制云图（仅最小网格）
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
    fprintf(fid, '\n求解器: MATLAB 内置反斜杠 (\\), 自动调用 Intel MKL PARDISO / SuiteSparse\n');
    fprintf(fid, '存储格式: CSC (Compressed Sparse Column)\n');
end

%% ---------- 网格生成辅助函数 ----------
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

%% ---------- 三角形单元刚度矩阵和载荷向量（修正版）----------
function [Ke, Re] = tri3_poisson_ke_rhs(x, y, f_func)
    % x, y 为 3x1 列向量
    % 计算面积：使用 [1,1,1; x'; y'] 矩阵的行列式
    area = abs(det([1,1,1; x'; y'])) / 2;   % 修正：分号分隔行
    B = [y(2)-y(3), y(3)-y(1), y(1)-y(2);
         x(3)-x(2), x(1)-x(3), x(2)-x(1)] / (2*area);
    Ke = area * (B' * B);
    centroid_x = mean(x);
    centroid_y = mean(y);
    f_center = f_func(centroid_x, centroid_y);
    Re = (area / 3) * f_center * ones(3,1);
end

%% ---------- 右端项 f(x,y) ----------
function fval = poisson_f(x, y)
    fval = 2 * pi^2 * sin(pi*x) .* sin(pi*y);
end

%% ---------- 边界节点识别 ----------
function fixed_dof = find_boundary_nodes(coords)
    tol = 1e-10;
    x = coords(:,1);
    y = coords(:,2);
    on_boundary = (abs(x-0) < tol) | (abs(x-1) < tol) | (abs(y-0) < tol) | (abs(y-1) < tol);
    fixed_dof = find(on_boundary);
end

%% ---------- 稀疏矩阵 CSR 转换示例（仅用于展示，未被调用）----------
function [rowptr, colind, values] = sparse_to_csr(A)
    [i, j, val] = find(A);
    n = size(A,1);
    rowptr = zeros(n+1, 1);
    for k = 1:length(i)
        rowptr(i(k)+1) = rowptr(i(k)+1) + 1;
    end
    rowptr = cumsum(rowptr);
    rowptr = [0; rowptr(1:end-1)];
    colind = j - 1;
    values = val;
end