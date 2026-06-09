%% ========================================================================
% 总体刚度矩阵组装与桁架结构求解程序（2D/3D）
% 2.3作业模块 + 使用MATLAB内置求解器（\）与奇异矩阵处理
% 包含：两个验证算例 + 附加题（三维桁架）
% 输出：屏幕显示 + results.txt 文件
% 求解器：非奇异时用反斜杠，奇异时自动切换为伪逆
% ========================================================================

clear; clc; close all;

% 打开输出文件
fid = fopen('results.txt', 'w');
fprintf(fid, '总体刚度矩阵组装与桁架结构求解结果 (MATLAB内置求解器+奇异处理)\n');
fprintf(fid, '================================================================\n\n');

%% --------------------------- 1. 前处理模块 ---------------------------
function model = preprocess_example1()
    model.title = '一维两单元杆结构';
    model.nsd = 1;
    model.ndof = 2;
    model.nnp = 3;
    model.nel = 2;
    model.nen = 2;
    model.x = [0; 1; 2];
    model.y = [0; 0; 0];
    model.IEN = [1, 2; 2, 3];
    model.E = [100, 200];
    model.A = [1, 1];
    model.fixed_dof = [1, 2];
    model.fixed_val = [0, 0];
    model.force_dof = [5];
    model.force_val = [10];
    model.n_dof = model.nnp * model.ndof;
end

function model = preprocess_example2()
    model.title = '二维两杆桁架';
    model.nsd = 2;
    model.ndof = 2;
    model.nnp = 3;
    model.nel = 2;
    model.nen = 2;
    model.x = [1; 0; 1];
    model.y = [0; 0; 1];
    model.IEN = [1, 3; 2, 3];
    model.E = [1, 1];
    model.A = [1, 1];
    model.fixed_dof = [1,2,3,4];
    model.fixed_val = [0,0,0,0];
    model.force_dof = [5];
    model.force_val = [10];
    model.n_dof = model.nnp * model.ndof;
end

%% --------------------------- 2. 单元分析模块 --------------------------------
function [L, c, s, Ke] = truss2d_element_stiffness(x1, y1, x2, y2, E, A)
    dx = x2 - x1; dy = y2 - y1;
    L = sqrt(dx^2 + dy^2);
    if L < eps, error('单元长度为零'); end
    c = dx / L; s = dy / L;
    k = E * A / L;
    Ke = k * [c^2, c*s, -c^2, -c*s;
              c*s, s^2, -c*s, -s^2;
              -c^2, -c*s, c^2, c*s;
              -c*s, -s^2, c*s, s^2];
end

function [epsilon, sigma, N] = truss2d_element_stress(x1, y1, x2, y2, E, A, de)
    [L, c, s, ~] = truss2d_element_stiffness(x1, y1, x2, y2, E, A);
    B = 1/L * [-c, -s, c, s];
    de = de(:);
    epsilon = B * de;
    sigma = E * epsilon;
    N = sigma * A;
end

%% --------------------------- 3. 组装模块 -----------------------------------
function LM = generate_LM(IEN, ndof)
    nel = size(IEN, 1);
    nen = 2;
    LM = zeros(nel, nen * ndof);
    for e = 1:nel
        n1 = IEN(e, 1);
        n2 = IEN(e, 2);
        dof1 = (n1-1)*ndof + 1 : n1*ndof;
        dof2 = (n2-1)*ndof + 1 : n2*ndof;
        LM(e, :) = [dof1, dof2];
    end
end

function K_global = assemble_global_stiffness(n_dof, nel, Ke_list, LM)
    K_global = zeros(n_dof, n_dof);
    for e = 1:nel
        dof = LM(e, :);
        Ke = Ke_list{e};
        K_global(dof, dof) = K_global(dof, dof) + Ke;
    end
end

%% --------------------------- 4. 求解模块（处理奇异矩阵）--------------------
function [d, reaction] = solve_displacement_reaction(K_global, F, fixed_dof, fixed_val)
    n = size(K_global, 1);
    free_dof = setdiff(1:n, fixed_dof);
    K_ff = K_global(free_dof, free_dof);
    K_fe = K_global(free_dof, fixed_dof);
    F_f = F(free_dof);
    d_E = fixed_val(:);
    
    rhs = F_f - K_fe * d_E;
    
    % 检查矩阵是否奇异（条件数过大）
    if rcond(K_ff) < 1e-12
        % 奇异矩阵：使用伪逆（无警告，给出最小范数解）
        d_F = pinv(K_ff) * rhs;
    else
        % 非奇异：使用高效的反斜杠
        d_F = K_ff \ rhs;
    end
    
    d = zeros(n, 1);
    d(free_dof) = d_F;
    d(fixed_dof) = d_E;
    R = K_global * d - F;
    reaction = R(fixed_dof);
end

%% --------------------------- 5. 后处理模块 ---------------------------------
function postprocess(model, d, reactions, Ke_list, LM, fid)
    fprintf(fid, '\n========== 后处理结果 ==========\n');
    fprintf(fid, '节点位移 (u, v):\n');
    for i = 1:model.nnp
        u = d((i-1)*2+1);
        v = d((i-1)*2+2);
        fprintf(fid, '节点%d: u=%.6f, v=%.6f\n', i, u, v);
    end
    fprintf(fid, '\n约束反力 (自由度: 数值):\n');
    for j = 1:length(model.fixed_dof)
        dof = model.fixed_dof(j);
        val = reactions(j);
        fprintf(fid, '自由度%2d: %.6f\n', dof, val);
    end
    fprintf(fid, '\n单元计算结果:\n');
    fprintf(fid, '单元  长度        c        s         应变           应力            轴力\n');
    for e = 1:model.nel
        n1 = model.IEN(e,1); n2 = model.IEN(e,2);
        x1 = model.x(n1); y1 = model.y(n1);
        x2 = model.x(n2); y2 = model.y(n2);
        de = d(LM(e,:))';
        [L, c, s, ~] = truss2d_element_stiffness(x1, y1, x2, y2, model.E(e), model.A(e));
        [eps, sig, N] = truss2d_element_stress(x1, y1, x2, y2, model.E(e), model.A(e), de);
        fprintf(fid, '%2d   %8.4f %8.4f %8.4f %12.6e %12.6f %12.6f\n', ...
                e, L, c, s, eps, sig, N);
    end
end

%% --------------------------- 6. 性质验证函数 ---------------------------------
function check_stiffness_properties(K, title_str, fid)
    fprintf(fid, '\n--- %s 总体刚度矩阵性质 ---\n', title_str);
    if norm(K - K') < 1e-10
        fprintf(fid, '✓ 对称性验证通过。\n');
    else
        fprintf(fid, '✗ 不对称。\n');
    end
    if rank(K) < size(K,1)
        fprintf(fid, '✓ 矩阵奇异 (秩=%d < %d)，存在刚体位移。\n', rank(K), size(K,1));
    else
        fprintf(fid, '✗ 矩阵非奇异。\n');
    end
    if all(diag(K) >= -1e-10)
        fprintf(fid, '✓ 所有对角元非负。\n');
    else
        fprintf(fid, '✗ 存在负对角元。\n');
    end
    nz = nnz(K);
    total = numel(K);
    fprintf(fid, '稀疏性：非零元素 %.2f%%\n', nz/total*100);
end

%% --------------------------- 主程序：运行验证算例 ---------------------------
fprintf(fid, '================== 总体刚度矩阵组装与求解验证 ==================\n\n');

%% 算例1
fprintf(fid, '------------------ 算例1：一维两单元杆结构 ------------------\n');
model1 = preprocess_example1();
Ke_list1 = cell(model1.nel, 1);
for e = 1:model1.nel
    n1 = model1.IEN(e,1); n2 = model1.IEN(e,2);
    x1 = model1.x(n1); y1 = model1.y(n1);
    x2 = model1.x(n2); y2 = model1.y(n2);
    [~, ~, ~, Ke] = truss2d_element_stiffness(x1, y1, x2, y2, model1.E(e), model1.A(e));
    Ke_list1{e} = Ke;
end
LM1 = generate_LM(model1.IEN, model1.ndof);
K1 = assemble_global_stiffness(model1.n_dof, model1.nel, Ke_list1, LM1);
fprintf(fid, '总体刚度矩阵 K (6x6):\n');
for i = 1:6
    fprintf(fid, '  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n', K1(i,:));
end
check_stiffness_properties(K1, '算例1 (施加边界条件前)', fid);

F1 = zeros(model1.n_dof,1);
F1(model1.force_dof) = model1.force_val;
[d1, react1] = solve_displacement_reaction(K1, F1, model1.fixed_dof, model1.fixed_val);
postprocess(model1, d1, react1, Ke_list1, LM1, fid);
fprintf(fid, '\n理论位移: d2=0.1, d3=0.15; 实际: u2=%.6f, u3=%.6f\n', d1(3), d1(5));
fprintf(fid, '理论反力: 节点1水平反力 -10; 实际: %.6f\n', react1(1));

fprintf(fid, '求解策略: 非奇异用\\，奇异用pinv。\n');

%% 算例2
fprintf(fid, '\n------------------ 算例2：二维两杆桁架结构 ------------------\n');
model2 = preprocess_example2();
Ke_list2 = cell(model2.nel,1);
for e = 1:model2.nel
    n1 = model2.IEN(e,1); n2 = model2.IEN(e,2);
    x1 = model2.x(n1); y1 = model2.y(n1);
    x2 = model2.x(n2); y2 = model2.y(n2);
    [~, ~, ~, Ke] = truss2d_element_stiffness(x1, y1, x2, y2, model2.E(e), model2.A(e));
    Ke_list2{e} = Ke;
end
LM2 = generate_LM(model2.IEN, model2.ndof);
K2 = assemble_global_stiffness(model2.n_dof, model2.nel, Ke_list2, LM2);
fprintf(fid, '总体刚度矩阵 K (6x6):\n');
for i = 1:6
    fprintf(fid, '  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', K2(i,:));
end
check_stiffness_properties(K2, '算例2 (施加边界条件前)', fid);

F2 = zeros(model2.n_dof,1);
F2(model2.force_dof) = model2.force_val;
[d2, react2] = solve_displacement_reaction(K2, F2, model2.fixed_dof, model2.fixed_val);
postprocess(model2, d2, react2, Ke_list2, LM2, fid);

u3_theory = 38.284271; v3_theory = -10.000000;
fprintf(fid, '\n理论节点3位移: u=%.6f, v=%.6f\n', u3_theory, v3_theory);
fprintf(fid, '实际节点3位移: u=%.6f, v=%.6f\n', d2(5), d2(6));
fprintf(fid, '理论单元1应力: -10.000, 单元2应力: 14.142136\n');
fprintf(fid, '实际见上表，应一致。\n');

fprintf(fid, '求解策略: 非奇异用\\，奇异用pinv。\n');

%% 总体刚度矩阵第 j 列物理意义
fprintf(fid, '\n========== 总体刚度矩阵第 j 列物理意义 ==========\n');
j = 5;
e_j = zeros(model2.n_dof,1); e_j(j) = 1;
F_col = K2 * e_j;
fprintf(fid, '取第 %d 个自由度的单位位移向量 e_%d\n', j, j);
fprintf(fid, 'K2 * e_%d 得到第 %d 列：\n', j, j);
fprintf(fid, '  [%8.4f; %8.4f; %8.4f; %8.4f; %8.4f; %8.4f]\n', F_col);
fprintf(fid, '解释：总体刚度矩阵的第 %d 列表示当第 %d 个自由度产生单位位移（其余位移为零）时，\n', j, j);
fprintf(fid, '需要在所有自由度上施加的节点力向量。其中 k_{ij} 是第 j 个单位位移引起的第 i 个自由度方向的力。\n');

%% 附加题：三维桁架
fprintf(fid, '\n========== 附加题：三维空间桁架算例 ==========\n');
nodes_3d = [0,0,0; 1,0,0; 0,1,0; 0,0,1];
IEN_3d = [1,4; 2,4; 3,4];
E_3d = 1e6; A_3d = 0.01;
ndof_3d = 3; nnp_3d = size(nodes_3d,1);
n_dof_3d = nnp_3d * ndof_3d;
nel_3d = size(IEN_3d,1);

Ke_list_3d = cell(nel_3d,1);
for e = 1:nel_3d
    n1 = IEN_3d(e,1); n2 = IEN_3d(e,2);
    x1 = nodes_3d(n1,:); x2 = nodes_3d(n2,:);
    [~, ~, Ke] = truss3d_element_stiffness(x1, x2, E_3d, A_3d);
    Ke_list_3d{e} = Ke;
end

LM_3d = generate_LM(IEN_3d, ndof_3d);
K3 = assemble_global_stiffness(n_dof_3d, nel_3d, Ke_list_3d, LM_3d);
fixed_dof_3d = [1,2,3]; fixed_val_3d = [0,0,0];
F3 = zeros(n_dof_3d,1);
F3(12) = -1000;
[d3, R3] = solve_displacement_reaction(K3, F3, fixed_dof_3d, fixed_val_3d);
fprintf(fid, '三维桁架节点位移 (m):\n');
for i = 1:nnp_3d
    idx = (i-1)*3+1 : i*3;
    fprintf(fid, '节点%d: u=%.6f, v=%.6f, w=%.6f\n', i, d3(idx(1)), d3(idx(2)), d3(idx(3)));
end
fprintf(fid, '\n单元轴力 (N):\n');
for e = 1:nel_3d
    n1 = IEN_3d(e,1); n2 = IEN_3d(e,2);
    x1 = nodes_3d(n1,:); x2 = nodes_3d(n2,:);
    dof = LM_3d(e,:);
    de = d3(dof);
    [~, ~, N] = truss3d_element_stress(x1, x2, E_3d, A_3d, de);
    fprintf(fid, '单元%d (%d-%d): N = %.2f N\n', e, n1, n2, N);
end

fprintf(fid, '\n程序运行完毕。结果已保存至 results.txt\n');
fclose(fid);
fprintf(1, '程序运行完毕。请查看 results.txt 获取详细输出。\n');

%% --------------------------- 三维杆单元辅助函数 ---------------------------
function [L, dir_cos, Ke3D] = truss3d_element_stiffness(x1, x2, E, A)
    delta = x2 - x1;
    L = norm(delta);
    if L < eps, error('单元长度为零'); end
    dir_cos = delta / L;
    cx = dir_cos(1); cy = dir_cos(2); cz = dir_cos(3);
    k = E * A / L;
    C = [cx^2, cx*cy, cx*cz;
         cx*cy, cy^2, cy*cz;
         cx*cz, cy*cz, cz^2];
    Ke3D = k * [C, -C; -C, C];
end

function [epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de)
    [L, dir_cos, ~] = truss3d_element_stiffness(x1, x2, E, A);
    B = 1/L * [-dir_cos, dir_cos];
    de = de(:);
    epsilon = B * de;
    sigma = E * epsilon;
    N = sigma * A;
end