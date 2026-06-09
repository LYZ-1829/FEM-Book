%% ========================================================================
% 三维杆单元有限元程序（修正版 + 文件输出）
% 包含：单元刚度矩阵、应力计算、性质验证、物理意义、空间桁架组装求解
% 输出：屏幕显示 + results.txt
% ========================================================================

clear; clc; close all;

% 打开结果文件
fid = fopen('results.txt', 'w');
if fid == -1
    error('无法创建 results.txt 文件');
end

%% ----------------------------- 任务3：验证算例 ---------------------------------
fprintf(1, '================== 三维杆单元程序验证 ==================\n\n');
fprintf(fid, '================== 三维杆单元程序验证 ==================\n\n');

% 算例1：一维杆单元 (沿x轴)
fprintf(1, '------------------ 算例1：一维杆单元 ------------------\n');
fprintf(fid, '------------------ 算例1：一维杆单元 ------------------\n');
x1 = [0,0,0]; x2 = [2,0,0];
E1 = 200e9; A1 = 1e-4;
de1 = [0,0,0, 1e-3,0,0]';

[L1, dc1, Ke1] = truss3d_element_stiffness(x1, x2, E1, A1);
fprintf(1, '单元长度 L = %.6f m (理论 2.000 m)\n', L1);
fprintf(fid, '单元长度 L = %.6f m (理论 2.000 m)\n', L1);
fprintf(1, '方向余弦 = (%.4f, %.4f, %.4f) (理论 1,0,0)\n', dc1(1), dc1(2), dc1(3));
fprintf(fid, '方向余弦 = (%.4f, %.4f, %.4f) (理论 1,0,0)\n', dc1(1), dc1(2), dc1(3));
fprintf(1, '刚度矩阵 Ke1 (仅显示非零):\n');
fprintf(fid, '刚度矩阵 Ke1 (仅显示非零):\n');
print_matrix(fid, Ke1, 'Ke1');

[eps1, sig1, N1] = truss3d_element_stress(x1, x2, E1, A1, de1);
fprintf(1, '应变 = %.4e (理论 5.000e-4)\n', eps1);
fprintf(fid, '应变 = %.4e (理论 5.000e-4)\n', eps1);
fprintf(1, '应力 = %.2f MPa (理论 100.00 MPa)\n', sig1/1e6);
fprintf(fid, '应力 = %.2f MPa (理论 100.00 MPa)\n', sig1/1e6);
fprintf(1, '轴力 = %.2f N (理论 10000.00 N)\n', N1);
fprintf(fid, '轴力 = %.2f N (理论 10000.00 N)\n', N1);

% 算例2：空间任意方向杆单元
fprintf(1, '\n------------------ 算例2：空间任意方向杆单元 ------------------\n');
fprintf(fid, '\n------------------ 算例2：空间任意方向杆单元 ------------------\n');
x1 = [0,0,0]; x2 = [1,2,2];
E2 = 210e9; A2 = 2e-4;
de2 = [0,0,0, 1e-3,2e-3,2e-3]';

[L2, dc2, Ke2] = truss3d_element_stiffness(x1, x2, E2, A2);
fprintf(1, '单元长度 L = %.6f m (理论 3.000 m)\n', L2);
fprintf(fid, '单元长度 L = %.6f m (理论 3.000 m)\n', L2);
fprintf(1, '方向余弦 = (%.4f, %.4f, %.4f) (理论 1/3,2/3,2/3)\n', dc2(1), dc2(2), dc2(3));
fprintf(fid, '方向余弦 = (%.4f, %.4f, %.4f) (理论 1/3,2/3,2/3)\n', dc2(1), dc2(2), dc2(3));

[eps2, sig2, N2] = truss3d_element_stress(x1, x2, E2, A2, de2);
fprintf(1, '应变 = %.4e (理论 1.000e-3)\n', eps2);
fprintf(fid, '应变 = %.4e (理论 1.000e-3)\n', eps2);
fprintf(1, '应力 = %.2f MPa (理论 210.00 MPa)\n', sig2/1e6);
fprintf(fid, '应力 = %.2f MPa (理论 210.00 MPa)\n', sig2/1e6);
fprintf(1, '轴力 = %.2f N (理论 42000.00 N)\n', N2);
fprintf(fid, '轴力 = %.2f N (理论 42000.00 N)\n', N2);

%% ----------------------------- 任务3&4：刚度矩阵性质验证 -----------------------------
fprintf(1, '\n========== 刚度矩阵性质验证 (使用算例2的 Ke2) ==========\n');
fprintf(fid, '\n========== 刚度矩阵性质验证 (使用算例2的 Ke2) ==========\n');
% 对称性
if norm(Ke2 - Ke2') < 1e-10
    fprintf(1, '✓ 对称性验证通过。\n');
    fprintf(fid, '✓ 对称性验证通过。\n');
else
    fprintf(1, '✗ 不对称。\n');
    fprintf(fid, '✗ 不对称。\n');
end
% 奇异性
detVal = det(Ke2);
fprintf(1, '行列式 = %.2e (应近似为0，奇异矩阵)\n', detVal);
fprintf(fid, '行列式 = %.2e (应近似为0，奇异矩阵)\n', detVal);
% 特征值
eigVal = eig(Ke2);
fprintf(1, '特征值：\n');
fprintf(fid, '特征值：\n');
for i = 1:length(eigVal)
    fprintf(1, '  %12.6e\n', eigVal(i));
    fprintf(fid, '  %12.6e\n', eigVal(i));
end
if all(eigVal >= -1e-10)
    fprintf(1, '✓ 所有特征值非负，半正定。\n');
    fprintf(fid, '✓ 所有特征值非负，半正定。\n');
    zeroCount = sum(abs(eigVal) < 1e-10);
    fprintf(1, '  零特征值个数 = %d (理论为1，对应刚体位移模式)\n', zeroCount);
    fprintf(fid, '  零特征值个数 = %d (理论为1，对应刚体位移模式)\n', zeroCount);
end

% 刚体平移位移不产生内力
de_rigid = [0.1, 0.2, 0.3, 0.1, 0.2, 0.3]';
Fe_rigid = Ke2 * de_rigid;
fprintf(1, '刚体平移产生的节点力 (应全为零):\n');
fprintf(fid, '刚体平移产生的节点力 (应全为零):\n');
for i = 1:length(Fe_rigid)
    fprintf(1, '  %12.6e\n', Fe_rigid(i));
    fprintf(fid, '  %12.6e\n', Fe_rigid(i));
end
if norm(Fe_rigid) < 1e-10
    fprintf(1, '✓ 刚体平移不产生内力，验证通过。\n');
    fprintf(fid, '✓ 刚体平移不产生内力，验证通过。\n');
end

% 物理意义：取第3个自由度 (w1) 为单位位移
fprintf(1, '\n========== 任务4：刚度矩阵物理意义验证 ==========\n');
fprintf(fid, '\n========== 任务4：刚度矩阵物理意义验证 ==========\n');
j = 3;
de_unit = zeros(6,1); de_unit(j)=1;
Fe = Ke2 * de_unit;
fprintf(1, '自由度 j=%d (节点1的z向位移) 为单位位移时，节点力列阵为 Ke 的第 %d 列：\n', j, j);
fprintf(fid, '自由度 j=%d (节点1的z向位移) 为单位位移时，节点力列阵为 Ke 的第 %d 列：\n', j, j);
for i = 1:length(Fe)
    fprintf(1, '  %12.6e\n', Fe(i));
    fprintf(fid, '  %12.6e\n', Fe(i));
end
fprintf(1, 'Ke 的第 %d 列：\n', j);
fprintf(fid, 'Ke 的第 %d 列：\n', j);
for i = 1:6
    fprintf(1, '  %12.6e\n', Ke2(i,j));
    fprintf(fid, '  %12.6e\n', Ke2(i,j));
end
fprintf(1, '解释：k_{ij} 表示第 j 个自由度产生单位位移时，在第 i 个自由度方向需施加的力。\n');
fprintf(fid, '解释：k_{ij} 表示第 j 个自由度产生单位位移时，在第 i 个自由度方向需施加的力。\n');

%% ----------------------------- 附加题：组装整体刚度矩阵并求解空间桁架 -----------------------------
fprintf(1, '\n========== 附加题：简单空间桁架求解 ==========\n');
fprintf(fid, '\n========== 附加题：简单空间桁架求解 ==========\n');
% 定义桁架结构：4个节点，5个杆件（三维四角锥桁架）
nodes = [0,0,0;   % 节点1
         1,0,0;   % 节点2
         0,1,0;   % 节点3
         0,0,1];  % 节点4
elements = [1,2, 1e-4, 200e9;
            1,3, 1e-4, 200e9;
            1,4, 1e-4, 200e9;
            2,3, 1e-4, 200e9;
            3,4, 1e-4, 200e9];
fixed_dofs = [1,2,3];  % 节点1固定
F = zeros(12,1);
F(12) = -1000;    % 节点4的z方向

% 组装整体刚度矩阵
nNodes = size(nodes,1);
nDof = nNodes * 3;
K_global = zeros(nDof, nDof);
for iElem = 1:size(elements,1)
    n1 = elements(iElem,1);
    n2 = elements(iElem,2);
    A_e = elements(iElem,3);
    E_e = elements(iElem,4);
    x1 = nodes(n1,:);
    x2 = nodes(n2,:);
    [~, ~, Ke] = truss3d_element_stiffness(x1, x2, E_e, A_e);
    dof = [ (n1-1)*3+1 : (n1-1)*3+3, (n2-1)*3+1 : (n2-1)*3+3 ];
    K_global(dof, dof) = K_global(dof, dof) + Ke;
end

% 求解位移
free_dofs = setdiff(1:nDof, fixed_dofs);
K_ff = K_global(free_dofs, free_dofs);
F_f = F(free_dofs);
U_f = K_ff \ F_f;
U = zeros(nDof,1);
U(free_dofs) = U_f;

% 显示结果
fprintf(1, '节点位移 (m)：\n');
fprintf(fid, '节点位移 (m)：\n');
for i = 1:nNodes
    fprintf(1, '节点%d: u=%.6f, v=%.6f, w=%.6f\n', i, U((i-1)*3+1), U((i-1)*3+2), U((i-1)*3+3));
    fprintf(fid, '节点%d: u=%.6f, v=%.6f, w=%.6f\n', i, U((i-1)*3+1), U((i-1)*3+2), U((i-1)*3+3));
end

% 计算各杆件内力
fprintf(1, '\n杆件内力 (N，拉力为正)：\n');
fprintf(fid, '\n杆件内力 (N，拉力为正)：\n');
for iElem = 1:size(elements,1)
    n1 = elements(iElem,1);
    n2 = elements(iElem,2);
    A_e = elements(iElem,3);
    E_e = elements(iElem,4);
    x1 = nodes(n1,:);
    x2 = nodes(n2,:);
    dof1 = (n1-1)*3+1 : (n1-1)*3+3;
    dof2 = (n2-1)*3+1 : (n2-1)*3+3;
    de = [U(dof1); U(dof2)];
    [~, ~, N] = truss3d_element_stress(x1, x2, E_e, A_e, de);
    fprintf(1, '杆件%d (%d-%d): N = %.2f N\n', iElem, n1, n2, N);
    fprintf(fid, '杆件%d (%d-%d): N = %.2f N\n', iElem, n1, n2, N);
end

fprintf(1, '\n程序运行完毕。结果已保存至 results.txt\n');
fprintf(fid, '\n程序运行完毕。结果已保存至 results.txt\n');
fclose(fid);

%% ----------------------------- 辅助函数定义（放在文件末尾）-------------------------
function [L, dir_cos, Ke] = truss3d_element_stiffness(x1, x2, E, A)
    delta = x2 - x1;
    L = norm(delta);
    if L < eps
        error('错误：两个节点坐标重合，无法构成杆单元。');
    end
    dir_cos = delta / L;
    cx = dir_cos(1); cy = dir_cos(2); cz = dir_cos(3);
    k = E * A / L;
    C = [cx^2, cx*cy, cx*cz;
         cx*cy, cy^2, cy*cz;
         cx*cz, cy*cz, cz^2];
    Ke = k * [ C, -C; -C, C];
end

function [epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de)
    [L, dir_cos, ~] = truss3d_element_stiffness(x1, x2, E, A);
    B = 1/L * [-dir_cos, dir_cos];
    de = de(:);
    epsilon = B * de;
    sigma = E * epsilon;
    N = sigma * A;
end

function print_matrix(fid, mat, name)
    fprintf(1, '%s =\n', name);
    fprintf(fid, '%s =\n', name);
    for i = 1:size(mat,1)
        fprintf(1, '  ');
        fprintf(fid, '  ');
        for j = 1:size(mat,2)
            fprintf(1, '%12.6e ', mat(i,j));
            fprintf(fid, '%12.6e ', mat(i,j));
        end
        fprintf(1, '\n');
        fprintf(fid, '\n');
    end
end