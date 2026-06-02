%% ========================================================================
% 三维杆单元（空间桁架单元）刚度矩阵与应力计算程序
% 作者：根据作业要求编写
% 功能：
%   1. 计算单元长度、方向余弦、全局刚度矩阵
%   2. 根据节点位移计算应变、应力、轴力
%   3. 验证刚度矩阵性质（对称、奇异、半正定）
%   4. 演示刚度矩阵各列的物理意义
% ========================================================================

%% 清理环境
clear; clc; close all;

%% ----------------------------- 任务2：核心函数 ---------------------------------
% 函数定义在文件末尾（或单独保存为 .m 文件）

%% ----------------------------- 任务3：算例验证 ---------------------------------
fprintf('================== 三维杆单元程序验证 ==================\n\n');

%% 算例1：沿 x 轴的一维杆单元
fprintf('------------------ 算例1：一维杆单元（沿x轴）------------------\n');
x1 = [0, 0, 0];
x2 = [2, 0, 0];
E = 200e9;      % 200 GPa
A = 1.0e-4;     % m^2
de1 = [0, 0, 0, 1.0e-3, 0, 0]';  % 节点位移，单位 m

% 计算刚度矩阵
[L, dir_cos, Ke1] = truss3d_element_stiffness(x1, x2, E, A);
fprintf('单元长度 L = %.6f m (理论: 2.000 m)\n', L);
fprintf('方向余弦 (cx,cy,cz) = (%.4f, %.4f, %.4f) (理论: 1,0,0)\n', dir_cos(1), dir_cos(2), dir_cos(3));

% 显示刚度矩阵（简化显示非零行）
fprintf('刚度矩阵 Ke (6x6)：\n');
disp(Ke1);
% 检查是否只与 x 方向自由度有关（非零元素仅出现在 [1,1], [1,4], [4,1], [4,4]）
idx_nonzero = find(abs(Ke1) > 1e-10);
[row,col] = ind2sub([6,6], idx_nonzero);
fprintf('非零元素位置 (行,列):\n');
disp(unique([row,col], 'rows'));

% 计算应力
[epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de1);
fprintf('轴向应变 epsilon = %.4e (理论: 5.000e-4)\n', epsilon);
fprintf('轴向应力 sigma  = %.2f MPa (理论: 100.00 MPa)\n', sigma/1e6);
fprintf('轴力 N         = %.2f N (理论: 10000.00 N)\n', N);

% 算例1性质验证
fprintf('\n--- 算例1刚度矩阵性质验证 ---\n');
check_stiffness_properties(Ke1);

%% 算例2：空间任意方向杆单元
fprintf('\n------------------ 算例2：空间任意方向杆单元 ------------------\n');
x1 = [0, 0, 0];
x2 = [1, 2, 2];
E = 210e9;      % 210 GPa
A = 2.0e-4;     % m^2
de2 = [0, 0, 0, 1.0e-3, 2.0e-3, 2.0e-3]';  % 节点位移

[L, dir_cos, Ke2] = truss3d_element_stiffness(x1, x2, E, A);
fprintf('单元长度 L = %.6f m (理论: 3.000 m)\n', L);
fprintf('方向余弦 (cx,cy,cz) = (%.4f, %.4f, %.4f) (理论: 0.3333,0.6667,0.6667)\n', ...
    dir_cos(1), dir_cos(2), dir_cos(3));

% 计算应力
[epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de2);
fprintf('轴向应变 epsilon = %.4e (理论: 1.000e-3)\n', epsilon);
fprintf('轴向应力 sigma  = %.2f MPa (理论: 210.00 MPa)\n', sigma/1e6);
fprintf('轴力 N         = %.2f N (理论: 42000.00 N)\n', N);

% 性质验证
fprintf('\n--- 算例2刚度矩阵性质验证 ---\n');
check_stiffness_properties(Ke2);

% 额外验证：刚体平移位移不产生内力
de_rigid = [0.1, 0.2, 0.3, 0.1, 0.2, 0.3]';  % 刚体平移（两端位移相同）
Fe_rigid = Ke2 * de_rigid;
fprintf('\n刚体平移位移产生的节点力列阵（应全为零）:\n');
disp(Fe_rigid);
if norm(Fe_rigid) < 1e-10
    fprintf('✓ 刚体平移位移不产生内力，验证通过。\n');
else
    fprintf('✗ 刚体平移位移产生了非零内力，请检查。\n');
end

%% ----------------------------- 任务4：刚度矩阵物理意义验证 -----------------------------
fprintf('\n================== 任务4：刚度矩阵物理意义验证 ==================\n');
% 任选一个自由度 j，例如 j = 3 (w1, 即节点1的z方向位移)
j = 3;
de_unit = zeros(6,1);
de_unit(j) = 1;   % 第j个自由度位移为1，其余为0
Fe = Ke2 * de_unit;   % 对应刚度矩阵的第j列

fprintf('取自由度 j = %d (节点1的z方向位移) 为单位位移，其他位移为0\n', j);
fprintf('得到的节点力列阵 Fe = Ke * e_j 为:\n');
disp(Fe);
fprintf('该向量恰好是刚度矩阵的第 %d 列:\n', j);
disp(Ke2(:, j));
% 解释
fprintf('\n解释：单元刚度矩阵 Ke 的第 j 列元素的物理意义是：\n');
fprintf('当第 j 个自由度产生单位位移（其余自由度位移为 0）时，\n');
fprintf('在各个自由度方向上所需要施加的节点力。\n');
fprintf('即 k_{ij} 表示：第 j 个自由度单位位移引起的第 i 个自由度方向的力。\n');

%% ========================================================================
% 函数定义部分
% ========================================================================

%% 函数1：计算三维杆单元刚度矩阵
function [L, dir_cos, Ke] = truss3d_element_stiffness(x1, x2, E, A)
    % 输入：
    %   x1, x2 : 节点坐标，1x3向量 [x,y,z]
    %   E      : 弹性模量 (Pa)
    %   A      : 横截面积 (m^2)
    % 输出：
    %   L            : 单元长度 (m)
    %   dir_cos      : 方向余弦 [cx, cy, cz]
    %   Ke           : 6x6全局坐标系下的单元刚度矩阵 (N/m)
    
    % 检查退化单元（两点重合）
    if norm(x2 - x1) < eps
        error('错误：两个节点坐标重合，无法构成杆单元。');
    end
    
    % 计算长度和方向余弦
    delta = x2 - x1;
    L = norm(delta);
    dir_cos = delta / L;   % [cx, cy, cz]
    cx = dir_cos(1); cy = dir_cos(2); cz = dir_cos(3);
    
    % 计算刚度矩阵系数
    k = E * A / L;   % 轴向刚度系数
    % 构造方向余弦的乘积矩阵
    C = [cx^2, cx*cy, cx*cz;
         cx*cy, cy^2, cy*cz;
         cx*cz, cy*cz, cz^2];
    
    % 全局刚度矩阵（6x6，分块形式）
    Ke = k * [ C, -C;
              -C,  C];
    % 注意：MATLAB中矩阵分块直接拼接
end

%% 函数2：计算应变、应力、轴力
function [epsilon, sigma, N] = truss3d_element_stress(x1, x2, E, A, de)
    % 输入：
    %   x1,x2 : 节点坐标
    %   E,A   : 材料属性
    %   de    : 单元节点位移向量 [u1,v1,w1,u2,v2,w2]' (m)
    % 输出：
    %   epsilon : 轴向应变 (无量纲)
    %   sigma   : 轴向应力 (Pa)
    %   N       : 轴力 (N)，受拉为正
    
    % 计算单元长度和方向余弦
    delta = x2 - x1;
    L = norm(delta);
    if L < eps
        error('错误：两个节点坐标重合，无法计算应变。');
    end
    dir_cos = delta / L;
    
    % 构造应变-位移矩阵 B (1x6)
    % B = [-cx/L, -cy/L, -cz/L, cx/L, cy/L, cz/L]? 注意：杆单元应变 = (u2-u1)/L * 方向余弦点积
    % 实际公式：epsilon = ( (u2-u1)*cx + (v2-v1)*cy + (w2-w1)*cz ) / L
    % 等价于 epsilon = [-dir_cos, dir_cos] / L * de
    % 但更常见的 B = 1/L * [-cx, -cy, -cz, cx, cy, cz]
    B = 1/L * [-dir_cos, dir_cos];
    
    % 轴向应变
    epsilon = B * de;   % 标量
    
    % 轴向应力 (胡克定律)
    sigma = E * epsilon;
    
    % 轴力 = 应力 * 面积
    N = sigma * A;
end

%% 辅助函数：检查刚度矩阵性质（对称、奇异、特征值非负）
function check_stiffness_properties(Ke)
    % 对称性检查
    if norm(Ke - Ke') < 1e-10
        fprintf('✓ 刚度矩阵对称性验证通过。\n');
    else
        fprintf('✗ 刚度矩阵不对称。\n');
    end
    
    % 奇异性检查：行列式是否为零（或接近零）
    det_val = det(Ke);
    if abs(det_val) < 1e-10
        fprintf('✓ 刚度矩阵行列式 ≈ 0 (det = %.2e)，奇异矩阵，符合预期。\n', det_val);
    else
        fprintf('✗ 刚度矩阵行列式不为零 (det = %.2e)，非奇异，请检查。\n', det_val);
    end
    
    % 特征值非负（半正定性）
    eigvals = eig(Ke);
    fprintf('刚度矩阵特征值：\n');
    disp(eigvals);
    if all(eigvals >= -1e-10)  % 允许微小负值（数值误差）
        fprintf('✓ 所有特征值非负（半正定），符合预期。\n');
        % 检查零特征值的个数（刚体模式，三维杆单元有1个零特征值）
        zero_count = sum(abs(eigvals) < 1e-10);
        fprintf('  零特征值个数 = %d (理论应为 1，对应轴向刚体位移)\n', zero_count);
    else
        fprintf('✗ 存在负特征值，不是半正定矩阵。\n');
    end
end