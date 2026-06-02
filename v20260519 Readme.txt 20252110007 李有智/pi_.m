clear; clc; close all;

% Polygon-based approximation of pi for a unit circle:
% L_n = n * 2 * sin(pi / n),  d = 2,  so  pi_n = L_n / d = n * sin(pi / n)
outDir = fileparts(mfilename('fullpath'));
n = 2.^(0:8).';
h = 1 ./ n;
piExact = pi;
piApprox = n .* sin(pi ./ n);
errApprox = abs(piExact - piApprox);
orderApprox = log(errApprox(2:end) ./ errApprox(1:end-1)) ./ log(h(2:end) ./ h(1:end-1));

[piWynnVals, idxWynn] = wynnEpsilonExtrapolated(piApprox);
piWynn = NaN(size(piApprox));
piWynn(idxWynn) = piWynnVals;
errWynn = abs(piExact - piWynn);
orderWynn = log(errWynn(idxWynn(2:end)) ./ errWynn(idxWynn(1:end-1))) ./ ...
    log(h(idxWynn(2:end)) ./ h(idxWynn(1:end-1)));

tbl = table(n, h, piApprox, errApprox, piWynn, errWynn, ...
    'VariableNames', {'n', 'h', 'pi_n', 'error_pi_n', 'pi_wynn', 'error_pi_wynn'});
writetable(tbl, fullfile(outDir, 'pi_results.csv'));
save(fullfile(outDir, 'pi_results.mat'), 'n', 'h', 'piApprox', 'errApprox', ...
    'orderApprox', 'piWynn', 'errWynn', 'orderWynn', 'tbl');

fid = fopen(fullfile(outDir, 'pi_results.txt'), 'w');
fprintf(fid, 'Polygon approximation of pi for a unit circle\n');
fprintf(fid, 'Formula: pi_n = n * sin(pi / n)\n\n');
fprintf(fid, '%6s %12s %22s %16s %22s %16s\n', ...
    'n', 'h', 'pi_n', 'error_pi_n', 'pi_wynn', 'error_pi_wynn');
for k = 1:numel(n)
    if ~isnan(piWynn(k))
        fprintf(fid, '%6d %12.6g %22.15f %16.8e %22.15f %16.8e\n', ...
            n(k), h(k), piApprox(k), errApprox(k), piWynn(k), errWynn(k));
    else
        fprintf(fid, '%6d %12.6g %22.15f %16.8e %22s %16s\n', ...
            n(k), h(k), piApprox(k), errApprox(k), 'NaN', 'NaN');
    end
end
fprintf(fid, '\nLocal order for pi_n      : ');
fprintf(fid, '%.2f ', orderApprox);
fprintf(fid, '\nLocal order for Wynn-eps  : ');
fprintf(fid, '%.2f ', orderWynn);
fprintf(fid, '\n');
fclose(fid);

fig = figure('Color', 'w', 'Position', [100, 100, 900, 560]);
ax = axes('Parent', fig);
loglog(ax, h, errApprox, '-v', ...
    'Color', [0.20, 0.45, 0.75], 'MarkerFaceColor', 'none', 'LineWidth', 1.4, 'MarkerSize', 7);
hold(ax, 'on');
loglog(ax, h(idxWynn), errWynn(idxWynn), '-^', ...
    'Color', [0.85, 0.40, 0.10], 'MarkerFaceColor', 'none', 'LineWidth', 1.4, 'MarkerSize', 7);
grid(ax, 'on');
box(ax, 'on');
xlim(ax, [1e-3, 1e0]);
ylim(ax, [1e-15, 1e5]);
xlabel(ax, 'h = 1 / n', 'FontSize', 12);
ylabel(ax, 'e_n', 'FontSize', 12);
title(ax, 'Convergence of Polygon Approximation for \pi', 'FontSize', 13);
text(ax, 1.1e-1, 2.5e-4, 'e_n = |\pi - \pi_n|', 'FontSize', 12, ...
    'Color', [0.10, 0.10, 0.10]);

for k = 1:numel(orderApprox)
    text(ax, h(k + 1) * 0.95, 1.2e4, sprintf('%.2f', orderApprox(k)), ...
        'Color', [0.10, 0.10, 0.10], 'FontSize', 11);
end

for k = 1:numel(orderWynn)
    text(ax, h(idxWynn(k + 1)) * 1.03, 3.0e-15, sprintf('%.2f', orderWynn(k)), ...
        'Color', [0.10, 0.10, 0.10], 'FontSize', 11);
end

text(ax, 1.6e-2, 2.0e1, sprintf('slope: %.2f', mean(orderApprox(end-3:end))), ...
    'Color', [0.20, 0.10, 0.10], 'FontSize', 12);
text(ax, 1.0e-2, 8.0e-13, sprintf('slope: %.2f', orderWynn(end)), ...
    'Color', [0.20, 0.10, 0.10], 'FontSize', 12);

exportgraphics(fig, fullfile(outDir, 'pi_convergence.png'), 'Resolution', 300);
savefig(fig, fullfile(outDir, 'pi_convergence.fig'));

disp(tbl);
fprintf('\nLocal order for direct approximation: ');
fprintf('%.2f ', orderApprox);
fprintf('\nLocal order for Wynn-epsilon: ');
fprintf('%.2f ', orderWynn);
fprintf('\nBest direct pi approximation       : %.15f (n = %d)\n', piApprox(end), n(end));
fprintf('Best Wynn-epsilon approximation    : %.15f (n = %d)\n', ...
    piWynn(idxWynn(end)), n(idxWynn(end)));

function [sHat, idxHat] = wynnEpsilonExtrapolated(s)
    s = s(:);
    m = numel(s);
    E = NaN(m + 2, m + 2);
    E(:, 1) = 0;
    E(2:m + 1, 2) = s;
    nHat = floor((m - 1) / 2);
    sHat = NaN(nHat, 1);
    idxHat = NaN(nHat, 1);
    count = 0;

    for col = 3:(m + 2)
        for row = 2:(m - col + 3)
            delta = E(row + 1, col - 1) - E(row, col - 1);
            if abs(delta) < eps
                E(row, col) = NaN;
            else
                E(row, col) = E(row + 1, col - 2) + 1 / delta;
            end
        end

        if mod(col, 2) == 0
            count = count + 1;
            sHat(count) = E(2, col);
            idxHat(count) = 2 * count + 1;
        end
    end
end
