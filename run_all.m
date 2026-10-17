function run_all()
%RUN_ALL  Reproduce every figure in the paper.
%
%   Each figN script begins with "clear". Because MATLAB's run() executes a
%   script in the CALLER's workspace, calling run() directly from a loop wipes
%   the loop variables themselves. Each script is therefore launched inside the
%   isolated workspace of the local function runOne() below, where the clear is
%   harmless.
%
%   Expect a long run: the full set is several hours on a laptop. Run the
%   scripts individually while you are still checking things.

scripts = {'fig1_nmse_floor','fig2_training_amplitude','fig3_antenna_scaling', ...
           'fig4_gain_mismatch','fig5_kappa_tracking','fig6_probe_design', ...
           'fig7_saleh','fig8_element_spread','fig9_scaling'};

for i = 1:numel(scripts)
    fprintf('\n===== %s =====\n', scripts{i});
    t = tic;
    try
        runOne(scripts{i});
        fprintf('(%.1f s)\n', toc(t));
    catch ME
        fprintf(2, 'FAILED after %.1f s: %s\n', toc(t), ME.message);
    end
end
end

% -------------------------------------------------------------------------
function runOne(name)
% Isolated workspace: the script's "clear" only affects this function.
run(name);
end
