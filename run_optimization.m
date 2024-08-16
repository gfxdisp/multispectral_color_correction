% define the wavelength range for calculation 
low = 390;
high = 780;
dim = high - low + 1;
lambda = low:high;
lambda2 = low:2:high;
lambda4 = low:4:high;
lambda5 = low:5:high;
lambda10 = low:10:high;

% define the color matching function
% 360-830 step=1
cmf = load('xyz1931.mat');
cmf = crop_spectrum(cmf, 'xyz', low, high);

% PMCC reflectance self-measured
% 380-780 step=1
pmcc = load('pmcc.mat');
pmcc = crop_spectrum(pmcc, 'reflectance', low, high);

% SFU reflectance
% 380-780 step=4
sfu = load('SFU_1993.mat');
sfu = crop_spectrum(sfu, 'all', 380, 780);
sfu = interp_spectrum(sfu, 380:4:780, 380:780);
sfu = sfu(:, low-380+1:high-380+1);

% define camera
cam = 'A7R3';
if strcmp(cam, 'A7R3')
    % define the self-calculated Sony A7R3 SSF 
    % 380-780 step=2
    camera = load('SSF_240506/A7R3_390_780_a0.001_r40_iter1.mat');
    ssf = crop_spectrum(camera, 'rgb', low, high);
    ssf = interp_spectrum(ssf, lambda2, lambda);
    % define CCM weight regularization
    beta = 1.0; % [0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
elseif strcmp(cam, 'IDS')
    % define the self-calculated IDS SSF 
    % 380-780 step=2
    camera = load('SSF_240506/IDS_390_780_a0.001_r40_iter5.mat');
    ssf = crop_spectrum(camera, 'rgb', low, high);
    ssf = interp_spectrum(ssf, lambda2, lambda);
    % define CCM weight regularization
    beta = 0.2; % [0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
else
    disp('wrong camera')
    pause
end

% define the illuminant
% 300-780 step=1
illum_std = load('illum_D65.mat');
illum_std = crop_spectrum(illum_std, 'power', low, high);

% define the LED channel response
% 350-1000 step=1
path = 'measurement_240515';
channels = 18;
samples = 20;
spectralon_ratio = 0.91;
responses = zeros(channels, 651);
for i = 1:channels
    for j = 1:samples
        filename = fullfile(path, 'channels', sprintf('%d/%02d.mat', i, j-1));
        tmp = load(filename);
        responses(i, :) = responses(i, :) + tmp.L ;
    end
end

bases = 683.002 * responses(:, low-350+1:high-350+1) / samples / spectralon_ratio;
w0 = [0.4759, 0.5237, 0.0000, 0.2429, 0.2514, 0.5726, 1.0000, 0.0000, 0.5239, 0.4102, 0.2986, 0.2335, 0.2552, 0.4199, 0.3793, 0.2984, 0.5311, 0.3095];
illum_apr = w0 * bases;

% Make standard illuminant the same luminance as the approximated illuminant
s2lum = sprad2color('lum1924');
illum_std = illum_std * (s2lum.get_color(lambda, illum_apr) / s2lum.get_color(lambda, illum_std));

% XYZ values from specbos measurement under the approximated illuminant
spec_path = fullfile(path, 'd65_specbos');
specbos_xyz_int = zeros(30, 3);
specbos_xyz_val = zeros(30, 3);
for i = 1:30
    filename = fullfile(spec_path, sprintf('p%02d_lambda_L.mat', i));
    tmp = load(filename);
    specbos_xyz_val(i, 1) = tmp.x / tmp.y * tmp.Y;
    specbos_xyz_val(i, 2) = tmp.Y;
    specbos_xyz_val(i, 3) = (1 - tmp.x - tmp.y) * tmp.Y / tmp.y;
    L = crop_lambda(tmp, 'L', low, high);
    s2lum.get_color(lambda, L);
    specbos_xyz_int(i, :) = 683.002 * L * cmf';
end

% Luther condition CCM to match camera SSF to CMF: ssf' * M = cmf'
M_lu = ssf' \ cmf';

% LCC from integrated camera RGBs to CMF XYZs under the approximated illuminant
wp_idx = 25; % white point index
RGBs = pmcc * (illum_apr' .* ssf');
XYZs = pmcc * (illum_apr' .* cmf');
M_apr_int = RGBs \ XYZs;

apr_cmf = illum_apr' .* cmf';
apr_ssf = illum_apr' .* ssf';
spec_apr_lu = apr_ssf * M_lu;
spec_apr_int = apr_ssf * M_apr_int;

% define how many illuminants for optimization
k = 1;

% define the noise param and SNR weight in loss function
noise = noise_param(); 
gamma = 0.1; % [0, 0.001, 0.002, 0.003, 0.006, 0.01, 0.02, 0.03, 0.06, 0.1, 0.2, 0.3, 0.6, 1.0, 2.0, 3.0, 6.0]

% define random seeds for running the optimization for multiple times
seeds = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
n_s = length(seeds);

ll_SNR  = -inf * ones(1, n_s);
ll_data = inf * ones(1, n_s);
ll_fval = inf * ones(1, n_s);
ll_M    = cell(1, n_s);
ll_w    = cell(1, n_s);

for ss = 1:n_s
    tic
    rng(seeds(ss));
    
    loss = @(x) spec_match_loss(x, k, bases', ssf', apr_cmf, beta, gamma, noise);
    un = channels + 9;
    x0 = zeros(1, un*k);
    for i = 1:k
        x0((i-1)*channels+1 : i*channels) = rand(1, channels);
        x0(channels*k+(i-1)*9+1 : channels*k+i*9) = rand(3, 3);
    end
    lb = [repmat(zeros(1, channels), 1, k), repmat(-ones(1, 9) * Inf, 1, k)];
    ub = [repmat(ones(1, channels), 1 ,k), repmat(ones(1, 9) * Inf, 1, k)];
    
    options = optimoptions('fmincon', 'UseParallel', true, 'MaxFunctionEvaluations', 500000, 'MaxIterations', 5000);
    [res, fval] = fmincon(loss, x0, [], [], [], [], lb, ub, [], options);

    [data, E_SNR] = spec_match(res, k, bases', ssf', apr_cmf, noise);
    w_opt = reshape(res(1:channels*k), [1, channels, k]);
    M_opt = reshape(permute(reshape(res((channels*k+1):un*k), [3, 3, k]), [1, 3, 2]), [3*k, 3]);
    
    toc
    ll_fval(ss) = fval;
    ll_data(ss) = data;
    ll_SNR(ss) = E_SNR;
    ll_M{ss} = M_opt;
    ll_w{ss} = w_opt;
end

% choose the optimal solutions from different seeds
[~, idx_min] = min(ll_fval);
w_opt = ll_w{idx_min};
M_opt = ll_M{idx_min};
opt = sum(w_opt .* bases', 2); % dim x 1 x k
opt_ssf = reshape(opt .* ssf', [dim, 3*k]);
spec_opt = opt_ssf * M_opt;
RGBs = pmcc * opt_ssf;
XYZs = pmcc * apr_cmf;
M_opt_int = RGBs \ XYZs;


de_apr_sim_lu = calculate_de(pmcc * spec_apr_lu, specbos_xyz_int, wp_idx);
fprintf('PMCC simulation CIEDE2000 under Approx illuminant using Luther matrix \n');
fprintf('mean: %5.3f\n', mean(de_apr_sim_lu));
fprintf('max: %5.3f\n\n', max(de_apr_sim_lu));

de_apr_sim_int = calculate_de(pmcc * spec_apr_int, specbos_xyz_int, wp_idx);
fprintf('PMCC simulation CIEDE2000 under Approx illuminant using fitted matrix from integration \n');
fprintf('mean: %5.3f\n', mean(de_apr_sim_int));
fprintf('max: %5.3f\n\n', max(de_apr_sim_int));

de_opt_sim = calculate_de(pmcc * spec_opt, specbos_xyz_int, wp_idx);
fprintf('PMCC simulation CIEDE2000 under optimized illuminant using optimized matrix \n');
fprintf('mean: %5.3f\n', mean(de_opt_sim));
fprintf('max: %5.3f\n\n', max(de_opt_sim));

de_opt_sim_int = calculate_de(RGBs * M_opt_int, specbos_xyz_int, wp_idx);
fprintf('PMCC simulation CIEDE2000 under optimized illuminant using fitted matrix from integration \n');
fprintf('mean: %5.3f\n', mean(de_opt_sim_int));
fprintf('max: %5.3f\n\n', max(de_opt_sim_int));
    

% SFU evaluation
tmp = sum(sfu, 2);
[~, sfu_wp_idx] = max(tmp); % find the white point

sfu_XYZs_tgt = sfu * apr_cmf;
sfu_RGBs_apr_lu = sfu * spec_apr_lu;
sfu_RGBs_apr_int = sfu * spec_apr_int;
sfu_RGBs_opt = sfu * spec_opt;
sfu_RGBs_opt_int = sfu * opt_ssf * M_opt_int;

sfu_de_apr_lu = calculate_de(sfu_RGBs_apr_lu, sfu_XYZs_tgt, sfu_wp_idx);
fprintf('SFU simulation CIEDE2000 under Approx illuminant using Luther matrix \n');
fprintf('mean: %5.3f\n', mean(sfu_de_apr_lu));
fprintf('max: %5.3f\n', max(sfu_de_apr_lu));
fprintf('95 %%: %5.3f\n\n', get95(sfu_de_apr_lu));

sfu_de_apr_int = calculate_de(sfu_RGBs_apr_int, sfu_XYZs_tgt, sfu_wp_idx);
fprintf('SFU simulation CIEDE2000 under Approx illuminant using fitted matrix from integration \n');
fprintf('mean: %5.3f\n', mean(sfu_de_apr_int));
fprintf('max: %5.3f\n', max(sfu_de_apr_int));
fprintf('95%%: %5.3f\n\n', get95(sfu_de_apr_int));

sfu_de_opt = calculate_de(sfu_RGBs_opt, sfu_XYZs_tgt, sfu_wp_idx);
fprintf('SFU simulation CIEDE2000 under optimized illuminant using optimized matrix \n');
fprintf('mean: %5.3f\n', mean(sfu_de_opt));
fprintf('max: %5.3f\n', max(sfu_de_opt));
fprintf('95%%: %5.3f\n\n', get95(sfu_de_opt));

sfu_de_opt_int = calculate_de(sfu_RGBs_opt_int, sfu_XYZs_tgt, sfu_wp_idx);
fprintf('SFU simulation CIEDE2000 under optimized illuminant using fitted matrix from integration \n');
fprintf('mean: %5.3f\n', mean(sfu_de_opt_int));
fprintf('max: %5.3f\n', max(sfu_de_opt_int));
fprintf('95%%: %5.3f\n\n', max(sfu_de_opt_int));


function res = get95(arr)
    b = sort(arr);
    idx = round(0.95 * length(arr));
    res = b(idx);
end

function de_vec = calculate_de(XYZs_pred, XYZs_gt, wp_idx)
    ratio = median(XYZs_gt(:, 2) ./ XYZs_pred(:, 2));
    lab_pred = xyz2lab(XYZs_pred * ratio, XYZs_gt(wp_idx, :));
    lab_gt = xyz2lab(XYZs_gt, XYZs_gt(wp_idx, :));
    de_vec = CIE2000deltaE(lab_pred, lab_gt);
end

function cropdata = crop_spectrum(mat, name, lo, hi)
    wave = mat.wavelength;
    ind_lo = find(wave == lo);
    ind_hi = find(wave == hi);
    data = mat.(name);
    cropdata = data(:, ind_lo:ind_hi);
end

function cropdata = crop_lambda(mat, name, lo, hi)
    wave = mat.lambda;
    ind_lo = find(wave == lo);
    ind_hi = find(wave == hi);
    data = mat.(name);
    cropdata = data(:, ind_lo:ind_hi);
end

function mat_interp = interp_spectrum(mat, x, x_interp)
    num = size(mat, 1);
    dim = size(x_interp, 2);
    mat_interp = zeros(num, dim);
    for i = 1:num
        mat_interp(i, :) = interp1(x, mat(i, :), x_interp, 'linear');
    end
end

function loss = spec_match_loss(x, k, bases, C_rgb, illum_cmf, beta, gamma, noise)
    % Dimensions: [lambda, band, illuminant]
    D = size(bases, 1);
    C = size(bases, 2);
    w = reshape(x(1:(C*k)), [1, C, k]);
    M = reshape(permute(reshape(x((C*k+1):end), [3, 3, k]), [1, 3, 2]), [3*k, 3]);
    illuminants = sum(w .* bases, 2); % [D, 1, k]
        
    % uniformly sample linear scale
    [Y, X, Z] = meshgrid(0.1:0.1:1, 0.1:0.1:1, 0.1:0.1:1);
    V = cat(2, X(:), Y(:), Z(:));

    V = repmat(V, 1, k);
    s = repmat(noise.con, 1, k) .* ones(size(V)) / (2^14-1);
    aV = repmat(noise.lin, 1, k) .* ones(size(V)) .* V;
    a2V = repmat((noise.lin.^2), 1, k) .* ones(size(V)) .* V;
    nume = (aV * M) .^ 2;
    denom = (a2V + s) * (M .^ 2);
    % Expected signal-to-noise-ratio
    SNR = 10 * (log10(2^14-1) + log10(nume ./ denom));
    E_SNR = mean(SNR, 'all');
    data = norm(reshape( illuminants .* C_rgb, [D, 3*k] ) * M - illum_cmf, 'fro');
    loss = data + beta * norm(M, 'fro') - gamma * E_SNR;
    % loss = data - beta * mean(illuminants(:)) - gamma * E_SNR;
end

function [data, E_SNR] = spec_match(x, k, bases, C_rgb, illum_cmf, noise)
    % Dimensions: [lambda, band, illuminant]
    D = size(bases, 1);
    C = size(bases, 2);
    w = reshape(x(1:(C*k)), [1, C, k]);
    M = reshape(permute(reshape(x((C*k+1):end), [3, 3, k]), [1, 3, 2]), [3*k, 3]);
    illuminants = sum(w .* bases, 2); % [D, 1, k]
        
    % uniformly sample linear scale
    [Y, X, Z] = meshgrid(0.1:0.1:1, 0.1:0.1:1, 0.1:0.1:1);
    V = cat(2, X(:), Y(:), Z(:));

    V = repmat(V, 1, k);
    s = repmat(noise.con, 1, k) .* ones(size(V)) / (2^14-1);
    aV = repmat(noise.lin, 1, k) .* ones(size(V)) .* V;
    a2V = repmat((noise.lin.^2), 1, k) .* ones(size(V)) .* V;
    nume = (aV * M) .^ 2;
    denom = (a2V + s) * (M .^ 2);
    % Expected signal-to-noise-ratio
    SNR = 10 * (log10(2^14-1) + log10(nume ./ denom));
    E_SNR = mean(SNR, 'all');
    data = norm(reshape( illuminants .* C_rgb, [D, 3*k] ) * M - illum_cmf, 'fro');
end

function noise = noise_param()
    % E(Y(p))   = phi(p) * t * g * kc
    % Var(Y(p)) = E(Y(p)) * g * kc + (sr^2 * g^2 * kc^2 + sa^2 * kc^2)
    % value range from 0 to 2^14-1
    gain = 1;   % ISO = 100
    kc = [0.422, 0.384, 0.389];
    sr = 0.705; % sigma_read
    sa = 3.028; % sigma_adc
    noise.con = (sr^2 * gain^2 * kc.^2) + (sa^2 * kc.^2);
    noise.lin = gain * kc;
    noise.qua = 0;
end
