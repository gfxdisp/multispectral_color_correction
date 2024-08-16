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

% define the illuminant
% 390-780 step=1
illum_std = load('illum_D65.mat');
illum_std = crop_spectrum(illum_std, 'power', low, high);

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

% approximate SFU using PMCC(24 chromatic + 1 achromatic)
if ~exist('sfu_coef', 'var')
    sfu_coef = fit_reflectance(sfu, pmcc(1:25, :));
    sfu_recons = sfu_coef * pmcc(1:25, :);
    target = sfu * (illum_std' .* cmf');
    pred = sfu_recons * (illum_std' .* cmf');
    DE = calculate_de(pred, target, 55);
    good_idx = (DE' <= 1) & (min(sfu_recons, [], 2) >= 0);
    sfu_good = sfu(good_idx, :);
    sfu_good_recons = sfu_recons(good_idx, :);
    sfu_good_coef = sfu_coef(good_idx, :);
end

% define camera
cam = 'A7R3';
if strcmp(cam, 'A7R3')
    % define the self-calculated Sony A7R3 SSF 
    % 380-780 step=2
    camera = load('SSF_240506/A7R3_390_780_a0.001_r40_iter1.mat');
    ssf = crop_spectrum(camera, 'rgb', low, high);
    ssf = interp_spectrum(ssf, lambda2, lambda);
elseif strcmp(cam, 'IDS')
    % define the self-calculated IDS SSF 
    % 380-780 step=2
    camera = load('SSF_240506/IDS_390_780_a0.001_r40_iter5.mat');
    ssf = crop_spectrum(camera, 'rgb', low, high);
    ssf = interp_spectrum(ssf, lambda2, lambda);
else
    disp('wrong camera')
    pause
end

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
specbos_xyz_lab = xyz2lab(specbos_xyz_int, specbos_xyz_int(wp_idx, :));

% read image patches under the approximated illuminant
imgs = dir(fullfile(path, cam, 'd65', '*.exr'));
images = cell(1, length(imgs));
im_apr = zeros(length(imgs), 3);
radius = 40; % radius of a patch to take image mean value
N = (2*radius+1) ^ 2;
im_val_apr = zeros(length(imgs), 3, N);
for i = 1:length(imgs)
    images{i} = pfs_read_image(fullfile(imgs(i).folder, imgs(i).name));
    patch = images{i}(400-radius:400+radius,400-radius:400+radius, :);
    im_apr(i, :) = mean(patch, [1, 2]);
    im_val_apr(i, :, :) = reshape(patch, N, 3)';
end
m = mean(im_val_apr, 3);
v = var(im_val_apr, 0, 3);
im_apr_snr = 10 * log10(m.^2 ./ v);
disp('SNR of the Approx illuminant');
disp(mean(im_apr_snr, 'all'));
gt_apr = pmcc * (illum_apr' .* ssf'); % integration result should be similar to im_apr after rescaling
rat_apr = (im_apr ./ im_apr(wp_idx, :)) ./ (gt_apr ./ gt_apr(wp_idx, :));
dif_apr = mean(abs(rat_apr - 1), 'all');

% Luther condition CCM to match camera SSF to CMF: ssf' * M = cmf'
M_lu = ssf' \ cmf';

% LCC from integrated camera RGBs to CMF XYZs under the approximated illuminant
wp_idx = 25; % white point index
RGBs = pmcc * (illum_apr' .* ssf');
XYZs = pmcc * (illum_apr' .* cmf');
M_apr_int = RGBs \ XYZs;

% LCC from real image RGBs to specbos XYZs under the approximated illuminant
M_apr_real = im_apr \ specbos_xyz_int;

de_apr_img_lu = calculate_de(im_apr * M_lu, specbos_xyz_int, wp_idx);
fprintf('PMCC real image CIEDE2000 under Approx illuminant using Luther matrix \n');
fprintf('mean: %5.3f\n', mean(de_apr_img_lu));
fprintf('max: %5.3f\n\n', max(de_apr_img_lu));

de_apr_img_int = calculate_de(im_apr * M_apr_int, specbos_xyz_int, wp_idx);
fprintf('PMCC real image CIEDE2000 under Approx illuminant using fitted matrix from integration \n');
fprintf('mean: %5.3f\n', mean(de_apr_img_int));
fprintf('max: %5.3f\n\n', max(de_apr_img_int));

de_apr_img_real = calculate_de(im_apr * M_apr_real, specbos_xyz_int, wp_idx);
fprintf('PMCC real image CIEDE2000 under Approx illuminant using fitted matrix from real image \n');
fprintf('mean: %5.3f\n', mean(de_apr_img_real));
fprintf('max: %5.3f\n\n', max(de_apr_img_real));

% read and color correct image patches under optimized lights
k = 1;
if strcmp(cam, 'A7R3')
    load(fullfile(path, [cam, '_k', num2str(k) ,'_b1.00_g0.10_390_780.mat']));
    index = [1, 1, 6];
elseif strcmp(cam, 'IDS')
    load(fullfile(path, [cam, '_k', num2str(k) ,'_b0.20_g0.10_390_780.mat']));
    index = [1, 14, 17];
else
    disp('wrong camera')
    pause
end
im_val_opt = zeros(length(imgs), k*3, N);
idx = index(k);
M_opt = ll_M{idx};
w_opt = ll_w{idx};
im_opt = zeros(length(imgs), 3, k);
rat_opt = zeros(length(imgs), 3, k);
dif_opt = zeros(1, k);
im_opt_xyz = zeros(length(imgs), 3);
for kk = 1:k
    imgs = dir(fullfile(path, cam, sprintf('k%d%d', k, kk), '*.exr'));
    images = cell(1, length(imgs));
    for i = 1:length(imgs)
        images{i} = pfs_read_image(fullfile(imgs(i).folder, imgs(i).name));
        patch = images{i}(400-radius:400+radius,400-radius:400+radius, :);
        im_opt(i, :, kk) = mean(patch, [1, 2]);
        im_val_opt(i, kk*3-2:kk*3, :) = reshape(patch, N, 3)';
    end
    gt_kk = pmcc * ((w_opt(:, :, kk) * bases)' .* ssf'); % integration result should be similar to im_all{kk}
    rat_opt(:, :, kk) = squeeze(im_opt(:, :, kk) ./ im_opt(wp_idx, :, kk)) ./ (gt_kk ./ gt_kk(wp_idx, :));
    dif_opt(kk) = mean(abs(rat_opt(:, :, kk) - 1), 'all');
    im_opt_xyz = im_opt_xyz + im_opt(:, :, kk) * M_opt(kk*3-2:kk*3, :);
end

im_val_opt = reshape(im_val_opt, length(imgs), 3*k, 1, N);
M_opt_reshape = reshape(M_opt, 1, size(M_opt, 1), size(M_opt, 2), 1);
im_val_opt_xyz = im_val_opt .* M_opt_reshape;
im_val_opt_xyz = squeeze(sum(im_val_opt_xyz, 2));
m = mean(im_val_opt_xyz, 3);
v = var(im_val_opt_xyz, 0, 3);
im_opt_snr = 10 * log10(m.^2 ./ v);
disp('SNR of the optimized illuminant');
disp(mean(im_opt_snr, 'all'));

% LCC from real image RGBs to specbos XYZs under optimized illuminants
im_opt = reshape(im_opt, length(imgs), 3*k);
M_opt_real = im_opt \ specbos_xyz_int;

de_opt_img = calculate_de(im_opt_xyz, specbos_xyz_int, wp_idx);
fprintf('PMCC real image CIEDE2000 under optimized illuminant using optimized matrix \n');
fprintf('mean: %5.3f\n', mean(de_opt_img));
fprintf('max: %5.3f\n\n', max(de_opt_img));

de_opt_img_real = calculate_de(im_opt * M_opt_real, specbos_xyz_int, wp_idx);
fprintf('PMCC real image CIEDE2000 under optimized illuminant using fitted matrix from real image \n');
fprintf('mean: %5.3f\n', mean(de_opt_img_real));
fprintf('max: %5.3f\n\n', max(de_opt_img_real));


% synthetic SFU images basaed on PMCC image value reconstruction
tmp = sum(sfu_good, 2);
[~, sfu_wp_idx] = max(tmp); % find the white point
sfu_specbos_xyz = sfu_good_coef * specbos_xyz_int(1:25, :);
sfu_im_apr = sfu_good_coef * im_apr(1:25, :);
sfu_im_opt = sfu_good_coef * im_opt(1:25, :);

sfu_de_apr_img_lu = calculate_de(sfu_im_apr * M_lu, sfu_specbos_xyz, sfu_wp_idx);
fprintf('SFU synthetic CIEDE2000 under Approx illuminant using Luther matrix \n');
fprintf('mean: %5.3f\n', mean(sfu_de_apr_img_lu));
fprintf('max: %5.3f\n', max(sfu_de_apr_img_lu));
fprintf('95%%: %5.3f\n\n', get95(sfu_de_apr_img_lu));

sfu_de_apr_img_real = calculate_de(sfu_im_apr * M_apr_real, sfu_specbos_xyz, sfu_wp_idx);
fprintf('SFU synthetic CIEDE2000 under Approx illuminant using fitted matrix from real image \n');
fprintf('mean: %5.3f\n', mean(sfu_de_apr_img_real));
fprintf('max: %5.3f\n', max(sfu_de_apr_img_real));
fprintf('95%%: %5.3f\n\n', get95(sfu_de_apr_img_real));

sfu_de_opt_img = calculate_de(sfu_im_opt * M_opt, sfu_specbos_xyz, sfu_wp_idx);
fprintf('SFU synthetic CIEDE2000 under optimized illuminant using optimized matrix \n');
fprintf('mean: %5.3f\n', mean(sfu_de_opt_img));
fprintf('max: %5.3f\n', max(sfu_de_opt_img));
fprintf('95%%: %5.3f\n\n', get95(sfu_de_opt_img));

sfu_de_opt_img_real = calculate_de(sfu_im_opt * M_opt_real, sfu_specbos_xyz, sfu_wp_idx);
fprintf('SFU synthetic CIEDE2000 under optimized illuminant using fitted matrix from real image \n');
fprintf('mean: %5.3f\n', mean(sfu_de_opt_img_real));
fprintf('max: %5.3f\n', max(sfu_de_opt_img_real));
fprintf('95%%: %5.3f\n\n', get95(sfu_de_opt_img_real));



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

function weight_mat = fit_reflectance(target, source)
    % for every row in the target, find the closest fit from 4 rows in source
    tic
    N = size(target, 1);
    k = size(source, 1);
    weight_mat = zeros(N, k);
    for i = 1:N
        best = 1e5;
        for u = 1:k-3
            for v = u+1:k-2
                for p = v+1:k-1
                    for q = p+1:k
                        A = [source(u,:)' source(v,:)' source(p,:)' source(q,:)'];
                        b = target(i,:)';
                        x = A \ b;
                        error = norm(A * x - b);
                        if error < best
                            best = error;
                            best_idx = [u,v,p,q];
                            best_x = x;
                        end
                    end
                end
            end
        end
        for j = 1:4
            weight_mat(i, best_idx(j)) = best_x(j);
        end
    end
    toc
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