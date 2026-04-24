function smaps = process_smaps(smaps_raw, emaps, fov_gre, fov, ...
    Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, Nvcoils, SENSEmethod, threshold_mask)
%PROCESS_SMAPS  Mask, crop, and resize raw sensitivity maps to EPI dimensions.
%
%   smaps = process_smaps(smaps_raw, emaps, fov_gre, fov,
%               Nx_gre, Ny_gre, Nz_gre, Nx, Ny, Nz, Nvcoils,
%               SENSEmethod, threshold_mask)
%
%   This function consolidates the three-step sensitivity-map post-processing
%   pipeline that was previously duplicated in main.m and cg_sense.m:
%     1. Build a support mask from the last eigenvalue map.
%     2. Crop the z dimension to match the EPI field of view.
%     3. Interpolate (x,y,z) to match the EPI acquisition grid.
%
%   Inputs
%   ------
%   smaps_raw      Unprocessed sensitivity maps  [Nx_gre × Ny_gre × Nz_gre × Nvcoils]
%   emaps          Eigenvalue maps returned by makeSmaps  [... × Neigs]
%   fov_gre        GRE field of view, [fx fy fz] in metres
%   fov            EPI field of view, [fx fy fz] in metres
%   Nx_gre, Ny_gre, Nz_gre   GRE grid dimensions
%   Nx, Ny, Nz     EPI grid dimensions
%   Nvcoils        Number of virtual coils after compression
%   SENSEmethod    'bart' or 'pisco'  (controls eigenmap sign convention)
%   threshold_mask Scalar threshold applied to the last eigenvalue map.
%                  Voxels below this value are retained; others are zeroed.
%
%   Output
%   ------
%   smaps          Processed sensitivity maps  [Nx × Ny × Nz × Nvcoils]

    %% ── Input validation ─────────────────────────────────────────────────────
    validateattributes(smaps_raw, {'numeric'}, ...
        {'ndims', 4, 'size', [Nx_gre, Ny_gre, Nz_gre, Nvcoils]}, ...
        'process_smaps', 'smaps_raw');
    validateattributes(emaps, {'numeric'}, {'nonempty'}, ...
        'process_smaps', 'emaps');
    % Spatial dimensions of emaps must match the GRE grid — number of
    % trailing dimensions (coils, eigenvalue index, etc.) is method-dependent.
    assert(size(emaps, 1) == Nx_gre && size(emaps, 2) == Ny_gre && size(emaps, 3) == Nz_gre, ...
        'process_smaps: emaps spatial dims [%d %d %d] do not match GRE grid [%d %d %d].', ...
        size(emaps,1), size(emaps,2), size(emaps,3), Nx_gre, Ny_gre, Nz_gre);
    validateattributes(threshold_mask, {'numeric'}, {'scalar', 'real', 'finite'}, ...
        'process_smaps', 'threshold_mask');
    assert(ismember(SENSEmethod, {'bart', 'pisco'}), ...
        'process_smaps: SENSEmethod must be ''bart'' or ''pisco'', got ''%s''.', ...
        SENSEmethod);
    assert(all(fov_gre(3) >= fov(3)), ...
        'process_smaps: EPI z-FoV (%.3f m) exceeds GRE z-FoV (%.3f m).', ...
        fov(3), fov_gre(3));

    %% ── 1. Eigenvalue support mask ───────────────────────────────────────────
    % BART ESPIRiT returns emaps = (1 - true_eig), so invert before thresholding.
    if strcmp(SENSEmethod, 'bart')
        emaps = 1 - emaps;
    end

    % Last eigenvalue map carries the most information about signal support.
    eig_last = emaps(:, :, :, end);
    eig_mask = double(eig_last < threshold_mask);  % 1 inside object, 0 outside

    smaps = smaps_raw .* eig_mask;

    %% ── 2. Crop z to match EPI FoV ───────────────────────────────────────────
    z_frac  = (fov_gre(3) - fov(3)) / fov_gre(3) / 2;
    z_start = round(z_frac * Nz_gre + 1);
    z_end   = round(Nz_gre - z_frac * Nz_gre);

    if z_start < 1 || z_end > Nz_gre || z_start > z_end
        error('process_smaps: computed z crop [%d, %d] is out of range [1, %d].', ...
            z_start, z_end, Nz_gre);
    end

    smaps = smaps(:, :, z_start:z_end, :);

    %% ── 3. Interpolate to EPI grid ───────────────────────────────────────────
    smaps_new = zeros(Nx, Ny, Nz, Nvcoils, 'like', smaps_raw);
    for coil = 1:Nvcoils
        smaps_new(:, :, :, coil) = imresize3(smaps(:, :, :, coil), [Nx, Ny, Nz]);
    end
    smaps = smaps_new;
end