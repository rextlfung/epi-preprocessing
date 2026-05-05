function [Img_recon] = cg_sense(KData_zf, Sens_maps, num_iter)
    % CG_SENSE  Iterative conjugate-gradient SENSE reconstruction.
    %
    %   Img_recon = cg_sense(KData_zf, Sens_maps, num_iter) reconstructs a
    %   single temporal frame from its zero-filled, multi-coil k-space using
    %   the standard CG-SENSE algorithm (Pruessmann et al., MRM 2001).  The
    %   implementation is generalised to arbitrary spatial dimensionality:
    %   coils are assumed to occupy the last dimension while all preceding
    %   dimensions are treated as spatial.
    %
    %   The cost function minimised is the unregularised least-squares:
    %       min_x  || E x - b ||_2^2
    %   where E = F_mask * S (masked FFT times coil maps) and b = KData_zf.
    %
    %   Inputs
    %     KData_zf  - Zero-filled k-space  [Nx, Ny, Nz, ..., NumCoils]
    %     Sens_maps - Coil sensitivity maps [Nx, Ny, Nz, ..., NumCoils]
    %     num_iter  - Maximum CG iterations (stops early if residual < 1e-10)
    %
    %   Outputs
    %     Img_recon - Reconstructed image   [Nx, Ny, Nz, ...]
    %
    %   References
    %     Pruessmann et al., MRM 2001. https://doi.org/10.1002/mrm.1241
    %     Conjugate gradient method: https://en.wikipedia.org/wiki/Conjugate_gradient_method
    %
    %   See also: run_cg_sense, process_smaps

    % Identify dimensions
    num_dims = ndims(KData_zf);
    coil_dim = num_dims; 
    spatial_dims = 1:(num_dims - 1); 
    
    % Extract the sampling mask from just the first coil
    % We use a dynamic cell array to index (:, :, ..., 1) for any dimension
    idx = repmat({':'}, 1, num_dims);
    idx{coil_dim} = 1;
    Mask = abs(KData_zf(idx{:})) > 0;
    
    % Define the Forward Operator (E): Image -> k-space
    E = @(X) Mask .* fftc(Sens_maps .* X, spatial_dims);

    % Define the Adjoint Operator (Eh): k-space -> Image
    Eh = @(Y) sum(conj(Sens_maps) .* ifftc(Mask .* Y, spatial_dims), coil_dim);

    % --- Initialization for Conjugate Gradient ---
    B = Eh(KData_zf); 
    X = zeros(size(B), 'like', B); 
    
    R = B; 
    P = R; 
    rsold = sum(abs(R(:)).^2);
    
    % --- Conjugate Gradient Loop ---
    fprintf('Starting %dD CG-SENSE reconstruction...\n', length(spatial_dims));

    prev_len = 0;
    for iter = 1:num_iter
        EP = E(P);
        EHEP = Eh(EP);

        alpha = rsold / sum(real(conj(P(:)) .* EHEP(:)));

        X = X + alpha * P;
        R = R - alpha * EHEP;

        rsnew = sum(abs(R(:)).^2);
        msg = sprintf('Iteration %d/%d - Residual: %e', iter, num_iter, rsnew);
        fprintf([repmat('\b', 1, prev_len), '%s'], msg);
        prev_len = length(msg);

        if sqrt(rsnew) < 1e-10
            fprintf('\n  Converged early at iteration %d\n', iter);
            break;
        end

        beta = rsnew / rsold;
        P = R + beta * P;
        rsold = rsnew;
    end
    fprintf('\n');
    
    Img_recon = X;
    fprintf('Reconstruction complete.\n');
end

% --- Generalized N-Dimensional Centered FFT Helpers ---
function Res = fftc(X, dims)
    % Performs centered FFT iteratively along specified dimensions
    Res = X;
    for d = dims
        Res = fftshift(fft(ifftshift(Res, d), [], d), d);
    end
end

function Res = ifftc(X, dims)
    % Performs centered IFFT iteratively along specified dimensions
    Res = X;
    for d = dims
        Res = fftshift(ifft(ifftshift(Res, d), [], d), d);
    end
end