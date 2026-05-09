function [Img_recon] = cg_sense(KData_zf, Sens_maps, num_iter, tol)
    % CG_SENSE  Iterative conjugate-gradient SENSE reconstruction.
    %
    %   Img_recon = cg_sense(KData_zf, Sens_maps, num_iter, tol)
    %
    %   Inputs
    %     KData_zf  - Zero-filled k-space  [Nx, Ny, Nz, ..., NumCoils]
    %     Sens_maps - Coil sensitivity maps [Nx, Ny, Nz, ..., NumCoils]
    %     num_iter  - Maximum CG iterations 
    %     tol       - Relative tolerance for stopping criteria (e.g., 1e-5)

    if nargin < 4 || isempty(tol)
        tol = 1e-9;
    end

    % Identify dimensions
    num_dims = ndims(KData_zf);
    coil_dim = num_dims; 
    spatial_dims = 1:(num_dims - 1); 
    
    % Extract the sampling mask from just the first coil
    idx = repmat({':'}, 1, num_dims);
    idx{coil_dim} = 1;
    Mask = abs(KData_zf(idx{:})) > 0;
    
    % Define the Forward Operator (E) and Adjoint Operator (Eh)
    E  = @(X) Mask .* fftc(Sens_maps .* X, spatial_dims);
    Eh = @(Y) sum(conj(Sens_maps) .* ifftc(Mask .* Y, spatial_dims), coil_dim);

    % --- Initialization for Conjugate Gradient ---
    B = Eh(KData_zf); 
    X = zeros(size(B), 'like', B); 
    
    R = B; 
    P = R; 
    rsold = sum(abs(R(:)).^2);
    rs0 = rsold;
    
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
        
        % Calculate relative residual
        rel_res = sqrt(rsnew / rs0);
        
        msg = sprintf('Iteration %d/%d - Rel Residual: %e', iter, num_iter, rel_res);
        fprintf([repmat('\b', 1, prev_len), '%s'], msg);
        prev_len = length(msg);

        if rel_res < tol
            fprintf('\n  Converged early at iteration %d\n', iter);
            break;
        end

        beta = rsnew / rsold;
        P = R + beta * P;
        rsold = rsnew;
    end
    fprintf('\nReconstruction complete.\n');
    
    Img_recon = X;
end

% --- Generalized N-Dimensional Centered & Unitary FFT Helpers ---
function Res = fftc(X, dims)
    % Performs unitary centered FFT iteratively along specified dimensions
    Res = X;
    for d = dims
        N = size(Res, d);
        Res = fftshift(fft(ifftshift(Res, d), [], d), d) / sqrt(N);
    end
end

function Res = ifftc(X, dims)
    % Performs unitary centered IFFT iteratively along specified dimensions
    Res = X;
    for d = dims
        N = size(Res, d);
        Res = fftshift(ifft(ifftshift(Res, d), [], d), d) * sqrt(N);
    end
end