function [u,us]=MCfieldOnWave(sigt,albedo,box_min,box_max,maxItr,pathsNum,lambda,smpFlg,sct_type,ampfunc,ampfunc0, ...
    l,v,dirl,dirv,refocus)
%% Input
isRefocus = (nargin > 15);

if(nargin < 14)
    dirl = l;
    dirv = v;
end

dim = size(box_min,1);

%% Calculate dims
MAX_DIM = 32;

% af_ang_vl
af_ang_vl_maxDim = max( ...
    [ndims(l.x), ndims(l.y), ndims(l.z), ...
     ndims(v.x), ndims(v.y), ndims(v.z)]);

af_ang_vl_lx_size = size(l.x); af_ang_vl_lx_size(end+1:af_ang_vl_maxDim) = 1;
af_ang_vl_ly_size = size(l.y); af_ang_vl_ly_size(end+1:af_ang_vl_maxDim) = 1;
af_ang_vl_lz_size = size(l.z); af_ang_vl_lz_size(end+1:af_ang_vl_maxDim) = 1;

af_ang_vl_vx_size = size(v.x); af_ang_vl_vx_size(end+1:af_ang_vl_maxDim) = 1;
af_ang_vl_vy_size = size(v.y); af_ang_vl_vy_size(end+1:af_ang_vl_maxDim) = 1;
af_ang_vl_vz_size = size(v.z); af_ang_vl_vz_size(end+1:af_ang_vl_maxDim) = 1;

af_ang_vl_dim = max([ ...
    af_ang_vl_lx_size; af_ang_vl_ly_size; af_ang_vl_lz_size; ...
    af_ang_vl_vx_size; af_ang_vl_vy_size; af_ang_vl_vz_size]);

af_ang_vl_dimProd = zeros(6,MAX_DIM);
af_ang_vl_dimProd(1 ,1:numel(af_ang_vl_vx_size)) = cumprod(af_ang_vl_vx_size);
af_ang_vl_dimProd(2 ,1:numel(af_ang_vl_vy_size)) = cumprod(af_ang_vl_vy_size);
af_ang_vl_dimProd(3 ,1:numel(af_ang_vl_vz_size)) = cumprod(af_ang_vl_vz_size);
af_ang_vl_dimProd(4 ,1:numel(af_ang_vl_lx_size)) = cumprod(af_ang_vl_lx_size);
af_ang_vl_dimProd(5 ,1:numel(af_ang_vl_ly_size)) = cumprod(af_ang_vl_ly_size);
af_ang_vl_dimProd(6 ,1:numel(af_ang_vl_lz_size)) = cumprod(af_ang_vl_lz_size);

zerosIdx = ([ones(6,1),diff(af_ang_vl_dimProd,[],2)] ~= 0);

af_ang_vl_dimProd = circshift(af_ang_vl_dimProd,1,2);
af_ang_vl_dimProd = zerosIdx .* af_ang_vl_dimProd;

% el + el_ms
el_maxDim = max( ...
    [ndims(l.x)   , ndims(l.y)   , ndims(l.z)   , ...
     ndims(dirl.x), ndims(dirl.y), ndims(dirl.z), ...
     ndims(lambda)]);
 el_maxDim = max(3,el_maxDim);
 
el_lx_size = size(l.x); el_lx_size(end+1:el_maxDim) = 1;
el_ly_size = size(l.y); el_ly_size(end+1:el_maxDim) = 1;
el_lz_size = size(l.z); el_lz_size(end+1:el_maxDim) = 1;

el_dirlx_size = size(dirl.x); el_dirlx_size(end+1:el_maxDim) = 1;
el_dirly_size = size(dirl.y); el_dirly_size(end+1:el_maxDim) = 1;
el_dirlz_size = size(dirl.z); el_dirlz_size(end+1:el_maxDim) = 1;

el_wavelength_size = size(lambda); el_wavelength_size(end+1:el_maxDim) = 1;

el_dim = max([ ...
    el_lx_size   ; el_ly_size   ; el_lz_size   ; ...
    el_dirlx_size; el_dirly_size; el_dirlz_size; ...
    el_wavelength_size]);

el_dim(3) = pathsNum;

el_dimProd = zeros(7,MAX_DIM);
el_dimProd(1 ,1:numel(el_lx_size))         = cumprod(el_lx_size)   ;
el_dimProd(2 ,1:numel(el_ly_size))         = cumprod(el_ly_size)   ;
el_dimProd(3 ,1:numel(el_lz_size))         = cumprod(el_lz_size)   ;
el_dimProd(4 ,1:numel(el_dirlx_size))      = cumprod(el_dirlx_size);
el_dimProd(5 ,1:numel(el_dirly_size))      = cumprod(el_dirly_size);
el_dimProd(6 ,1:numel(el_dirlz_size))      = cumprod(el_dirlz_size);
el_dimProd(7 ,1:numel(el_wavelength_size)) = cumprod(el_wavelength_size);

zerosIdx = ([ones(7,1),diff(el_dimProd,[],2)] ~= 0);

el_dimProd = circshift(el_dimProd,1,2);
el_dimProd = zerosIdx .* el_dimProd;

% ff_scattering
ff_scattering_maxDim = max( ...
    [af_ang_vl_maxDim, el_maxDim, ...
     ndims(dirl.x), ndims(dirl.y), ndims(dirl.z), ...
     ndims(dirv.x), ndims(dirv.y), ndims(dirv.z)]);

ff_scattering_vx_size = size(v.x); ff_scattering_vx_size(end+1:ff_scattering_maxDim) = 1;
ff_scattering_vy_size = size(v.y); ff_scattering_vy_size(end+1:ff_scattering_maxDim) = 1;
ff_scattering_vz_size = size(v.z); ff_scattering_vz_size(end+1:ff_scattering_maxDim) = 1;

ff_scattering_dirvx_size = size(dirv.x); ff_scattering_dirvx_size(end+1:ff_scattering_maxDim) = 1;
ff_scattering_dirvy_size = size(dirv.y); ff_scattering_dirvy_size(end+1:ff_scattering_maxDim) = 1;
ff_scattering_dirvz_size = size(dirv.z); ff_scattering_dirvz_size(end+1:ff_scattering_maxDim) = 1;

ff_scattering_af_ang_vl_dim = af_ang_vl_dim; ff_scattering_af_ang_vl_dim(end+1:ff_scattering_maxDim) = 1;
ff_scattering_el_dim = el_dim; ff_scattering_el_dim(end+1:ff_scattering_maxDim) = 1;

ff_wavelength_size = size(lambda); ff_wavelength_size(end+1:ff_scattering_maxDim) = 1;

ff_scattering_dim = max([ ff_scattering_af_ang_vl_dim ; ff_scattering_el_dim ; ff_wavelength_size ; ...
    ff_scattering_vx_size   ; ff_scattering_vy_size   ; ff_scattering_vz_size   ; ...
    ff_scattering_dirvx_size; ff_scattering_dirvy_size; ff_scattering_dirvz_size]);

ff_scattering_dimProd = zeros(9,MAX_DIM);
ff_scattering_dimProd(1 ,1:numel(ff_scattering_vx_size))       = cumprod(ff_scattering_vx_size)      ;
ff_scattering_dimProd(2 ,1:numel(ff_scattering_vy_size))       = cumprod(ff_scattering_vy_size)      ;
ff_scattering_dimProd(3 ,1:numel(ff_scattering_vz_size))       = cumprod(ff_scattering_vz_size)      ;
ff_scattering_dimProd(4 ,1:numel(ff_scattering_dirvx_size))    = cumprod(ff_scattering_dirvx_size)   ;
ff_scattering_dimProd(5 ,1:numel(ff_scattering_dirvy_size))    = cumprod(ff_scattering_dirvy_size)   ;
ff_scattering_dimProd(6 ,1:numel(ff_scattering_dirvz_size))    = cumprod(ff_scattering_dirvz_size)   ;
ff_scattering_dimProd(7 ,1:numel(ff_scattering_af_ang_vl_dim)) = cumprod(ff_scattering_af_ang_vl_dim);
ff_scattering_dimProd(8 ,1:numel(ff_scattering_el_dim))        = cumprod(ff_scattering_el_dim)       ;
ff_scattering_dimProd(9 ,1:numel(ff_wavelength_size))          = cumprod(ff_wavelength_size)         ;

zerosIdx = ([ones(9,1),diff(ff_scattering_dimProd,[],2)] ~= 0);
ff_scattering_dimProd = circshift(ff_scattering_dimProd,1,2);
ff_scattering_dimProd = zerosIdx .* ff_scattering_dimProd;
ff_scattering_dimProd = circshift(ff_scattering_dimProd,-2,2);

ff_scattering_threads_dim = ff_scattering_dim(3:end);
ff_u_size = ff_scattering_dim(4:end);

if(numel(ff_u_size) == 1)
    ff_u_size = [ff_u_size,1];
end

%% Check legal directions

% directions cannot be complex
l.x = real(l.x); l.y = real(l.y); l.z = real(l.z);
v.x = real(v.x); v.y = real(v.y); v.z = real(v.z);
dirl.x = real(dirl.x); dirl.y = real(dirl.y); dirl.z = real(dirl.z);
dirv.x = real(dirv.x); dirv.y = real(dirv.y); dirv.z = real(dirv.z);

% directions norm is 1
ff_legal_idx = 1;
ff_legal_idx = ff_legal_idx .* ((l.x.^2 + l.y.^2 + l.z.^2) > 0.9999999);
ff_legal_idx = ff_legal_idx .* ((l.x.^2 + l.y.^2 + l.z.^2) < 1.0000001);
ff_legal_idx = ff_legal_idx .* ((v.x.^2 + v.y.^2 + v.z.^2) > 0.9999999);
ff_legal_idx = ff_legal_idx .* ((v.x.^2 + v.y.^2 + v.z.^2) < 1.0000001);
ff_legal_idx = ff_legal_idx .* ((dirl.x.^2 + dirl.y.^2 + dirl.z.^2) > 0.9999999);
ff_legal_idx = ff_legal_idx .* ((dirl.x.^2 + dirl.y.^2 + dirl.z.^2) < 1.0000001);
ff_legal_idx = ff_legal_idx .* ((dirv.x.^2 + dirv.y.^2 + dirv.z.^2) > 0.9999999);
ff_legal_idx = ff_legal_idx .* ((dirv.x.^2 + dirv.y.^2 + dirv.z.^2) < 1.0000001);

% in case of correlation computing
ff_legal_idx = all(ff_legal_idx,2)*1;

permuteVec = 1:1:ndims(ff_legal_idx);
permuteVec(1:end-3) = permuteVec(4:end);
permuteVec(end-2:end) = [1,2,3];

ff_legal_idx = permute(ff_legal_idx,permuteVec);

ff_legal_idx = gpuArray(ff_legal_idx);
if(all(ff_legal_idx(:)))
    ff_legal_idx = 1;
end

if(~any(ff_legal_idx(:)))
    error('no legal direction')
end

%% Refocus
if(isRefocus)
    % refocus directions must be legal
    if(     ~isreal(refocus.focalDirectionsL.eval.x) || ...
            ~isreal(refocus.focalDirectionsL.eval.y) || ...
            ~isreal(refocus.focalDirectionsL.eval.z))
        error('refocused directions cannot be ilegal');
    end
    
    dirAbs = refocus.focalDirectionsL.eval.x.^2 + ...
             refocus.focalDirectionsL.eval.y.^2 + ...
             refocus.focalDirectionsL.eval.z.^2;
    dirAbs = dirAbs(:);
    
    if( any(dirAbs < 0.9999999 | dirAbs > 1.0000001))
        error('refocused directions cannot be ilegal');
    end
    
    kappaL = 1/refocus.mask_varL^2;
    kappaV = 1/refocus.mask_varL^2;
        
    ff_maxDim = refocus.ff_dims;
    l_x_size = size(l.x); l_x_size = l_x_size(4:end); l_x_size(end+1:ff_maxDim) = 1;
    l_y_size = size(l.y); l_y_size = l_y_size(4:end); l_y_size(end+1:ff_maxDim) = 1;
    l_z_size = size(l.z); l_z_size = l_z_size(4:end); l_z_size(end+1:ff_maxDim) = 1;
    ff_l_dims = max([l_x_size;l_y_size;l_z_size]);

    v_x_size = size(v.x); v_x_size = v_x_size(4:end); v_x_size(end+1:ff_maxDim) = 1;
    v_y_size = size(v.y); v_y_size = v_y_size(4:end); v_y_size(end+1:ff_maxDim) = 1;
    v_z_size = size(v.z); v_z_size = v_z_size(4:end); v_z_size(end+1:ff_maxDim) = 1;
    ff_v_dims = max([v_x_size;v_y_size;v_z_size]);
    
    if(any((ff_l_dims - ff_v_dims) == 0 & ff_l_dims ~= 1))
        error('in refocus, v and l dims must be separable')
    end
    
    w_l = (kappaL/(2*pi)) .* exp( -kappaL + ...
        kappaL * (l.x .* refocus.focalDirectionsL.eval.x + l.y .* refocus.focalDirectionsL.eval.y + l.z .* refocus.focalDirectionsL.eval.z) + ...
         + 1i * (2*pi ./ lambda) .* (l.x .* refocus.focalPointsL.eval.x + l.y .* refocus.focalPointsL.eval.y + l.z .* refocus.focalPointsL.eval.z)) .* ...
        sqrt(1 - (l.x.^2) ./ (l.x.^2 + l.y.^2 - 1) - (l.y.^2) ./ (l.x.^2 + l.y.^2 - 1)) .* (l.x.^2 + l.y.^2 < 1);
    
    w_v = (kappaV/(2*pi)) .* exp( -kappaV + ...
        kappaV * (v.x .* refocus.focalDirectionsV.eval.x + v.y .* refocus.focalDirectionsV.eval.y + v.z .* refocus.focalDirectionsV.eval.z) + ...
         + 1i * (2*pi ./ lambda) .* (v.x .* refocus.focalPointsV.eval.x + v.y .* refocus.focalPointsV.eval.y + v.z .* refocus.focalPointsV.eval.z)) .* ...
        sqrt(1 - (v.x.^2) ./ (v.x.^2 + v.y.^2 - 1) - (v.y.^2) ./ (v.x.^2 + v.y.^2 - 1)) .* (v.x.^2 + v.y.^2 < 1);
    
    w_l(isnan(w_l)) = 0;
    w_v(isnan(w_v)) = 0;
    
    w_l = conj(w_l);
    
    % v dims | l dims | common dims reshape
    ff_v_dims_num = find(ff_v_dims ~= 1);
    ff_l_dims_num = find(ff_l_dims ~= 1);
    permute_ff_dims = [ff_v_dims_num, ff_l_dims_num];
    
    maxDim = refocus.nf_dims + refocus.common_dims + refocus.ff_dims;
    permute_ff_dims = [permute_ff_dims,(numel(permute_ff_dims)+1):1:maxDim];
        
    nf_maxDim = refocus.nf_dims + refocus.common_dims;
        
    l_x_size = size(refocus.focalPointsL.eval.x); l_x_size = l_x_size((4+ff_maxDim):end); l_x_size(end+1:nf_maxDim) = 1;
    l_y_size = size(refocus.focalPointsL.eval.y); l_y_size = l_y_size((4+ff_maxDim):end); l_y_size(end+1:nf_maxDim) = 1;
    l_z_size = size(refocus.focalPointsL.eval.z); l_z_size = l_z_size((4+ff_maxDim):end); l_z_size(end+1:nf_maxDim) = 1;
    
    l_dirx_size = size(refocus.focalDirectionsL.eval.x); l_dirx_size = l_dirx_size((4+ff_maxDim):end); l_dirx_size(end+1:nf_maxDim) = 1;
    l_diry_size = size(refocus.focalDirectionsL.eval.y); l_diry_size = l_diry_size((4+ff_maxDim):end); l_diry_size(end+1:nf_maxDim) = 1;
    l_dirz_size = size(refocus.focalDirectionsL.eval.z); l_dirz_size = l_dirz_size((4+ff_maxDim):end); l_dirz_size(end+1:nf_maxDim) = 1;
    nf_l_dims = max([l_x_size;l_y_size;l_z_size;l_dirx_size;l_diry_size;l_dirz_size]);

    v_x_size = size(refocus.focalPointsV.eval.x); v_x_size = v_x_size((4+ff_maxDim):end); v_x_size(end+1:nf_maxDim) = 1;
    v_y_size = size(refocus.focalPointsV.eval.y); v_y_size = v_y_size((4+ff_maxDim):end); v_y_size(end+1:nf_maxDim) = 1;
    v_z_size = size(refocus.focalPointsV.eval.z); v_z_size = v_z_size((4+ff_maxDim):end); v_z_size(end+1:nf_maxDim) = 1;
    
    v_dirx_size = size(refocus.focalDirectionsV.eval.x); v_dirx_size = v_dirx_size((4+ff_maxDim):end); v_dirx_size(end+1:nf_maxDim) = 1;
    v_diry_size = size(refocus.focalDirectionsV.eval.y); v_diry_size = v_diry_size((4+ff_maxDim):end); v_diry_size(end+1:nf_maxDim) = 1;
    v_dirz_size = size(refocus.focalDirectionsV.eval.z); v_dirz_size = v_dirz_size((4+ff_maxDim):end); v_dirz_size(end+1:nf_maxDim) = 1;
    nf_v_dims = max([v_x_size;v_y_size;v_z_size;v_dirx_size;v_diry_size;v_dirz_size]);
    
    nf_dims = max([nf_l_dims;nf_v_dims]);
    
    wavelength_size = size(lambda); wavelength_size = wavelength_size((4+ff_maxDim):end); wavelength_size(end+1:nf_maxDim) = 1;
    common_dims = wavelength_size;
    
    cnf_dims = max([nf_dims; common_dims]);
    cnf_v_dims = max([nf_v_dims; common_dims]);
    cnf_l_dims = max([nf_l_dims; common_dims]);
    
    reshape_ff_dims = [prod(ff_v_dims), prod(ff_l_dims),prod(common_dims)];
    
    % permute wl
    w_l = reshape(w_l,size(w_l,2),prod(ff_l_dims),prod(nf_l_dims),prod(common_dims));
    w_l = permute(w_l,[2,1,3,4]);
    w_l = reshape(w_l,size(w_l,1),size(w_l,2) * size(w_l,3),size(w_l,4));
    
    % permute wv
    w_v = reshape(w_v,[size(w_v,2),prod(ff_v_dims),cnf_v_dims]);
    permuteVec = 1:1:ndims(w_v);
    permuteVec(1) = 2; permuteVec(2) = 1;
    w_v = permute(w_v,permuteVec);
    
    d_ff = (l.x(2) - l.x(1)) * (l.y(2) - l.y(1)) * (v.x(2) - v.x(1)) * (v.x(2) - v.x(1));

    refocus_maxDim = 2 + numel(cnf_dims);
        
    refocus_dimProd = zeros(2,MAX_DIM);
    refocus_dimProd(1 ,1:refocus_maxDim) = cumprod([prod(ff_v_dims),size(w_v,2),cnf_v_dims]);
    refocus_dimProd(2 ,1:refocus_maxDim) = cumprod([prod(ff_v_dims),size(w_v,2),cnf_l_dims]);

    zerosIdx = ([ones(2,1),diff(refocus_dimProd,[],2)] ~= 0);
    refocus_dimProd = circshift(refocus_dimProd,1,2);
    refocus_dimProd = zerosIdx .* refocus_dimProd;
    refocus_dimProd = circshift(refocus_dimProd,-2,2);
    
    w_l = gpuArray(complex(w_l));
    w_v = gpuArray(complex(w_v));
    
    u_wl = gpuArray(complex(zeros(prod(ff_v_dims),size(w_l,2),size(w_l,3))));

end

%% Gpu function
% af_ang_vl
gpuFunc.af_ang_vl = parallel.gpu.CUDAKernel('ff_core.ptx','ff_core.cu','af_ang_vl');
gpuFunc.af_ang_vl.GridSize = [ceil(prod(af_ang_vl_dim)/gpuFunc.af_ang_vl.MaxThreadsPerBlock) 1 1];
gpuFunc.af_ang_vl.ThreadBlockSize = [gpuFunc.af_ang_vl.MaxThreadsPerBlock 1 1];

setConstantMemory(gpuFunc.af_ang_vl,'box_min',box_min);
setConstantMemory(gpuFunc.af_ang_vl,'box_max',box_max);
setConstantMemory(gpuFunc.af_ang_vl,'sigt',sigt/2);
setConstantMemory(gpuFunc.af_ang_vl,'dims',uint32(dim));
setConstantMemory(gpuFunc.af_ang_vl,'scatteringType',uint32(sct_type));

if(sct_type == 2)
    setConstantMemory(gpuFunc.af_ang_vl,'nAmpfunc',uint32(ampfunc.evalAmp));
end

if(sct_type == 3)
    setConstantMemory(gpuFunc.af_ang_vl,'fw',ampfunc.forwardWeight);
    setConstantMemory(gpuFunc.af_ang_vl,'g_up',(1 - ampfunc.g * ampfunc.g)/(4*pi));
    setConstantMemory(gpuFunc.af_ang_vl,'g_down1',1 + ampfunc.g * ampfunc.g);
    setConstantMemory(gpuFunc.af_ang_vl,'g_down2',-2 * ampfunc.g);
end

setConstantMemory(gpuFunc.af_ang_vl,'af_ang_vl_dimProd',uint32(cumprod(af_ang_vl_dim)));
setConstantMemory(gpuFunc.af_ang_vl,'af_ang_vl_maxDim',uint32(af_ang_vl_maxDim));
setConstantMemory(gpuFunc.af_ang_vl,'af_ang_vl_inDimProd',uint32(af_ang_vl_dimProd.'));

% el
gpuFunc.el = parallel.gpu.CUDAKernel('ff_core.ptx','ff_core.cu','ff_calcEl');
gpuFunc.el.GridSize = [ceil(prod(el_dim)/gpuFunc.el.MaxThreadsPerBlock) 1 1];
gpuFunc.el.ThreadBlockSize = [gpuFunc.el.MaxThreadsPerBlock 1 1];

setConstantMemory(gpuFunc.el,'el_dimProd',uint32(cumprod(el_dim)));
setConstantMemory(gpuFunc.el,'el_maxDim',uint32(el_maxDim));
setConstantMemory(gpuFunc.el,'el_inDimProd',uint32(el_dimProd.'));

% ff scattering
gpuFunc.ff_scattering = parallel.gpu.CUDAKernel('ff_core.ptx','ff_core.cu','ff_single');
gpuFunc.ff_scattering.GridSize = [ceil(prod(ff_scattering_threads_dim)/gpuFunc.ff_scattering.MaxThreadsPerBlock) 1 1];
gpuFunc.ff_scattering.ThreadBlockSize = [gpuFunc.ff_scattering.MaxThreadsPerBlock 1 1];

setConstantMemory(gpuFunc.ff_scattering,'ff_scattering_dimProd',uint32(cumprod(ff_scattering_threads_dim)));
setConstantMemory(gpuFunc.ff_scattering,'ff_scattering_maxDim',uint32(ff_scattering_maxDim - 2));
setConstantMemory(gpuFunc.ff_scattering,'ff_scattering_inDimProd',uint32(ff_scattering_dimProd.'));
setConstantMemory(gpuFunc.ff_scattering,'ff_corrParamNum',uint32(size(l.x,2)));

% ff multiple scattering
gpuFunc.ff_mscattering = parallel.gpu.CUDAKernel('ff_core.ptx','ff_core.cu','ff_multiple');
gpuFunc.ff_mscattering.GridSize = [ceil(prod(ff_scattering_threads_dim)/gpuFunc.ff_mscattering.MaxThreadsPerBlock) 1 1];
gpuFunc.ff_mscattering.ThreadBlockSize = [gpuFunc.ff_mscattering.MaxThreadsPerBlock 1 1];

% refocus correlation
if(isRefocus)
    if(refocus.correlationActive)
        gpuFunc.refocus = parallel.gpu.CUDAKernel('ff_core.ptx','ff_core.cu','refocus_correlation');
    else
        gpuFunc.refocus = parallel.gpu.CUDAKernel('ff_core.ptx','ff_core.cu','refocus_field');
    end
    
    gpuFunc.refocus.GridSize = [ceil(prod(cnf_dims)/gpuFunc.refocus.MaxThreadsPerBlock) 1 1];
    gpuFunc.refocus.ThreadBlockSize = [gpuFunc.refocus.MaxThreadsPerBlock 1 1];

    setConstantMemory(gpuFunc.refocus,'refocus_dimProd',uint32(cumprod(cnf_dims)));
    setConstantMemory(gpuFunc.refocus,'refocus_maxDim',uint32(refocus_maxDim - 2));
    setConstantMemory(gpuFunc.refocus,'refocus_inDimProd',uint32(refocus_dimProd.'));

    setConstantMemory(gpuFunc.refocus,'refocus_ffElements',uint32(prod(ff_v_dims)));
end

%% Prepare for algorithm

% Box size
box_w = box_max-box_min;

% initiate parameters
if(sct_type == 2)
    tabulatedAmplitudeFunction = gpuArray(ampfunc.evalAmp);
else
    tabulatedAmplitudeFunction = gpuArray(0);
end

af_ang_vl = zeros(af_ang_vl_dim, 'gpuArray');
e_l0 = gpuArray(complex(zeros(el_dim)));
e_l0_ms = gpuArray(complex(zeros(el_dim)));

if(isRefocus && refocus.correlationActive)
    u0_ff = gpuArray(complex(zeros(ff_u_size)));
    u_ff = gpuArray(complex(zeros(ff_u_size)));
else
    u_ff = gpuArray(complex(zeros(ff_u_size)));
    us_ff = gpuArray(complex(zeros(ff_u_size)));
end

if(isRefocus)
    u_nf = gpuArray(complex(zeros(cnf_dims)));
    us_nf = gpuArray(complex(zeros(cnf_dims)));
end

% Pre-calculate single scattering rotation amplitude, only possible when
% both light and view are far field (otherwise it also dependent on the
% first scatter position)

af_ang_vl = feval(gpuFunc.af_ang_vl, af_ang_vl, tabulatedAmplitudeFunction, ...
    v.x, v.y, v.z, l.x, l.y, l.z);

% threshold to begin kill particles with low weight
killThr=0.2;
V=prod(box_w);

% direction 0
d0 = zeros(dim,1); d0(end) = 1;

%% Begin the main loop
for itr=1:maxItr
%     if(mod(itr,1e3) == 0)
%         disp([num2str(round((itr/maxItr)*100)),'%'])
%     end
    
    % Sample the first scatter
    % x: first scattering point
    % px: probability by which first point was sampled. Needed so that
    % importance sampling integration is weighted properly
    switch mod(smpFlg,10)
        case 1
            % uniform distribution
            x=rand(dim,1,pathsNum).*(box_w)+box_min;
            px=ones(1,1,pathsNum);
        case 2
            % exponential distribution
            [x,px] = expSmpX(box_min,box_max,d0,sigt);
%             [x,px] = expSmpX([-10;-35],[10;35],d0,sigt(2));
        case 3
            error('bad!!')
        case 6
            x = ampfunc0.inX(:,1,1,itr);
            px = ampfunc0.inPxpw(1,1,itr);
    end
    
    % First scattering direction
    % w - sampled direction.
    % w0p - probability of the sampled direction, needed to compute inportance sampling integral correctly. 
    switch floor(smpFlg/10)
        case 1
            % uniform distribution
            w=randn(dim,1,pathsNum); w=w./sqrt(sum(w.^2,1));
            w0p=ones(1,1,pathsNum)./sqrt(2^(dim-1)*pi);
        case 2
            % g0
            meanl = [0;0;1];
            w=smpampfunc_general(meanl, sct_type,ampfunc0,pathsNum);     
            w0p=(evalampfunc_general(meanl'*w,sct_type,ampfunc0,dim));
        case 3
            % multimode
%             [w,w0p] = smpampfunc_multimode(l_1,sct_type,ampfunc0);
        case 6
            w = ampfunc0.inw(:,1,1,itr);
            w0p = ampfunc0.inPxpw(2,1,itr);
    end

    [e_l0, e_l0_ms] = feval(gpuFunc.el, e_l0, e_l0_ms, lambda, tabulatedAmplitudeFunction, x, w, w0p, ...
        l.x, l.y, l.z, dirl.x, dirl.y, dirl.z);
    
    activatedPaths = true(pathsNum,1,'gpuArray');
      
    % number of scattering events
    pL=0;
    
    % intensity loss due to albedo
    weight=albedo;
    
    % begin paths sampling loop
    while 1

        pL=pL+1;

        % calculate the complex volumetric throughput for the last
        % scattering event in case of multiple scattering
        if (pL>1)
            constCont = sqrt(weight./px)*exp(2*pi*1i*rand);
            
            if(isRefocus && refocus.correlationActive)
                u_ff = feval(gpuFunc.ff_mscattering, u0_ff, e_l0_ms, lambda, activatedPaths, ...
                    tabulatedAmplitudeFunction, x, ow, constCont, ...
                    v.x, v.y, v.z, dirv.x, dirv.y, dirv.z) .* ff_legal_idx;
                
                u_ff_reshaped = reshape(permute(u_ff,permute_ff_dims),reshape_ff_dims);
                
                for commonParamNum = 1:1:size(w_l,3)
                    u_wl(:,:,commonParamNum) = u_ff_reshaped(:,:,commonParamNum) * w_l(:,:,commonParamNum);
                end
                
                u_nf = feval(gpuFunc.refocus, u_nf, reshape(u_wl,[prod(ff_v_dims),2,cnf_l_dims]), w_v);
            else
                u_ff = feval(gpuFunc.ff_mscattering, u_ff, e_l0_ms, lambda, activatedPaths, ...
                    tabulatedAmplitudeFunction, x, ow, constCont, ...
                    v.x, v.y, v.z, dirv.x, dirv.y, dirv.z) .* ff_legal_idx;
            end
        end

        % Update field with next-event estimation
        if (pL==1)
            constCont = sqrt(weight./px).*exp(2*pi*1i*rand(1,1,pathsNum));
            
            if(isRefocus && refocus.correlationActive)
                u_ff = feval(gpuFunc.ff_scattering, u0_ff, e_l0, lambda, af_ang_vl, x, constCont, ...
                    v.x, v.y, v.z, dirv.x, dirv.y, dirv.z) .* ff_legal_idx;
                
                u_ff_reshaped = reshape(permute(u_ff,permute_ff_dims),reshape_ff_dims);
                
                for commonParamNum = 1:1:size(w_l,3)
                    u_wl(:,:,commonParamNum) = u_ff_reshaped(:,:,commonParamNum) * w_l(:,:,commonParamNum);
                end
                
                us_nf = feval(gpuFunc.refocus, us_nf, reshape(u_wl,[prod(ff_v_dims),2,cnf_l_dims]), w_v);
            else
                us_ff = feval(gpuFunc.ff_scattering, us_ff, e_l0, lambda, af_ang_vl, x, constCont, ...
                    v.x, v.y, v.z, dirv.x, dirv.y, dirv.z) .* ff_legal_idx;
            end

        end

        % advance to the next scattering event
        if(smpFlg == 66)
            x = ampfunc0.inX(:,pL+1,1,itr);
        else
            d=-log(-rand+1)/(sigt);
            x=x+d*w;
        end

        % move to the next particle if the next scattering event is outside
        % the box
        if(dim == 3)
            outsideIdx = ...
                x(1,1,:) > box_max(1) | x(1,1,:) < box_min(1) | ...
                x(2,1,:) > box_max(2) | x(2,1,:) < box_min(2) | ...
                x(3,1,:) > box_max(3) | x(3,1,:) < box_min(3) ;
        else
            outsideIdx = ...
                x(1,1,:) > box_max(1) | x(1,1,:) < box_min(1) | ...
                x(2,1,:) > box_max(2) | x(2,1,:) < box_min(2) ;
        end
        
        activatedPaths(outsideIdx) = false;
        
        if(~any(activatedPaths))
            break
        end

        % albedo intensity reduction. If the intensity is too low, it is
        % possible that the particle will be killed
        if weight<killThr
            if rand>albedo
                break
            end
        else
            weight=weight*albedo;
        end

        % Sample new scatteing direction
        ow=w;
        
        if(smpFlg == 66)
            w = ampfunc0.inw(:,pL+1,1,itr);
        else
            w=smpampfunc_general(ow, sct_type, ampfunc, pathsNum);
        end

    end
  
end

if(isRefocus && refocus.correlationActive)
    u_nf = u_nf + us_nf;
else
    u_ff = u_ff + us_ff;
end

%% refocus nf field
if(isRefocus && ~refocus.correlationActive)
    u_wl = (reshape(permute(u_ff,permute_ff_dims),reshape_ff_dims) * w_l);
    u_wl = reshape(u_wl,[prod(ff_v_dims),1,nf_l_dims]);            
    u_nf = feval(gpuFunc.refocus, u_nf, u_wl, w_v);
    
    u_wl = (reshape(permute(us_ff,permute_ff_dims),reshape_ff_dims) * w_l);
    u_wl = reshape(u_wl,[prod(ff_v_dims),1,nf_l_dims]);            
    us_nf = feval(gpuFunc.refocus, us_nf, u_wl, w_v);
end


%% Normalization
if(isRefocus)
    if(~refocus.correlationActive)
        u  = u_nf  * sqrt(1/maxItr*V*sigt) * d_ff;
        us = us_nf * sqrt(1/maxItr*V*sigt) * d_ff;
    else
        u  = u_nf  * (1/maxItr*V*sigt) * d_ff^2;
        us = us_nf * (1/maxItr*V*sigt) * d_ff^2;
    end
else
    if(size(l.x,2) == 1)
        u  = u_ff  * sqrt(1/maxItr*V*sigt);
        us = us_ff * sqrt(1/maxItr*V*sigt);
    else
        u  = u_ff  * (1/maxItr*V*sigt);
        us = us_ff * (1/maxItr*V*sigt);
    end
end

u = gather(u);
us = gather(us);

end

