#include "nf_mex.h"

const float* getDataFromStruct(const mxArray *prhs[], ub32 paramNum, float paramType)
{
    return mxGetSingles(prhs[paramNum]);
}

const double* getDataFromStruct(const mxArray *prhs[], ub32 paramNum, double paramType)
{
    return mxGetDoubles(prhs[paramNum]);
}

const float* getData(const mxArray *in_param, float paramType)
{
    return mxGetSingles(in_param);
}

const double* getData(const mxArray *in_param, double paramType)
{
    return mxGetDoubles(in_param);
}

const float2* getComplexData(const mxArray *in_param, float paramType)
{
    return (float2*) mxGetComplexSingles(in_param);
}

const double2* getComplexData(const mxArray *in_param, double paramType)
{
    return (double2*) mxGetComplexDoubles(in_param);
}

const ub32 getDimsNum(const mxArray *mx_st)
{
    ub32 dimsNum = (ub32) mxGetNumberOfDimensions(mx_st);

    if(dimsNum == 2)
    {
        if(mxGetN(mx_st) == 1)
        {
            dimsNum = 1;
        }
    }

    return dimsNum;
}

void parseError(calc_size_return errVal)
{
    switch(errVal)
    {
        case ERR_K:
            mexErrMsgIdAndTxt("NF:inputSize:k","Wrong input size of k."); break;
        case ERR_LIGHT_P_X:
            mexErrMsgIdAndTxt("NF:inputSize:lpx","Wrong input size of x light position."); break;
        case ERR_LIGHT_P_Y:
            mexErrMsgIdAndTxt("NF:inputSize:lpy","Wrong input size of y light position."); break;
        case ERR_LIGHT_P_Z:
            mexErrMsgIdAndTxt("NF:inputSize:lpz","Wrong input size of z light position."); break;
        case ERR_LIGHT_D_X:
            mexErrMsgIdAndTxt("NF:inputSize:ldx","Wrong input size of x light direction."); break;
        case ERR_LIGHT_D_Y:
            mexErrMsgIdAndTxt("NF:inputSize:ldy","Wrong input size of y light direction."); break;
        case ERR_LIGHT_D_Z:
            mexErrMsgIdAndTxt("NF:inputSize:ldz","Wrong input size of z light direction."); break;
        case ERR_VIEW_P_X:
            mexErrMsgIdAndTxt("NF:inputSize:vpx","Wrong input size of x view position."); break;
        case ERR_VIEW_P_Y:
            mexErrMsgIdAndTxt("NF:inputSize:vpy","Wrong input size of y view position."); break;
        case ERR_VIEW_P_Z:
            mexErrMsgIdAndTxt("NF:inputSize:vpz","Wrong input size of z view position."); break;
        case ERR_VIEW_D_X:
            mexErrMsgIdAndTxt("NF:inputSize:vdx","Wrong input size of x view direction."); break;
        case ERR_VIEW_D_Y:
            mexErrMsgIdAndTxt("NF:inputSize:vdy","Wrong input size of y view direction."); break;
        case ERR_VIEW_D_Z:
            mexErrMsgIdAndTxt("NF:inputSize:vdz","Wrong input size of z view direction."); break;
        case ERR_LIGHT_P_X_2:
            mexErrMsgIdAndTxt("NF:inputSize:lpx2","Wrong input size of x light position of u2."); break;
        case ERR_LIGHT_P_Y_2:
            mexErrMsgIdAndTxt("NF:inputSize:lpy2","Wrong input size of y light position of u2."); break;
        case ERR_LIGHT_P_Z_2:
            mexErrMsgIdAndTxt("NF:inputSize:lpz2","Wrong input size of z light position of u2."); break;
        case ERR_LIGHT_D_X_2:
            mexErrMsgIdAndTxt("NF:inputSize:ldx2","Wrong input size of x light direction of u2."); break;
        case ERR_LIGHT_D_Y_2:
            mexErrMsgIdAndTxt("NF:inputSize:ldy2","Wrong input size of y light direction of u2."); break;
        case ERR_LIGHT_D_Z_2:
            mexErrMsgIdAndTxt("NF:inputSize:ldz2","Wrong input size of z light direction of u2."); break;
        case ERR_VIEW_P_X_2:
            mexErrMsgIdAndTxt("NF:inputSize:vpx2","Wrong input size of x view position of u2."); break;
        case ERR_VIEW_P_Y_2:
            mexErrMsgIdAndTxt("NF:inputSize:vpy2","Wrong input size of y view position of u2."); break;
        case ERR_VIEW_P_Z_2:
            mexErrMsgIdAndTxt("NF:inputSize:vpz2","Wrong input size of z view position of u2."); break;
        case ERR_VIEW_D_X_2:
            mexErrMsgIdAndTxt("NF:inputSize:vdx2","Wrong input size of x view direction of u2."); break;
        case ERR_VIEW_D_Y_2:
            mexErrMsgIdAndTxt("NF:inputSize:vdy2","Wrong input size of y view direction of u2."); break;
        case ERR_VIEW_D_Z_2:
            mexErrMsgIdAndTxt("NF:inputSize:vdz2","Wrong input size of z view direction of u2."); break;
    }
}

void c_to_matlab_array_size(const ub32 *c_ptr, mwSize *matlab_ptr, ub32 dims_num)
{
    for(ub32 curr_dim = 0; curr_dim < dims_num ; curr_dim++)
    {
        matlab_ptr[curr_dim] = c_ptr[curr_dim];
    }
}

template<typename T>
const void setData(input_data_triplet<T>* data_triplet, const mxArray *mx_st, 
        const char* param_name_x, const char* param_name_y, const char* param_name_z)
{
    const mwSize* dim_vec;
    ub32 dims_num;
    T paramType = 0.0;

    const mxArray *param_x = mxGetField(mx_st, 0, param_name_x);
    dims_num = getDimsNum(param_x);
    dim_vec = mxGetDimensions(param_x);
    for(ub32 dimNum = 0; dimNum < MAX_INPUT_DIM ; dimNum++)
    {
        if(dimNum < dims_num)
        {
            data_triplet->data_size_x[dimNum] = (ub32) dim_vec[dimNum];
        }
        else
        {
            data_triplet->data_size_x[dimNum] = 1;
        }
    }
    data_triplet->data_dims_num_x = dims_num;
    data_triplet->data_ptr_x = getData(param_x,paramType);

    const mxArray *param_y = mxGetField(mx_st, 0, param_name_y);
    dims_num = getDimsNum(param_y);
    dim_vec = mxGetDimensions(param_y);
    for(ub32 dimNum = 0; dimNum < MAX_INPUT_DIM ; dimNum++)
    {
        if(dimNum < dims_num)
        {
            data_triplet->data_size_y[dimNum] = (ub32) dim_vec[dimNum];
        }
        else
        {
            data_triplet->data_size_y[dimNum] = 1;
        }
    }
    data_triplet->data_dims_num_y = dims_num;
    data_triplet->data_ptr_y = getData(param_y,paramType);

    const mxArray *param_z = mxGetField(mx_st, 0, param_name_z);
    dims_num = getDimsNum(param_z);
    dim_vec = mxGetDimensions(param_z);
    for(ub32 dimNum = 0; dimNum < MAX_INPUT_DIM ; dimNum++)
    {
        if(dimNum < dims_num)
        {
            data_triplet->data_size_z[dimNum] = (ub32) dim_vec[dimNum];
        }
        else
        {
            data_triplet->data_size_z[dimNum] = 1;
        }
    }
    data_triplet->data_dims_num_z = dims_num;
    data_triplet->data_ptr_z = getData(param_z,paramType);
}

template<typename T>
const void setData(input_data<T>* data_in, const mxArray *mx_st, 
        const char* param_name)
{
    const mwSize* dim_vec;
    ub32 dims_num;
    T paramType = 0.0;

    const mxArray *param = mxGetField(mx_st, 0, param_name);
    dims_num = getDimsNum(param);
    dim_vec = mxGetDimensions(param);
    for(ub32 dimNum = 0; dimNum < MAX_INPUT_DIM ; dimNum++)
    {
        if(dimNum < dims_num)
        {
            data_in->data_size[dimNum] = (ub32) dim_vec[dimNum];
        }
        else
        {
            data_in->data_size[dimNum] = 1;
        }
    }
    data_in->data_dims_num = dims_num;
    data_in->data_ptr = getData(param,paramType);
}

template<typename T, typename T2, typename T3>
input_st<T,T2,T3>* initMex(int nrhs, const mxArray *prhs[])
{
    T paramType = 0.0;
    ub32 corrFlag = (mxGetField(prhs[ILLUMINATION], 0, "P1_2") != NULL);

    input_st<T,T2,T3> *dataHandler = (input_st<T,T2,T3> *)malloc(sizeof(input_st<T,T2,T3>));
    
    dataHandler->cuda_device_num = *mxGetUint32s(mxGetField(prhs[SIMULATION], 0, "gpuNum"));
    dataHandler->precision = (ub32) mxIsDouble(mxGetField(prhs[MEDIUM], 0, "sigt"));
    dataHandler->is_correlation = corrFlag;
    dataHandler->iterations_num = *mxGetUint32s(mxGetField(prhs[SIMULATION], 0, "iterations"));
    
    ub32 vector_dims_num = (ub32) (mxGetN(mxGetField(prhs[MEDIUM], 0, "box_min")) *  mxGetM(mxGetField(prhs[MEDIUM], 0, "box_min")));
    
    if(vector_dims_num != 3)
    {
        mexErrMsgIdAndTxt("NF:dims:number","vector dims number in NF must be 3D.");
    }

    dataHandler->box_min.x = *getData(mxGetField(prhs[MEDIUM], 0, "box_min"),paramType);
    dataHandler->box_min.y = *(getData(mxGetField(prhs[MEDIUM], 0, "box_min"),paramType) + 1);
    dataHandler->box_min.z = *(getData(mxGetField(prhs[MEDIUM], 0, "box_min"),paramType) + 2);

    dataHandler->box_max.x = *getData(mxGetField(prhs[MEDIUM], 0, "box_max"),paramType);
    dataHandler->box_max.y = *(getData(mxGetField(prhs[MEDIUM], 0, "box_max"),paramType) + 1);
    dataHandler->box_max.z = *(getData(mxGetField(prhs[MEDIUM], 0, "box_max"),paramType) + 2);

    dataHandler->sigt = *getData(mxGetField(prhs[MEDIUM], 0, "sigt"),paramType);

    // scattering
    dataHandler->scattering_input.type = *mxGetUint32s(mxGetField(prhs[SCATTERING], 0, "type"));
    if(dataHandler->scattering_input.type > 3)
    {
        mexErrMsgIdAndTxt("NF:scattering:flag","Invalid scattering flag.");
    }
    if(dataHandler->scattering_input.type == 1)
    {
        dataHandler->scattering_input.tabulated_values = 0;
    }
    if(dataHandler->scattering_input.type == 2)
    {
        dataHandler->scattering_input.tabulated_entries = (ub32)(
                mxGetN(mxGetField(prhs[SCATTERING], 0, "f")) * 
                mxGetM(mxGetField(prhs[SCATTERING], 0, "f")));
        
        if(dataHandler->scattering_input.tabulated_entries > 1048576)
        {
            mexErrMsgIdAndTxt("NF:scattering:tabulated","maximal tabulated entries allowed is 1048576.");
        }
        
        dataHandler->scattering_input.tabulated_values = getComplexData(mxGetField(prhs[SCATTERING], 0, "f"), paramType);
    }
    if(dataHandler->scattering_input.type == 3)
    {
        dataHandler->scattering_input.g = *getData(mxGetField(prhs[SCATTERING], 0, "g"),paramType);
        dataHandler->scattering_input.tabulated_values = 0;
        
        if(fabs(dataHandler->scattering_input.g) > 1)
        {
            mexErrMsgIdAndTxt("NF:scattering:g","g value must be smaller than 1.");
        }
    }
    
    // movmf
    
    // get scattering mixture size (all mixtures have only 1 dim)
    dataHandler->scattering_vMF_mixture.mixture_size = (ub32) mxGetM(mxGetField(prhs[MOVMF], 0, "alpha"));

    dataHandler->scattering_vMF_mixture.mixture_alpha = getData(mxGetField(prhs[MOVMF], 0, "alpha"),paramType);
    dataHandler->scattering_vMF_mixture.mixture_mu = getData(mxGetField(prhs[MOVMF], 0, "mu3"),paramType);
    dataHandler->scattering_vMF_mixture.mixture_c = getData(mxGetField(prhs[MOVMF], 0, "c"),paramType);

    // aperture
    
    dataHandler->aperture_kappa_l = *getData(mxGetField(prhs[APERTURE], 0, "kappa_l"),paramType);
    dataHandler->aperture_c_l = *getData(mxGetField(prhs[APERTURE], 0, "c_l"),paramType);
    dataHandler->aperture_kappa_v = *getData(mxGetField(prhs[APERTURE], 0, "kappa_v"),paramType);
    dataHandler->aperture_c_v = *getData(mxGetField(prhs[APERTURE], 0, "c_v"),paramType);

     // sampling
    dataHandler->iS.position_type = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "position_type"));
    dataHandler->iS.direction_type = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "direction_type"));
    dataHandler->iS.is_same_beam = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "same_beam"));

    if(dataHandler->iS.position_type == 4)
    {
        dataHandler->iS.z0_sample_num = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "z0_sample_num"));
        dataHandler->iS.z_sample_num = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "z_sample_num"));
        dataHandler->iS.z0_samples = getData(mxGetField(prhs[SAMPLING], 0, "z0_samples"),paramType);
        
        if(dataHandler->iS.z_sample_num > 1024)
        {
            mexErrMsgIdAndTxt("NF:IS:z_smp","maximal value of z samples number is 1024, which is 1048576 samples.");
        }
    }
    
    dataHandler->iS.is_min_known = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "is_min_known"));
    
    if(dataHandler->iS.is_min_known)
    {
        dataHandler->iS.min_px0 = *getData(mxGetField(prhs[SAMPLING], 0, "min_px0"),paramType);
        dataHandler->iS.min_pw0 = *getData(mxGetField(prhs[SAMPLING], 0, "min_pw0"),paramType);
    }
    else
    {
        dataHandler->iS.sample_rounds = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "test_rounds"));
        dataHandler->iS.min_percent = *getData(mxGetField(prhs[SAMPLING], 0, "min_probability_percent"),paramType);
    }

    // pixel data
    setData<T>(&dataHandler->k, prhs[WAVENUMBER], "k");
	setData<T>(&dataHandler->light_point    , prhs[ILLUMINATION], "P1", "P2", "P3");
    setData<T>(&dataHandler->light_direction, prhs[ILLUMINATION], "D1", "D2", "D3");
    setData<T>(&dataHandler->view_point     , prhs[VIEW]        , "P1", "P2", "P3");
    setData<T>(&dataHandler->view_direction , prhs[VIEW]        , "D1", "D2", "D3");

    if(corrFlag)
    {
        setData<T>(&dataHandler->k_2, prhs[WAVENUMBER], "k_2");
        setData<T>(&dataHandler->light_point_2    , prhs[ILLUMINATION], "P1_2", "P2_2", "P3_2");
        setData<T>(&dataHandler->light_direction_2, prhs[ILLUMINATION], "D1_2", "D2_2", "D3_2");
        setData<T>(&dataHandler->view_point_2     , prhs[VIEW]        , "P1_2", "P2_2", "P3_2");
        setData<T>(&dataHandler->view_direction_2 , prhs[VIEW]        , "D1_2", "D2_2", "D3_2");
    }

	return dataHandler;
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Check input
    if(nrhs != 9) {
        mexErrMsgIdAndTxt("NF:arrayProduct:nrhs",
                          "9 inputs required.");
        
        /*
         * Inputs:
         * 0 - Simulation
         * 1 - Medium
         * 2 - Aperture
         * 3 - sampling
         * 4 - Scattering
         * 5 - Mixture of vMF
         * 6 - Wavenumber
         * 7 - Illumination
         * 8 - View
         */
    }

    // true for double, false for floats
    bool precisionType = mxIsDouble(mxGetField(prhs[MEDIUM], 0, "sigt"));
    mwSize dims[MAX_INPUT_DIM];
    
    mxArray *u, *us, *norm_factor, *total_iterations, *min_px0, *min_pw0;

    total_iterations = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    norm_factor = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
    
    if(precisionType)
    {
        input_st<double,double2,double3>* dataHandler;
        dataHandler = initMex<double,double2,double3>(nrhs,prhs);
        
        calc_size_return errVal = calculate_size<double,double2,double3>(dataHandler);
        parseError(errVal);
        
        c_to_matlab_array_size(dataHandler->out_size,dims,dataHandler->out_dims_num);

        u = mxCreateNumericArray(dataHandler->out_dims_num, dims,
            mxDOUBLE_CLASS, mxCOMPLEX);
        
        us = mxCreateNumericArray(dataHandler->out_dims_num, dims,
            mxDOUBLE_CLASS, mxCOMPLEX);
        
        min_px0 = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        min_pw0 = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);

        nf_run<double,double2,double3>(dataHandler,(double2 *)mxGetComplexDoubles(u),(double2 *) mxGetComplexDoubles(us),
                (double *) mxGetDoubles(norm_factor), (ub64 *) mxGetUint64s(total_iterations),
                (double *) mxGetDoubles(min_px0), (double *) mxGetDoubles(min_pw0));

        free(dataHandler);
    }
    else
    {
        input_st<float,float2,float3>* dataHandler;
        dataHandler = initMex<float,float2,float3>(nrhs,prhs);
        
        calc_size_return errVal = calculate_size<float,float2,float3>(dataHandler);
        parseError(errVal);
        
        c_to_matlab_array_size(dataHandler->out_size,dims,dataHandler->out_dims_num);
                
        u = mxCreateNumericArray(dataHandler->out_dims_num, dims,
            mxSINGLE_CLASS, mxCOMPLEX);
        
        us = mxCreateNumericArray(dataHandler->out_dims_num, dims,
            mxSINGLE_CLASS, mxCOMPLEX);
        
        min_px0 = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
        min_pw0 = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);
        
        nf_run<float,float2,float3>(dataHandler,(float2 *)mxGetComplexSingles(u),(float2 *)mxGetComplexSingles(us),
                (double *) mxGetDoubles(norm_factor), (ub64 *) mxGetUint64s(total_iterations),
                (float *) mxGetSingles(min_px0), (float *) mxGetSingles(min_pw0) );

        free(dataHandler);
    }
    
    plhs[0] = u;
    plhs[1] = us;
    plhs[2] = norm_factor;
    plhs[3] = total_iterations;
    plhs[4] = min_px0;
    plhs[5] = min_pw0;
}
