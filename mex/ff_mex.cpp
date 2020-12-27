#include "ff_mex.h"

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
            mexErrMsgIdAndTxt("FF:inputSize:k","Wrong input size of k."); break;
        case ERR_LIGHT_P_X:
            mexErrMsgIdAndTxt("FF:inputSize:lpx","Wrong input size of x light direction."); break;
        case ERR_LIGHT_P_Y:
            mexErrMsgIdAndTxt("FF:inputSize:lpy","Wrong input size of y light direction."); break;
        case ERR_LIGHT_P_Z:
            mexErrMsgIdAndTxt("FF:inputSize:lpz","Wrong input size of z light direction."); break;
        case ERR_LIGHT_D_X:
            mexErrMsgIdAndTxt("FF:inputSize:ldx","Wrong input size of x light attenuation direction."); break;
        case ERR_LIGHT_D_Y:
            mexErrMsgIdAndTxt("FF:inputSize:ldy","Wrong input size of y light attenuation direction."); break;
        case ERR_LIGHT_D_Z:
            mexErrMsgIdAndTxt("FF:inputSize:ldz","Wrong input size of z light attenuation direction."); break;
        case ERR_VIEW_P_X:
            mexErrMsgIdAndTxt("FF:inputSize:vpx","Wrong input size of x view direction."); break;
        case ERR_VIEW_P_Y:
            mexErrMsgIdAndTxt("FF:inputSize:vpy","Wrong input size of y view direction."); break;
        case ERR_VIEW_P_Z:
            mexErrMsgIdAndTxt("FF:inputSize:vpz","Wrong input size of z view direction."); break;
        case ERR_VIEW_D_X:
            mexErrMsgIdAndTxt("FF:inputSize:vdx","Wrong input size of x view attenuation direction."); break;
        case ERR_VIEW_D_Y:
            mexErrMsgIdAndTxt("FF:inputSize:vdy","Wrong input size of y view attenuation direction."); break;
        case ERR_VIEW_D_Z:
            mexErrMsgIdAndTxt("FF:inputSize:vdz","Wrong input size of z view attenuation direction."); break;
        case ERR_LIGHT_P_X_2:
            mexErrMsgIdAndTxt("FF:inputSize:lpx2","Wrong input size of x light direction of u2."); break;
        case ERR_LIGHT_P_Y_2:
            mexErrMsgIdAndTxt("FF:inputSize:lpy2","Wrong input size of y light direction of u2."); break;
        case ERR_LIGHT_P_Z_2:
            mexErrMsgIdAndTxt("FF:inputSize:lpz2","Wrong input size of z light direction of u2."); break;
        case ERR_LIGHT_D_X_2:
            mexErrMsgIdAndTxt("FF:inputSize:ldx2","Wrong input size of x light attenuation direction of u2."); break;
        case ERR_LIGHT_D_Y_2:
            mexErrMsgIdAndTxt("FF:inputSize:ldy2","Wrong input size of y light attenuation direction of u2."); break;
        case ERR_LIGHT_D_Z_2:
            mexErrMsgIdAndTxt("FF:inputSize:ldz2","Wrong input size of z light attenuation direction of u2."); break;
        case ERR_VIEW_P_X_2:
            mexErrMsgIdAndTxt("FF:inputSize:vpx2","Wrong input size of x view direction of u2."); break;
        case ERR_VIEW_P_Y_2:
            mexErrMsgIdAndTxt("FF:inputSize:vpy2","Wrong input size of y view direction of u2."); break;
        case ERR_VIEW_P_Z_2:
            mexErrMsgIdAndTxt("FF:inputSize:vpz2","Wrong input size of z view direction of u2."); break;
        case ERR_VIEW_D_X_2:
            mexErrMsgIdAndTxt("FF:inputSize:vdx2","Wrong input size of x view attenuation direction of u2."); break;
        case ERR_VIEW_D_Y_2:
            mexErrMsgIdAndTxt("FF:inputSize:vdy2","Wrong input size of y view attenuation direction of u2."); break;
        case ERR_VIEW_D_Z_2:
            mexErrMsgIdAndTxt("FF:inputSize:vdz2","Wrong input size of z view attenuation direction of u2."); break;
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
const void setData(input_data_triplet<T>* data_triplet, const mxArray *mx_st, bool is3D, 
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

    if(is3D)
    {
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
}

template<typename T>
const void setData(input_data<T>* data_in, const mxArray *mx_st, bool is3D, 
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
    ub32 corrFlag = (mxGetField(prhs[ILLUMINATION], 0, "X_2") != NULL);
    
    input_st<T,T2,T3> *dataHandler = (input_st<T,T2,T3> *)malloc(sizeof(input_st<T,T2,T3>));
    
    dataHandler->cuda_device_num = *mxGetUint32s(mxGetField(prhs[SIMULATION], 0, "gpuNum"));
    dataHandler->precision = (ub32) mxIsDouble(mxGetField(prhs[MEDIUM], 0, "sigt"));
    dataHandler->vector_dims_num = (ub32) (mxGetN(mxGetField(prhs[MEDIUM], 0, "box_min")) *  mxGetM(mxGetField(prhs[MEDIUM], 0, "box_min")));
    dataHandler->is_correlation = corrFlag;
    dataHandler->iterations_num = *mxGetUint32s(mxGetField(prhs[SIMULATION], 0, "iterations"));
    dataHandler->is_cbs = *mxGetUint32s(mxGetField(prhs[SIMULATION], 0, "cbs"));
    
    bool is3D = (dataHandler->vector_dims_num == 3);

    dataHandler->box_min.x = *getData(mxGetField(prhs[MEDIUM], 0, "box_min"),paramType);
    dataHandler->box_min.y = *(getData(mxGetField(prhs[MEDIUM], 0, "box_min"),paramType) + 1);
    if(is3D)
    {
        dataHandler->box_min.z = *(getData(mxGetField(prhs[MEDIUM], 0, "box_min"),paramType) + 2);
    }

    dataHandler->box_max.x = *getData(mxGetField(prhs[MEDIUM], 0, "box_max"),paramType);
    dataHandler->box_max.y = *(getData(mxGetField(prhs[MEDIUM], 0, "box_max"),paramType) + 1);
    if(is3D)
    {
        dataHandler->box_max.z = *(getData(mxGetField(prhs[MEDIUM], 0, "box_max"),paramType) + 2);
    }

    dataHandler->sigt = *getData(mxGetField(prhs[MEDIUM], 0, "sigt"),paramType);

    // scattering
    dataHandler->scattering_input.type = *mxGetUint32s(mxGetField(prhs[SCATTERING], 0, "type"));
    if(dataHandler->scattering_input.type > 3)
    {
        mexErrMsgIdAndTxt("FF:scattering:type","Invalid scattering type.");
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
            mexErrMsgIdAndTxt("FF:scattering:tabulated","maximal tabulated entries allowed is 1048576.");
        }
        
        dataHandler->scattering_input.tabulated_values = getComplexData(mxGetField(prhs[SCATTERING], 0, "f"), paramType);
    }
    if(dataHandler->scattering_input.type == 3)
    {
        dataHandler->scattering_input.g = *getData(mxGetField(prhs[SCATTERING], 0, "g"),paramType);
        dataHandler->scattering_input.tabulated_values = 0;
        
        if(fabs(dataHandler->scattering_input.g) > 1)
        {
            mexErrMsgIdAndTxt("FF:scattering:g","g value must be smaller than 1.");
        }
    }
    
    // sampling
    dataHandler->iS.position_type = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "position_type"));
    dataHandler->iS.direction_type = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "direction_type"));
    dataHandler->iS.is_mean = *mxGetUint32s(mxGetField(prhs[SAMPLING], 0, "mean_l"));

    
    if(dataHandler->iS.direction_type == 2)
    {
        if(mxGetField(prhs[SAMPLING], 0, "f0") != NULL)
        {
            dataHandler->iS.tabulated_entries = (ub32)(
                mxGetN(mxGetField(prhs[SAMPLING], 0, "f0")) * 
                mxGetM(mxGetField(prhs[SAMPLING], 0, "f0")));
            
            if(dataHandler->iS.tabulated_entries > 1048576)
            {
                mexErrMsgIdAndTxt("FF:sampling:tabulated","maximal tabulated entries allowed is 1048576.");
            }
            
            dataHandler->iS.tabulated_values = getComplexData(mxGetField(prhs[SAMPLING], 0, "f0"), paramType);
        }
        else
        {
            dataHandler->iS.tabulated_entries = dataHandler->scattering_input.tabulated_entries;
            dataHandler->iS.tabulated_values = dataHandler->scattering_input.tabulated_values;
        }
    }
    
    if(dataHandler->iS.direction_type == 3)
    {
        if(mxGetField(prhs[SAMPLING], 0, "g0") != NULL)
        {
            dataHandler->iS.g0 = *getData(mxGetField(prhs[SAMPLING], 0, "g0"),paramType);
        }
        else
        {
            dataHandler->iS.g0 = dataHandler->scattering_input.g;
        }
    }
    
    // pixel data
    setData<T>(&dataHandler->k, prhs[WAVENUMBER], is3D, "k");
	setData<T>(&dataHandler->light_direction             , prhs[ILLUMINATION], is3D, "X"    , "Y"    , "Z"    );
    setData<T>(&dataHandler->light_attenuation_direction , prhs[ILLUMINATION], is3D, "DIR_X", "DIR_Y", "DIR_Z");
    setData<T>(&dataHandler->view_direction              , prhs[VIEW]        , is3D, "X"    , "Y"    , "Z"    );
    setData<T>(&dataHandler->view_attenuation_direction  , prhs[VIEW]        , is3D, "DIR_X", "DIR_Y", "DIR_Z");

    if(corrFlag)
    {
        setData<T>(&dataHandler->k_2, prhs[WAVENUMBER], is3D, "k_2");
        setData<T>(&dataHandler->light_direction_2             , prhs[ILLUMINATION], is3D, "X_2"    , "Y_2"    , "Z_2"    );
        setData<T>(&dataHandler->light_attenuation_direction_2 , prhs[ILLUMINATION], is3D, "DIR_X_2", "DIR_Y_2", "DIR_Z_2");
        setData<T>(&dataHandler->view_direction_2              , prhs[VIEW]        , is3D, "X_2"    , "Y_2"    , "Z_2"    );
        setData<T>(&dataHandler->view_attenuation_direction_2  , prhs[VIEW]        , is3D, "DIR_X_2", "DIR_Y_2", "DIR_Z_2");
    }

	return dataHandler;
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Check input
    if(nrhs != 7) {
        mexErrMsgIdAndTxt("FF:arrayProduct:nrhs",
                          "7 inputs required.");
        
        /*
         * Inputs:
         * 0 - Simulation
         * 1 - Medium
         * 2 - sampling
         * 3 - Scattering
         * 4 - Wavenumber
         * 5 - Illumination
         * 6 - View
         */
    }

    // true for double, false for floats
    bool precisionType = mxIsDouble(mxGetField(prhs[MEDIUM], 0, "sigt"));
    mwSize dims[MAX_INPUT_DIM];
    
    mxArray *u, *us, *norm_factor, *total_iterations;

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
        
        ff_run<double,double2,double3>(dataHandler,(double2 *)mxGetComplexDoubles(u),(double2 *) mxGetComplexDoubles(us),
                (double *) mxGetDoubles(norm_factor), (ub64 *) mxGetUint64s(total_iterations) );

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
                
        ff_run<float,float2,float3>(dataHandler,(float2 *)mxGetComplexSingles(u),(float2 *)mxGetComplexSingles(us),
                (double *) mxGetDoubles(norm_factor), (ub64 *) mxGetUint64s(total_iterations) );

        free(dataHandler);
    }
    
    plhs[0] = u;
    plhs[1] = us;
    plhs[2] = norm_factor;
    plhs[3] = total_iterations;
}
