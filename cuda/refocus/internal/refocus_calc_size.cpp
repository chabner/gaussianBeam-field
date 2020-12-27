#include "refocus_interface.h"

// debug
// #include "mex.h"

void calculate_current_size(ub32 *current_size, calc_size_return *retVal,
        ub32 entry_size, calc_size_return err_val)
{
    if(entry_size != 1)
    {
        if(*current_size == 1)
        {
            *current_size = entry_size;
        }
        else
        {
            if(*current_size != entry_size)
            {
                *retVal = err_val;
            }
        }
    }
}

template<typename T,typename T2,typename T3>
calc_size_return calculate_size(input_st<T,T2,T3> *data_in)
{
    calc_size_return retVal = VALID_SIZE;
    
    ub32 out_dims_num = 0;
    out_dims_num = (out_dims_num > data_in->k.data_dims_num ? out_dims_num : data_in->k.data_dims_num);
    out_dims_num = (out_dims_num > data_in->light_point.data_dims_num_x ? out_dims_num : data_in->light_point.data_dims_num_x);
    out_dims_num = (out_dims_num > data_in->light_point.data_dims_num_y ? out_dims_num : data_in->light_point.data_dims_num_y);
    out_dims_num = (out_dims_num > data_in->light_point.data_dims_num_z ? out_dims_num : data_in->light_point.data_dims_num_z);
    out_dims_num = (out_dims_num > data_in->light_direction.data_dims_num_x ? out_dims_num : data_in->light_direction.data_dims_num_x);
    out_dims_num = (out_dims_num > data_in->light_direction.data_dims_num_y ? out_dims_num : data_in->light_direction.data_dims_num_y);
    out_dims_num = (out_dims_num > data_in->light_direction.data_dims_num_z ? out_dims_num : data_in->light_direction.data_dims_num_z);
    out_dims_num = (out_dims_num > data_in->view_point.data_dims_num_x ? out_dims_num : data_in->view_point.data_dims_num_x);
    out_dims_num = (out_dims_num > data_in->view_point.data_dims_num_y ? out_dims_num : data_in->view_point.data_dims_num_y);
    out_dims_num = (out_dims_num > data_in->view_point.data_dims_num_z ? out_dims_num : data_in->view_point.data_dims_num_z);
    out_dims_num = (out_dims_num > data_in->view_direction.data_dims_num_x ? out_dims_num : data_in->view_direction.data_dims_num_x);
    out_dims_num = (out_dims_num > data_in->view_direction.data_dims_num_y ? out_dims_num : data_in->view_direction.data_dims_num_y);
    out_dims_num = (out_dims_num > data_in->view_direction.data_dims_num_z ? out_dims_num : data_in->view_direction.data_dims_num_z);
    
    if(data_in->is_correlation)
    {
        out_dims_num = (out_dims_num > data_in->k_2.data_dims_num ? out_dims_num : data_in->k_2.data_dims_num);
        out_dims_num = (out_dims_num > data_in->light_point_2.data_dims_num_x ? out_dims_num : data_in->light_point_2.data_dims_num_x);
        out_dims_num = (out_dims_num > data_in->light_point_2.data_dims_num_y ? out_dims_num : data_in->light_point_2.data_dims_num_y);
        out_dims_num = (out_dims_num > data_in->light_point_2.data_dims_num_z ? out_dims_num : data_in->light_point_2.data_dims_num_z);
        out_dims_num = (out_dims_num > data_in->light_direction_2.data_dims_num_x ? out_dims_num : data_in->light_direction_2.data_dims_num_x);
        out_dims_num = (out_dims_num > data_in->light_direction_2.data_dims_num_y ? out_dims_num : data_in->light_direction_2.data_dims_num_y);
        out_dims_num = (out_dims_num > data_in->light_direction_2.data_dims_num_z ? out_dims_num : data_in->light_direction_2.data_dims_num_z);
        out_dims_num = (out_dims_num > data_in->view_point_2.data_dims_num_x ? out_dims_num : data_in->view_point_2.data_dims_num_x);
        out_dims_num = (out_dims_num > data_in->view_point_2.data_dims_num_y ? out_dims_num : data_in->view_point_2.data_dims_num_y);
        out_dims_num = (out_dims_num > data_in->view_point_2.data_dims_num_z ? out_dims_num : data_in->view_point_2.data_dims_num_z);
        out_dims_num = (out_dims_num > data_in->view_direction_2.data_dims_num_x ? out_dims_num : data_in->view_direction_2.data_dims_num_x);
        out_dims_num = (out_dims_num > data_in->view_direction_2.data_dims_num_y ? out_dims_num : data_in->view_direction_2.data_dims_num_y);
        out_dims_num = (out_dims_num > data_in->view_direction_2.data_dims_num_z ? out_dims_num : data_in->view_direction_2.data_dims_num_z);
    }
    
    for(ub32 dim_num = 0; dim_num < out_dims_num; dim_num++)
    {
        ub32 current_size = 1;
        
        calculate_current_size(&current_size,&retVal,data_in->k.data_size[dim_num],ERR_K);
        
        calculate_current_size(&current_size,&retVal,data_in->light_point.data_size_x[dim_num]    ,ERR_LIGHT_P_X);
        calculate_current_size(&current_size,&retVal,data_in->light_point.data_size_y[dim_num]    ,ERR_LIGHT_P_Y);
        calculate_current_size(&current_size,&retVal,data_in->light_point.data_size_z[dim_num]    ,ERR_LIGHT_P_Z);
        calculate_current_size(&current_size,&retVal,data_in->light_direction.data_size_x[dim_num],ERR_LIGHT_D_X);
        calculate_current_size(&current_size,&retVal,data_in->light_direction.data_size_y[dim_num],ERR_LIGHT_D_Y);
        calculate_current_size(&current_size,&retVal,data_in->light_direction.data_size_z[dim_num],ERR_LIGHT_D_Z);
        
        calculate_current_size(&current_size,&retVal,data_in->view_point.data_size_x[dim_num]    ,ERR_VIEW_P_X);
        calculate_current_size(&current_size,&retVal,data_in->view_point.data_size_y[dim_num]    ,ERR_VIEW_P_Y);
        calculate_current_size(&current_size,&retVal,data_in->view_point.data_size_z[dim_num]    ,ERR_VIEW_P_Z);
        calculate_current_size(&current_size,&retVal,data_in->view_direction.data_size_x[dim_num],ERR_VIEW_D_X);
        calculate_current_size(&current_size,&retVal,data_in->view_direction.data_size_y[dim_num],ERR_VIEW_D_Y);
        calculate_current_size(&current_size,&retVal,data_in->view_direction.data_size_z[dim_num],ERR_VIEW_D_Z);
        
        if(data_in->is_correlation)
        {
            calculate_current_size(&current_size,&retVal,data_in->k_2.data_size[dim_num],ERR_K);
            
            calculate_current_size(&current_size,&retVal,data_in->light_point_2.data_size_x[dim_num]    ,ERR_LIGHT_P_X_2);
            calculate_current_size(&current_size,&retVal,data_in->light_point_2.data_size_y[dim_num]    ,ERR_LIGHT_P_Y_2);
            calculate_current_size(&current_size,&retVal,data_in->light_point_2.data_size_z[dim_num]    ,ERR_LIGHT_P_Z_2);
            calculate_current_size(&current_size,&retVal,data_in->light_direction_2.data_size_x[dim_num],ERR_LIGHT_D_X_2);
            calculate_current_size(&current_size,&retVal,data_in->light_direction_2.data_size_y[dim_num],ERR_LIGHT_D_Y_2);
            calculate_current_size(&current_size,&retVal,data_in->light_direction_2.data_size_z[dim_num],ERR_LIGHT_D_Z_2);

            calculate_current_size(&current_size,&retVal,data_in->view_point_2.data_size_x[dim_num]    ,ERR_VIEW_P_X_2);
            calculate_current_size(&current_size,&retVal,data_in->view_point_2.data_size_y[dim_num]    ,ERR_VIEW_P_Y_2);
            calculate_current_size(&current_size,&retVal,data_in->view_point_2.data_size_z[dim_num]    ,ERR_VIEW_P_Z_2);
            calculate_current_size(&current_size,&retVal,data_in->view_direction_2.data_size_x[dim_num],ERR_VIEW_D_X_2);
            calculate_current_size(&current_size,&retVal,data_in->view_direction_2.data_size_y[dim_num],ERR_VIEW_D_Y_2);
            calculate_current_size(&current_size,&retVal,data_in->view_direction_2.data_size_z[dim_num],ERR_VIEW_D_Z_2);
        }
        
        if(retVal == VALID_SIZE)
        {
            data_in->out_size[dim_num] = current_size;
        }
        else
        {
            return retVal;
        }
    }
    
    data_in->out_dims_num = out_dims_num;
    
    // calculate the total number of elements for the output data
    data_in->total_elements = 1;
    
    for(ub32 dimElemNum = 0 ; dimElemNum < data_in->out_dims_num  ; dimElemNum++ )
    {
        data_in->total_elements *= data_in->out_size[dimElemNum];
    }
    
    return retVal;
}

template calc_size_return calculate_size<double, double2, double3>(input_st<double, double2, double3> *data_in);
template calc_size_return calculate_size<float, float2, float3>(input_st<float, float2, float3> *data_in);
