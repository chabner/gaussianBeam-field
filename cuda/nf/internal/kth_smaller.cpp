// Quickselect Algorithm code, taken from https://www.geeksforgeeks.org/quickselect-algorithm/

template<typename T>
void inner_swap(T *x, T *y)
{
    T tmp = *x;
    *x = *y;
    *y = tmp;
}

// Standard partition process of QuickSort(). 
// It considers the last element as pivot 
// and moves all smaller element to left of 
// it and greater elements to right 
template<typename T>
int partition(T *arr, int l, int r) 
{ 
    T x = arr[r];
    int i = l;
    for (int j = l; j <= r - 1; j++) { 
        if (arr[j] <= x) { 
            inner_swap<T>(arr + i, arr + j); 
            i++; 
        } 
    } 
    inner_swap<T>(arr + i, arr + r); 
    return i; 
} 
  
// This function returns k'th smallest  
// element in arr[l..r] using QuickSort  
// based method.  ASSUMPTION: ALL ELEMENTS 
// IN ARR[] ARE DISTINCT 
template<typename T>
T kthSmallest(T *arr, int l, int r, int k) 
{ 
    // If k is smaller than number of  
    // elements in array 
    if (k > 0 && k <= r - l + 1)
    { 
        // Partition the array around last  
        // element and get position of pivot  
        // element in sorted array 
        int index = partition<T>(arr, l, r); 
  
        // If position is same as k 
        if (index - l == k - 1) 
            return arr[index]; 
  
        // If position is more, recur  
        // for left subarray 
        if (index - l > k - 1)  
            return kthSmallest<T>(arr, l, index - 1, k); 
  
        // Else recur for right subarray 
        return kthSmallest<T>(arr, index + 1, r, k - index + l - 1); 
    } 
  
    // If k is more than number of  
    // elements in array 
    return (T) 0.0; 
} 

template<typename T>
T kthSmallest(T *arr, int n, int k)
{
    return kthSmallest<T>(arr,0,n-1,k);
}
