#ifndef SPGEMM_DEFINITIONS_H
#define SPGEMM_DEFINITIONS_H

#include <limits>

float comp_ratio_maxium = 0.56;
bool USE_COMPRESSION = true;

#define L2_CACHE_SIZE 512*1024

template<typename T>
T inf()
{
    if (std::is_floating_point<T>::value )
        return std::numeric_limits<T>::infinity();
    else
        return std::numeric_limits<T>::max();
}

#define TASK_PER_THREAD_SYM 30
#define TASK_PER_THREAD_NUM 30
#define TASK_PER_THREAD_OTHERS 7

#define STR_LEN 40
#define NUM_LEN 10

#endif
