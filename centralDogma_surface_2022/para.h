#ifndef PARA_H
#define PARA_H

namespace Para {
//-----MODEL PARAMETERS
inline unsigned sys_nrow{512};
inline unsigned sys_ncol{512};
inline unsigned grid_size{sys_nrow * sys_ncol};
inline unsigned scale{1}; // 10x pixel scale
inline double decay_probability{0.02};
inline double diffusion_probability{0.1};
inline double mutation_probability{0.5};
inline double alpha{0.9999 / 1.62};
inline double beta{0.5};
inline double gamma{0.5};

//------MISC
inline unsigned visualization{0};
inline unsigned mt_seed{0};
inline unsigned movie{0};
inline unsigned display_interval{1};
inline unsigned max_time{1000000};
} // namespace Para

#endif
