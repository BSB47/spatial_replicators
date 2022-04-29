#ifndef PARA_H
#define PARA_H

namespace Para {
//-----MODEL PARAMETERS
inline constexpr int sys_nrow{512};
inline constexpr int sys_ncol{512};
inline constexpr int grid_size{sys_nrow * sys_ncol};
inline constexpr int scale{10}; // 10x pixel scale
inline constexpr double decay_probability{0.02};
inline constexpr double diffusion_probability{0.1};
inline constexpr double mutation_probability{0.5};
inline constexpr double alpha{0.9999 / 1.27};
inline constexpr double beta{0.25};
inline constexpr double delta{0.5}; // beta but for complexes
inline constexpr double gamma{0.25};

//------MISC
inline unsigned visualization{};
inline unsigned movie{};
inline constexpr long display_interval{100};
inline constexpr long max_time{1000000};
} // namespace Para

#endif
