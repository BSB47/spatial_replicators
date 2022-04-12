#ifndef PARA_H
#define PARA_H

namespace Para {
inline constexpr int sys_nrow{512};
inline constexpr int sys_ncol{512};
// inline constexpr int cash_margin{10}; // Magic number 10 from Nobuto's code.
inline constexpr int grid_size{sys_nrow * sys_ncol};
inline constexpr double diffusion_rate{0.3};
inline constexpr int scale{10}; // 10x pixel scale
inline constexpr double decay_probability{0.5};
inline constexpr double diffusion_probability{0.1};
inline constexpr double alpha{0.9999 /
                              2.02}; // 2.02 is the maximum sum of rates, which
                                     // belongs to simple molecules
inline constexpr double beta{0.5};
inline constexpr long display_interval{5};
inline constexpr long max_time{1000};
} // namespace Para

#endif
