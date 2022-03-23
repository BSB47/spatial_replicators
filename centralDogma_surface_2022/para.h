#ifndef PARA_H
#define PARA_H

namespace Para {
inline constexpr int sys_nrow{512};
inline constexpr int sys_ncol{512};
// inline constexpr int cash_margin{10}; // Magic number 10 from Nobuto's code.
inline constexpr int grid_size{sys_nrow * sys_ncol};
inline constexpr double diffusion_rate{0.3};
inline constexpr int scale{10}; // 10x pixel scale
inline constexpr double decay_probability{0.1};
inline constexpr double alpha{
    0.05 / Para::decay_probability}; // alpha * d needs to be as large as
                                     // possible while being < 1.
inline constexpr long display_interval{1};
inline constexpr long max_time{1000};
} // namespace Para

#endif
