#ifndef PARA_H
#define PARA_H

namespace Para {
inline constexpr int sys_nrow{10};
inline constexpr int sys_ncol{10};
// inline constexpr int cash_margin{10}; // Magic number 10 from Nobuto's code.
inline constexpr int grid_size{sys_nrow * sys_ncol};
inline constexpr float diffusion_rate{0.3};
inline constexpr int scale{10}; // 10x pixel scale
} // namespace Para

#endif
