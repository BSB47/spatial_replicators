#ifndef PARA_H
#define PARA_H

namespace Para {
inline constexpr int sys_nrow{512};
inline constexpr int sys_ncol{512};
// inline constexpr int cash_margin{10}; // Magic number 10 from Nobuto's code.
inline constexpr int grid_size{sys_nrow * sys_ncol};
inline constexpr float diffusion_rate{0.3};
} // namespace Para

#endif
