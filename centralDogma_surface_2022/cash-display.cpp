#include "cash-display.hpp"

extern int scale;

bool CashDisplay::no_other_instantiation = true;

CashDisplay::CashDisplay(const int _window_row, const int _window_col,
                         const std::vector<CashPanelInfo> &_panel_info,
                         const int _scale, const unsigned char _margin_color)
    : window_row(_window_row), window_col(_window_col), window_is_open(false),
      png_is_open(false), margin_color(_margin_color) {
  /* CashDisplay must have only one instantiation at a time */
  if (no_other_instantiation) {
    no_other_instantiation = false;
  } else {
    std::cerr << "CashDisplay::CashDisplay() the number of instantiated "
                 "Display objects cannot be more than two at the same time."
              << std::endl;
    exit(-1);
  }

  /* Allocate memory in data and initialize it. */
  if (window_row == 0 || window_col == 0) {
    std::cerr << "CashDisplay::Data::Data() Error: nrow=ncol=0 is not allowed."
              << std::endl;
    return;
  }
  try {
    data = new unsigned char[window_row * window_col];
  } catch (std::bad_alloc) {
    std::cerr << "CashDisplay::Data::Data() Memory exausted" << std::endl;
    exit(-1);
  }
  memset(data, margin_color, sizeof(unsigned char) * window_row * window_col);

  /* An error check */
  if (_panel_info.empty()) {
    std::cerr << "CashDisplay::CashDisplay() the number of elements in "
                 "panel_info is zero."
              << std::endl;
    exit(-1);
  }

  /* Get each panel information */
  for (std::vector<CashPanelInfo>::const_iterator iter = _panel_info.begin();
       iter != _panel_info.end(); iter++) {
    if (iter->n_row <= 0 || iter->n_col <= 0 || iter->o_row < 0 ||
        iter->o_col < 0) {
      std::cerr << "CashDisplay::CashDisplay() the "
                << _panel_info.size() -
                       static_cast<std::vector<CashPanelInfo>::size_type>(
                           _panel_info.end() - iter)
                << "-th panel_info element has abnormal specification: n_row="
                << iter->n_row << " n_col=" << iter->n_col
                << " o_row=" << iter->o_row << " o_col=" << iter->o_col
                << std::endl;
      exit(-1);
    }

    if (iter->n_row + iter->o_row > window_row ||
        iter->n_col + iter->o_col > window_col) {
      std::cerr << "CashDisplay::CashDisplay() the "
                << _panel_info.size() -
                       static_cast<std::vector<CashPanelInfo>::size_type>(
                           _panel_info.end() - iter)
                << "-th panel_info element's specification violates the window "
                   "size: n_row="
                << iter->n_row << " n_col=" << iter->n_col
                << " o_row=" << iter->o_row << " o_col=" << iter->o_col
                << " window_row=" << window_row << " window_col=" << window_col
                << std::endl;
      exit(-1);
    }

    panels.push_back(*iter);
  }

  /* Set the size of pixel */
  if (_scale > 0) {
    scale = _scale;
  } else {
    std::cerr << "CashDisplay::CashDisplay(): Error, scale=" << _scale
              << std::endl;
    exit(-1);
  }

  /* Set default colors */
  ColorRGB(CashColor::BLACK, 0, 0, 0);
  ColorRGB(CashColor::WHITE, 255, 255, 255);
  ColorRGB(CashColor::RED, 255, 0, 0);
  ColorRGB(CashColor::GREEN, 0, 255, 0);
  ColorRGB(CashColor::BLUE, 0, 0, 255);
  ColorRGB(CashColor::YELLOW, 255, 255, 0);
  ColorRGB(CashColor::BROWN, 188, 143, 143);
  ColorRGB(CashColor::GRAY, 220, 220, 220);
  ColorRGB(CashColor::VIOLET, 148, 0, 211);
  ColorRGB(CashColor::CYAN, 0, 255, 255);
  ColorRGB(CashColor::MAGENTA, 255, 0, 255);
  ColorRGB(CashColor::ORANGE, 255, 165, 0);
  ColorRGB(CashColor::INDIGO, 114, 33, 188);
  ColorRGB(CashColor::MAROON, 103, 7, 72);
  ColorRGB(CashColor::TURQUOISE, 64, 224, 208);
  ColorRGB(CashColor::GREEN4, 0, 139, 0);
}

CashDisplay::~CashDisplay() {
  delete[] data;

  if (window_is_open)
    CloseDisplayImmediately();
  if (png_is_open)
    ClosePNG();

  no_other_instantiation = true;
}

void CashDisplay::color_yellow2yellow(const unsigned char color_ind_begin,
                                      const unsigned char length) {
  int r = 256, g = 255, b = 0;
  double x = 1530. / static_cast<double>(length - 1);
  int y = 0;
  for (int i = 0; i < 1531; ++i) {
    if (i <= 255) {
      r = r - 1;
    } else if (i <= 255 * 2) {
      b = b + 1;
    } else if (i <= 255 * 3) {
      g = g - 1;
    } else if (i <= 255 * 4) {
      r = r + 1;
    } else if (i <= 255 * 5) {
      b = b - 1;
    } else if (i <= 255 * 6) {
      g = g + 1;
    }
    if (i == static_cast<int>(y * x + 0.5)) {
      ColorRGB(y + color_ind_begin, r, g, b);
      y++;
    }
  }
}

void CashDisplay::color_yellow2red(const unsigned char color_ind_begin,
                                   const unsigned char length) {
  int r = 256, g = 255, b = 0;
  double x = 1275. / static_cast<double>(length - 1);
  int y = 0;
  for (int i = 0; i < 1276; ++i) {
    if (i <= 255) {
      r = r - 1;
    } else if (i <= 255 * 2) {
      b = b + 1;
    } else if (i <= 255 * 3) {
      g = g - 1;
    } else if (i <= 255 * 4) {
      r = r + 1;
    } else if (i <= 255 * 5) {
      b = b - 1;
    }
    if (i == static_cast<int>(y * x + 0.5)) {
      ColorRGB(y + color_ind_begin, r, g, b);
      y++;
    }
  }
}

void CashDisplay::put_histogram(
    const unsigned panel_ind, const std::vector<double> &frequencies,
    const double max_y, const unsigned char color_of_blank,
    const unsigned char color_of_background,
    const std::vector<unsigned char> &color_of_bar) {
  try {
    Assert<GeneralError>(!ASSERT::ERROR_CHECK || panel_ind < panels.size());
  } catch (GeneralError) {
    std::cerr << "CashDisplay::put_histogram(): Error, out-of-bound panel_ind="
              << panel_ind << std::endl;
    exit(-1);
  }

  try {
    Assert<GeneralError>(!ASSERT::ERROR_CHECK || (max_y > 0.));
  } catch (GeneralError) {
    std::cerr << "CashDisplay::put_histogram(): Error, max_y=" << max_y
              << std::endl;
    exit(-1);
  }

  int n_pixel_per_bin = panels[panel_ind].n_col / frequencies.size();

  if (ASSERT::ERROR_CHECK) {
    if (n_pixel_per_bin < 1) {
      std::cerr << "CashHistogram::CashHistogram() Error, the number of bins "
                   "must be smaller than n_row."
                << std::endl;
      exit(-1);
    }
  }

  if (ASSERT::ERROR_CHECK) {
    for (std::vector<double>::const_iterator iter = frequencies.begin();
         iter != frequencies.end(); iter++) {
      if (*iter < 0.) {
        std::cerr << "CashDisplay::put_histogram(): Error, frequencies have "
                     "negative value."
                  << std::endl;
        exit(-1);
      }
    }
  }

  int row, col, height;
  int col_limit, anti_height;
  std::vector<unsigned char>::const_iterator color_iter = color_of_bar.begin();

  /* Draw each bin */
  for (int bin = 0; bin < static_cast<int>(frequencies.size()); ++bin) {
    height = static_cast<unsigned>(
        frequencies[bin] / max_y * panels[panel_ind].n_row + 0.5);

    /* If height is out of boundary, set it to maximum */
    if (height > panels[panel_ind].n_row)
      height = panels[panel_ind].n_row;

    anti_height = panels[panel_ind].n_row - height;

    col_limit = (bin + 1) * n_pixel_per_bin + 1;
    /* Draw each column */
    for (col = bin * n_pixel_per_bin + 1; col < col_limit; ++col) {
      /* Draw each row */
      for (row = panels[panel_ind].n_row; row > 0; --row) {
        data_at(panel_ind, row, col) =
            (row > anti_height) ? (*color_iter) : color_of_blank;
      }
    }
    /* Incliment the colors. */
    if (color_iter + 1 < color_of_bar.end())
      ++color_iter;
  }

  /* Fill the rest by blank */
  for (col = static_cast<int>(frequencies.size()) * n_pixel_per_bin + 1;
       col <= panels[panel_ind].n_col; ++col)
    for (row = panels[panel_ind].n_row; row > 0; --row) {
      data_at(panel_ind, row, col) = color_of_background;
    }
}

bool CashDisplay::xy_window_to_rc_panel(int x, int y, unsigned &panel_ind,
                                        int &row, int &col) const {
  x = x / scale + 1;
  y = y / scale + 1;

  for (Panels::const_iterator iter = panels.begin(); iter != panels.end();
       iter++) {
    if (iter->o_row < y && iter->o_row + iter->n_row >= y && iter->o_col < x &&
        iter->o_col + iter->n_col >= x) {
      panel_ind = static_cast<unsigned>(iter - panels.begin());
      row = y - iter->o_row;
      col = x - iter->o_col;
      return true;
    }
  }

  /* if we reach here, there is no panel correspongind to this
     (x,y) coordinate. */
  return false;
}

/******************************************************/

int CashDisplay::get_n_row(const unsigned panel_ind) const {
  try {
    Assert<GeneralError>(!ASSERT::ERROR_CHECK || panel_ind < panels.size());
  } catch (GeneralError) {
    std::cerr << "CashDisplay::get_n_row(): Error, out-of-bound panel_ind="
              << panel_ind << std::endl;
    exit(-1);
  }
  return panels[panel_ind].n_row;
}

int CashDisplay::get_n_col(const unsigned panel_ind) const {
  try {
    Assert<GeneralError>(!ASSERT::ERROR_CHECK || panel_ind < panels.size());
  } catch (GeneralError) {
    std::cerr << "CashDisplay::get_n_col(): Error, out-of-bound panel_ind="
              << panel_ind << std::endl;
    exit(-1);
  }
  return panels[panel_ind].n_col;
}

int CashDisplay::get_o_row(const unsigned panel_ind) const {
  try {
    Assert<GeneralError>(!ASSERT::ERROR_CHECK || panel_ind < panels.size());
  } catch (GeneralError) {
    std::cerr << "CashDisplay::get_o_row(): Error, out-of-bound panel_ind="
              << panel_ind << std::endl;
    exit(-1);
  }
  return panels[panel_ind].o_row;
}

int CashDisplay::get_o_col(const unsigned panel_ind) const {
  try {
    Assert<GeneralError>(!ASSERT::ERROR_CHECK || panel_ind < panels.size());
  } catch (GeneralError) {
    std::cerr << "CashDisplay::get_o_row(): Error, out-of-bound panel_ind="
              << panel_ind << std::endl;
    exit(-1);
  }
  return panels[panel_ind].o_col;
}

void CashDisplay::open_window(char *window_title) {
  OpenDisplay(window_title, window_row, window_col);
  window_is_open = true;
}

void CashDisplay::open_png(char *directory_name) {
  OpenPNG(directory_name, window_row, window_col);
  png_is_open = true;
}

void CashDisplay::color_rgb(const unsigned char color_ind,
                            const unsigned char r, const unsigned char g,
                            const unsigned char b) {
  if (window_is_open || png_is_open) {
    std::cerr << "CashDisplay::color_rgb(): Error, call this method before "
                 "openning window or png"
              << std::endl;
    exit(-1);
  }

  ColorRGB(color_ind, r, g, b);
}
