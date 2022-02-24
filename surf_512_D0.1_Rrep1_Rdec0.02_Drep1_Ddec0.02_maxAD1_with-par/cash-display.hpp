#include <functional>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <cstring>

#include "cash.h"
#include "assert.hpp"

#ifndef CASH_DISPLAY
#define CASH_DISPLAY

/* CashDisplay allows multiple panels to be displayed in the window
   (there can be only one window). CashPanelInfo is used to specify each
   panel regarding its size and the offset. CashPanelInfor is given to
   CashDisplay, which uses the information and create window and draw
   pixels etc.  */
struct CashPanelInfo{
  /* These specifie the size of the panel */
  int n_row;
  int n_col;

  /* The coordinate of the offset in drawing this panel. It specifies
     the coordinate of the upper left corner [1+o_row][1+o_col]. */
  int o_row;
  int o_col;
};

/* Since cash.h uses macros, I have to undo it to use these names. */
#undef BLACK
#undef WHITE
#undef RED
#undef GREEN
#undef BLUE
#undef YELLOW
#undef MAGENTA
#undef GRAY

/* These are the set of predefined colors. One can use it, but don't
   have to use it. */
namespace CashColor{
  const unsigned char BLACK = 0;
  const unsigned char WHITE = 1;
  const unsigned char RED = 2;
  const unsigned char GREEN = 3;
  const unsigned char BLUE = 4;
  const unsigned char YELLOW = 5;
  const unsigned char BROWN = 6;
  const unsigned char GRAY = 7;
  const unsigned char VIOLET = 8;
  const unsigned char CYAN = 9;
  const unsigned char MAGENTA = 10;
  const unsigned char ORANGE = 11;
  const unsigned char INDIGO = 12;
  const unsigned char MAROON = 13;
  const unsigned char TURQUOISE = 14;
  const unsigned char GREEN4 = 15;
}


/*********************************************************/
/* CashDisplay can have only one instantiation at a time. When
   instantiated, it first defines colors and then open a window (colors
   have to be defined before opening a window). */
class CashDisplay{

  /* This is used to assert that there is no other instantiation of this
     object in the same rung of the program. This is necessary because
     of the design of CASH */
  static bool no_other_instantiation;
  
  /* The size of the window. Window should be large enough to contain
     all panels. */
  int window_row;
  int window_col;

  /* The vector holding information of panels. */
  typedef std::vector<CashPanelInfo> Panels;
  Panels panels;
  
  /* This is the actual data that is handed over to Cash when drawing
     the window or a png file. */
  unsigned char* data;

  /* Flags to know if the window and png is opened */
  bool window_is_open;
  bool png_is_open;

  /* The color of the magins between panels */
  unsigned char margin_color;

public:
  /* Constructer must get basic information about window and panels. The
     size of _panel_info is identical to the number of panels. Each
     element of _panel_info should have right information about that
     panel. */
  CashDisplay(const int _window_row,const int _window_col,const std::vector<CashPanelInfo>& _panel_info,const int _scale=1,const unsigned char _margin_color=0);
  ~CashDisplay();

  /* Setting colors */
  /* This simply calls ColorRGB() of CASH */
  void color_rgb(const unsigned char color_ind,const unsigned char r,const unsigned char g,const unsigned char b);
  /* This generates a particular kind of color gradient from yellow to
     yellow. Color_ind_begin specifies the begining of the color indexes
     that are used for the gradient. 'Length' specifies the number of
     colors used for the gradient. This color gradient is useful for
     drawing histogram with bars with many colors, where each color
     corresponds to individuals in the panel showing a simulation. */
  void color_yellow2yellow(const unsigned char color_ind_begin,const unsigned char length);
  void color_yellow2red(const unsigned char color_ind_begin,const unsigned char length);
  
  /* To open a window and make a png directory. They call OpenWindow()
     and OpenPNG() of CASH  */
  void open_window(char* window_title="cash_display");
  void open_png(char* directory_name="cash_movie");

  /* To get information about the number of panels. This is identical to
     the size of _panel_info given to the constructer.  */
  unsigned get_n_panels() const{return panels.size();}

  /* These give information about different panels. The information is
     identical to that given to the constructor via _panel_info. */
  int get_n_row(const unsigned panel_ind) const;
  int get_n_col(const unsigned panel_ind) const;
  int get_o_row(const unsigned panel_ind) const;
  int get_o_col(const unsigned panel_ind) const;

  /* Setting pixels of panels. This method only sets a pixel of panel,
     but does not draw it on the window. After setting pixels in panels,
     one must call draw_window() to actually draw pixels into the
     window. */
  inline void put_pixel(const unsigned panel_ind,const int row,const int col,const unsigned char color);


  /* This draws a histgram in a panel with given frequencies. The height
     of bars is specified by the argument "frequency". The maximum
     height is obtained when this value equals to the argument max_y. In
     other words, the bar that is longer than max_y is truncated. The
     color of every bar is spcified by the argument color_of_bar.  If
     one wants to use a single color for all bars, use the next version
     of the same method. */
  void put_histogram(const unsigned panel_ind,const std::vector<double>& frequencies,const double max_y,const unsigned char color_of_blank,const unsigned char color_of_background,const std::vector<unsigned char>& color_of_bar);
  /* This is a version of put_histogram() that set the color of all bars
     to one color. */
  inline void put_histogram(const unsigned panel_ind,const std::vector<double>& frequencies,const double max_y,const unsigned char color_of_blank,const unsigned char color_of_background,const unsigned char color_of_bar);

  /* This actually draw pixels into window and png. Pixels should be
     specified by put_pixel(). */
  inline void draw_window();
  inline void draw_png();
  /* This fills the whole window by setting every pixel to the specified
     color. This can also be used to modify the color of the margins
     between the panels. */
  //  inline void fill_window(const unsigned char color_ind);

  /* Mouse handling */
  /* Button 1 is to pause a simulation. The action of Button 2 and 3 can
     be set by the arguments. Note that (x,y) will be the x and y
     coordinates (in pixels) from the top-left of the window (0,0). This
     means that x and y will contain margin and are scaled by scale. Use
     xy_window_to_rc_panel() to get the coordinate in panels.

     To use this function in object oriented programing is rather
     trickily. Let us assume that one defines an objectthat contains a
     CellularAutomata and CashDisplay. Let me call this a simulation
     object. It seems natural that one defines button2() and button3() as
     a method of the simulation object, because this allows button2() and
     button3() to have direct access to CellularAutomata. But, according
     to C++, the pointer to a method of an object instance is different
     from the usual C function pointer. Hence, unfortunately, button2()
     and button3() defined as a method of the simulation object can not
     directly be given to handle_mouse_event(). To solve this, one must
     define a wrapper, a normal C function which calls either button2() or
     button3() "for the focal object instant", and one gives a pointer to
     this wrapper to handle_mouse_event.
  */
  void handle_mouse_event(void (*button2)(int x,int y),void (*button3)(int x,int y)){HandleMouseEvent(button2,button3);}

  /* This method converts (x,y) coordinate of the window into (row,col)
     coordinate of the panel [the upper left corner is (1,1)]. Panel_ind
     gives which panel this (x,y) coordinate belongs to. If there are
     overlapping panels (which is unrecommended), the panel with the
     smallest index has priority. If there is no panel that corresponds
     to this (x,y) coordinate, the method returns false; otherwise, it
     returns true.  */
  bool xy_window_to_rc_panel(int x,int y,unsigned& panel_ind,int& row,int& col)const ;

  /* This method set nframes of png.c to new_nframes. This can be used
     when you are continuing a simulation. After calling this, the next
     png file created will be new_nframes.png. */
  void reset_movie_frame(const long new_nframes){ResetNFrames(new_nframes);}

private:
  CashDisplay(const CashDisplay& rhs){}
  void operator=(const CashDisplay& rhs){}
  inline unsigned char& data_at(const unsigned panel_ind,const unsigned row,const unsigned col);
};

void CashDisplay::put_pixel(const unsigned panel_ind,const int row,const int col,const unsigned char color)
{
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||panel_ind < panels.size());
  }catch(GeneralError){
    std::cerr << "CashDisplay::put_pixel(): Error, out-of-bound panel_ind=" << panel_ind << std::endl;
    exit(-1);
  }
  
  try{
    Assert<GeneralError>(!ASSERT::ERROR_CHECK||(row>=1 && row<=panels[panel_ind].n_row && col>=1 && col<=panels[panel_ind].n_col));
  }catch(GeneralError){
    std::cerr << "CashDisplay::put_pixel(): Error, out-of-bound row=" << row << " col=" << col << std::endl;
    exit(-1);
  }

  data_at(panel_ind,row,col) = color;
}

/* 
   "data" will be used as if 2D array. Unlike Cash coordinate system,
   the size of data is [window_row]*[window_col]. Thus, the index rages
   from [0,0] to [window_row-1,window_col-1]. However,
   data_at(panel_ind,row,col) will get a value that ranges from [1,1] to
   [window_row,window_col]. Hence, we have to subtract 1 from row and
   col.
*/
unsigned char& CashDisplay::data_at(const unsigned panel_ind,const unsigned row,const unsigned col)
{
  return *(data + (panels[panel_ind].o_row + row - 1)*window_col + panels[panel_ind].o_col + col - 1);
}

void CashDisplay::draw_window()
{
  BlockDisplay(data,window_row,window_col,0,0,0);
}

void CashDisplay::draw_png()
{
  BlockPNG(data,window_row,window_col,0);
}

void CashDisplay::put_histogram(const unsigned panel_ind,const std::vector<double>& frequencies,const double max_y,const unsigned char color_of_blank,const unsigned char color_of_background,const unsigned char color_of_bar)
{
  std::vector<unsigned char> _color_of_bar(1,color_of_bar);
  put_histogram(panel_ind,frequencies,max_y,color_of_blank,color_of_background,_color_of_bar);
}

#endif
