/*
  CA2D is a very low level class for 2-dimensional Cellular
  Automata. It uses template to enable the use of user defined
  class as an automaton. 

  ------------------------------------------------------------
  Constructer:

  It has two parameters.

  nrow: the size of vertical dimension.
  ncol: the size of horizontal dimension.

  The argument to the template, <T>, is the type of the class
  that will be the type of the automata.

  When constructed, it allocates an <T>-array of size
  (nrow+2)*(ncol+2), where nrow and ncol are the grid
  size. During a simulation, the "inside" of the grid is
  supposed to be manipulated, i.e. the grid of size
  nrow*ncol. The other parts of the grid are considered as the
  boundaries of the grid.

  ------------------------------------------------------------
  Methods:

  It offers two methods to access the allocated automata in
  the array. There is no check for out-range access (CA2D is
  meant to be a very low level class).
  
  cell(row,col): 

  It returns a pointer to the automaton at coordinate
  (row,col). It also returns a pointer to the automaton in the
  boundaries if required.

  neigh_wrap(row,col,nei): 

  It returns a pointer to the automaton that is the nei-th
  neighbor of the automaton at (row,col). The boundaries are
  taken to be wrapped (toroidal).

  neigh_fix(row,col,nei): 


  It is the same as neigh_wrap(), except that the
  boundaries are taken to be fixed; i.e. it returns an
  automaton in the boundaries when required. When this method
  is used, the automata in the boundaries must be properly
  initilized by the user.

  xy_neigh_wrap(row,col,nei,nei_row,nei_col),
  xy_neigh_fix(row,col,nei,nei_row,nei_col):

  These methods allow to get the coordinate of the
  neighbor. The value of nei_row and nei_col will be set to
  the wanted coordinate.


*/

#include <iostream>
#include <unistd.h>

#ifndef CELLULARAUTOMATA
#define CELLULARAUTOMATA


template <class T> class CA2D {
 
private:
  unsigned nrow{};
  unsigned ncol{};
  unsigned ncol2{};
  
  /* "cells" is an array of Cells. Each cell contains an
     automaton. */
  T* cells;
  
public:
  inline CA2D(const unsigned a_nrow,const unsigned a_ncol);
  inline ~CA2D();
  inline T& cell(const unsigned row,const unsigned col);
  inline const T& cell(const unsigned row,const unsigned col) const;

  /* The following methods use an "index" to access a cell. An
     index is an offset to "cells". To use this functions, one
     has to know the index of the cell one wants to
     access. This can be done by index(). cell(index(row,col))
     is equivalent to cell(row,col).
     */
  inline T& cell(const unsigned ind);
  inline const T& cell(const unsigned ind) const;
  inline unsigned index (const unsigned row,const unsigned col) const;

  /* 1st-order neighbor retrieval */
  /* 516
     203
     748 */
  /* Wrapped boundaries (troidal) */
  inline T& neigh_wrap(const unsigned row,const unsigned col,const unsigned nei);
  inline const T& neigh_wrap(const unsigned row,const unsigned col,const unsigned nei) const;
  inline void xy_neigh_wrap(const unsigned row,const unsigned col,const unsigned nei,unsigned& nei_row,unsigned& nei_col) const;
  inline unsigned index_neigh_wrap(const unsigned row,const unsigned col,const unsigned nei) const;

  /* Fixed boundaries */
  inline T& neigh_fix(const unsigned row,const unsigned col,const unsigned nei);
  inline const T& neigh_fix(const unsigned row,const unsigned col,const unsigned nei) const;
  inline void xy_neigh_fix(const unsigned row,const unsigned col,const unsigned nei,unsigned& nei_row,unsigned& nei_col) const;
  inline unsigned index_neigh_fix(const unsigned row,const unsigned col,const unsigned nei) const;

  /* Miscellaneous functions */
  inline unsigned neighbor_converter(const unsigned nei);
  inline unsigned neigh_7_select(const unsigned nei,const unsigned nei_exclude);

  /* 2nd-order neighbor retrieval. */
  /*
                N
      [21][13][ 9][14][22]
      [15][ 5][ 1][ 6][16]
    W [10][ 2][ 0][ 3][11] E
      [17][ 7][ 4][ 8][18]
      [23][19][12][20][24]  
                S           
  */
  /* Wrapped boundaries (troidal) */
  inline T& neigh_wrap2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei);
  inline const T& neigh_wrap2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const;
  inline void xy_neigh_wrap2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei,unsigned& nei_row,unsigned& nei_col) const;
  inline unsigned index_neigh_wrap2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const;

  /* Overloaded methods for simplicity */
  T& neigh_wrap(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) {return neigh_wrap2(row,col,nei,neinei);}
  const T& neigh_wrap(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const {return neigh_wrap2(row,col,nei,neinei);}
  void xy_neigh_wrap(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei,unsigned& nei_row,unsigned& nei_col) const {xy_neigh_wrap2(row,col,nei,neinei,nei_row,nei_col);}
 
  /* Fixed boundaries */
  inline T& neigh_fix2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei);
  inline const T& neigh_fix2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const;
  inline void xy_neigh_fix2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei,unsigned& nei_row,unsigned& nei_col) const;
  inline unsigned index_neigh_fix2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const;

  /* Overloaded methods for simplicity */
  inline T& neigh_fix(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) {return neigh_fix2(row,col,nei,neinei);}
  inline const T& neigh_fix(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const {return neigh_fix2(row,col,nei,neinei);}
  inline void xy_neigh_fix(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei,unsigned& nei_row,unsigned& nei_col) const {xy_neigh_fix2(row,col,nei,neinei,nei_row,nei_col);}

  /* Miscellaneous functions */
  inline unsigned neighbor_converter(const unsigned nei,const unsigned neinei);
};

/* 
   Give to 'nei', a number between 1 to 7. Give to
   nei_exclude, the index of neighbor you want to exclude
   (because it has been chosen before).
*/
template <class T> unsigned CA2D<T>::neigh_7_select(const unsigned nei,const unsigned nei_exclude)
{
  return (nei >= nei_exclude)? (nei+1):nei;
}


/*
  When you want to know the nei-index refering to this cell
  from a view point of the neighbor of this cell whose
  nei-index is "nei".
*/
template <class T> unsigned CA2D<T>::neighbor_converter(const unsigned nei)
{
  switch(nei)
    {
    case 1: return 4;
    case 2: return 3;
    case 3: return 2;
    case 4: return 1;
    case 5: return 8;
    case 6: return 7;
    case 7: return 6;
    case 8: return 5;
    default: return 0;
    }
}

template <class T> unsigned CA2D<T>::index(const unsigned row,const unsigned col) const
{
  return row*ncol2 + col;
  //return row*512 + col;
}

template <class T> T& CA2D<T>::cell(const unsigned row,const unsigned col)
{
  return *(cells + index(row,col));
}

template <class T> const T& CA2D<T>::cell(const unsigned row,const unsigned col) const
{
  return *(cells + index(row,col));
}

template <class T> T& CA2D<T>::cell(const unsigned ind)
{
  return *(cells + ind);
}
template <class T> const T& CA2D<T>::cell(const unsigned ind) const
{
  return *(cells + ind);
}

template <class T> unsigned CA2D<T>::index_neigh_fix(const unsigned row,const unsigned col,const unsigned nei) const
{
    switch(nei)
    {
    case 0:/* C */
      return index(row,col);
    case 1:/* N */
      return index(row-1,col);
    case 2:/* W */
      return index(row,col-1);
    case 3:/* E */
      return index(row,col+1);
    case 4:/* S */
      return index(row+1,col);
    case 5:/* NW */
      return index(row-1,col-1);
    case 6:/* NE */
      return index(row-1,col+1);
    case 7:/* SW */
      return index(row+1,col-1);
    case 8:/* SE */
      return index(row+1,col+1);
    default:
      std::cerr << "CA2D::index_neigh_fix() Error, got an invalid nei value, which is " << nei << std::endl;
      exit(-1);
    }
}

template <class T> T& CA2D<T>::neigh_fix(const unsigned row,const unsigned col,const unsigned nei)
{
  return *(cells + index_neigh_fix(row,col,nei));
}

template <class T> const T& CA2D<T>::neigh_fix(const unsigned row,const unsigned col,const unsigned nei) const
{
  return *(cells + index_neigh_fix(row,col,nei));
}

template <class T> void CA2D<T>::xy_neigh_fix(const unsigned row,const unsigned col,const unsigned nei,unsigned& nei_row,unsigned& nei_col) const
{
  switch(nei)
    {
    case 0:/* C */
      nei_row=row;
      nei_col=col;
      return;
    case 1:/* N */
      nei_row=row-1;
      nei_col=col;
      return;
    case 2:/* W */
      nei_row=row;
      nei_col=col-1;
      return;
    case 3:/* E */
      nei_row=row;
      nei_col=col+1;
      return;
    case 4:/* S */
      nei_row=row+1;
      nei_col=col;
      return;
    case 5:/* NW */
      nei_row=row-1;
      nei_col=col-1;
      return;
    case 6:/* NE */
      nei_row=row-1;
      nei_col=col+1;
      return;
    case 7:/* SW */
      nei_row=row+1;
      nei_col=col-1;
      return;
    case 8:/* SE */
      nei_row=row+1;
      nei_col=col+1;
      return;
    default:
      std::cerr << "CA2D::xy_neigh_fix() Error, got an invalid nei value, which is " << nei << std::endl;
      exit(-1);
      return;
    }
}

template <class T> unsigned CA2D<T>::index_neigh_wrap(const unsigned row,const unsigned col,const unsigned nei) const
{
  switch(nei)
    {
    case 0:/* C */
      return index(row,col);
    case 1:/* N */
      return index(((row!=1)? (row-1):nrow),col);
    case 2:/* W */
      return index(row,((col!=1)? (col-1):ncol));
    case 3:/* E */
      return index(row,((col!=ncol)? (col+1):1));
    case 4:/* S */
      return index(((row!=nrow)? (row+1):1),col);
    case 5:/* NW */
      return index(((row!=1)? (row-1):nrow),((col!=1)? (col-1):ncol));
    case 6:/* NE */
      return index(((row!=1)? (row-1):nrow),((col!=ncol)? (col+1):1));
    case 7:/* SW */
      return index(((row!=nrow)? (row+1):1),((col!=1)? (col-1):ncol));
    case 8:/* SE */
      return index(((row!=nrow)? (row+1):1),((col!=ncol)? (col+1):1));
    default:
      std::cerr << "CA2D::index_neigh_wrap() Error, got an invalid nei value, which is " << nei << std::endl;
      exit(-1);
    }
}

template <class T> T& CA2D<T>::neigh_wrap(const unsigned row,const unsigned col,const unsigned nei)
{
  return *(cells + index_neigh_wrap(row,col,nei));
}

template <class T> const T& CA2D<T>::neigh_wrap(const unsigned row,const unsigned col,const unsigned nei) const
{
  return *(cells + index_neigh_wrap(row,col,nei));
}

template <class T> void CA2D<T>::xy_neigh_wrap(const unsigned row,const unsigned col,const unsigned nei,unsigned& nei_row,unsigned& nei_col) const
{
  switch(nei)
    {
    case 0:/* C */
      nei_row=row;
      nei_col=col;
      return;
    case 1:/* N */
      nei_row=((row!=1)? (row-1):nrow);
      nei_col=col;
      return;
    case 2:/* W */
      nei_row=row;
      nei_col=((col!=1)? (col-1):ncol);
      return;
    case 3:/* E */
      nei_row=row;
      nei_col=((col!=ncol)? (col+1):1);
      return;
    case 4:/* S */
      nei_row=((row!=nrow)? (row+1):1);
      nei_col=col;
      return;
    case 5:/* NW */
      nei_row=((row!=1)? (row-1):nrow);
      nei_col=((col!=1)? (col-1):ncol);
      return;
    case 6:/* NE */
      nei_row=((row!=1)? (row-1):nrow);
      nei_col=((col!=ncol)? (col+1):1);
      return;
    case 7:/* SW */
      nei_row=((row!=nrow)? (row+1):1);
      nei_col=((col!=1)? (col-1):ncol);
      return;
    case 8:/* SE */
      nei_row=((row!=nrow)? (row+1):1);
      nei_col=((col!=ncol)? (col+1):1);
      return;
    default:
      std::cerr << "CA2D::xy_neigh_wrap() Error, got an invalid nei value, which is " << nei << std::endl;
      return;
    }
}

template <class T> CA2D<T>::CA2D(const unsigned a_nrow,const  unsigned a_ncol)
  : nrow(a_nrow), 
    ncol(a_ncol),
    ncol2(a_ncol+2)
{
  if(nrow==0 || ncol==0){
    std::cerr << "CA2D() Error: nrow=ncol=0 is not allowed." << std::endl;
    return;
  }

  try{
    cells = new T [(nrow+2)*(ncol+2)];
  }catch(std::bad_alloc){
    std::cerr << "CA2D::CA2D() Memory exausted" << std::endl;
    exit(-1);
  }
}

template <class T>  CA2D<T>::~CA2D()
{
  delete[] cells;
}

template <class T> unsigned CA2D<T>::index_neigh_fix2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const
{
  /* We expect that nei=[0,8] and neinei=[0,8], so they take 3 bits.
     Thus, we left-shift nei 4 bits and add it to neinei by |
     operator. This gives a unique number for each pair (nei,neinei). We
     then assign the neighbor index to these unique numbers, which can
     be done by another program. The following is the output of that
     program. */
  switch((nei << 4) | neinei)
    {
    case 0: case 20: case 35: case 50: case 65: case 88: case 103: case 118: case 133: goto nei0;
    case 1: case 16: case 38: case 53: case 83: case 98: goto nei1;
    case 2: case 23: case 32: case 69: case 84: case 113: goto nei2;
    case 3: case 24: case 48: case 70: case 100: case 129: goto nei3;
    case 4: case 40: case 55: case 64: case 115: case 130: goto nei4;
    case 5: case 18: case 33: case 80: goto nei5;
    case 6: case 19: case 49: case 96: goto nei6;
    case 7: case 36: case 66: case 112: goto nei7;
    case 8: case 52: case 67: case 128: goto nei8;
    case 17: case 86: case 101: goto nei9;
    case 34: case 87: case 117: goto nei10;
    case 51: case 104: case 134: goto nei11;
    case 68: case 120: case 135: goto nei12;
    case 21: case 81: goto nei13;
    case 22: case 97: goto nei14;
    case 37: case 82: goto nei15;
    case 54: case 99: goto nei16;
    case 39: case 114: goto nei17;
    case 56: case 131: goto nei18;
    case 71: case 116: goto nei19;
    case 72: case 132: goto nei20;
    case 85: goto nei21;
    case 102: goto nei22;
    case 119: goto nei23;
    case 136: goto nei24;
    default:
      std::cerr << "CA2D::index_neigh_fix2() Error, got an invalid argument value(s), nei= " << nei << " neinei="<< neinei << std::endl;
      exit(-1);
    }
  
 nei0:/* C */
  return index(                            row,                             col);
 nei1:/* N */
  return index(                          row-1,                             col);
 nei2:/* W */
  return index(                            row,                           col-1);
 nei3:/* E */
  return index(                            row,                           col+1);
 nei4:/* S */
  return index(                          row+1,                             col);
 nei5:/* NW */
  return index(                          row-1,                           col-1);
 nei6:/* NE */
  return index(                          row-1,                           col+1);
 nei7:/* SW */
  return index(                          row+1,                           col-1);
 nei8:/* SE */
  return index(                          row+1,                           col+1);
 nei9:/* Nn */
  return index(            ((row>2)?(row-2):0),                             col);
 nei10:/* Ww */
  return index(                            row,             ((col>2)?(col-2):0));
 nei11:/* Ee */
  return index(                            row,  ((col!=ncol)?(col+2):(ncol+1)));
 nei12:/* Ss */
  return index( ((row!=nrow)?(row+2):(nrow+1)),                             col);
 nei13:/* NnW */
  return index(            ((row>2)?(row-2):0),                           col-1);
 nei14:/* NnE */
  return index(            ((row>2)?(row-2):0),                           col+1);
 nei15:/* NWw */
  return index(                          row-1,             ((col>2)?(col-2):0));
 nei16:/* NEe */
  return index(                          row-1,  ((col!=ncol)?(col+2):(ncol+1)));
 nei17:/* SWw */
  return index(                          row+1,             ((col>2)?(col-2):0));
 nei18:/* SEe */
  return index(                          row+1,  ((col!=ncol)?(col+2):(ncol+1)));
 nei19:/* Ssw */
  return index( ((row!=nrow)?(row+2):(nrow+1)),                           col-1);
 nei20:/* SsE */
  return index( ((row!=nrow)?(row+2):(nrow+1)),                           col+1);
 nei21:/* NnWw */
  return index(            ((row>2)?(row-2):0),             ((col>2)?(col-2):0));
 nei22:/* NnEe */
  return index(            ((row>2)?(row-2):0),  ((col!=ncol)?(col+2):(ncol+1)));
 nei23:/* SsWw */
  return index( ((row!=nrow)?(row+2):(nrow+1)),             ((col>2)?(col-2):0));
 nei24:/* SsEe */
  return index( ((row!=nrow)?(row+2):(nrow+1)),  ((col!=ncol)?(col+2):(ncol+1)));
}

template <class T> T& CA2D<T>::neigh_fix2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei)
{
  return *(cells + index_neigh_fix2(row,col,nei,neinei));
}

template <class T> const T& CA2D<T>::neigh_fix2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const
{
  return *(cells + index_neigh_fix2(row,col,nei,neinei));
}

template <class T> void CA2D<T>::xy_neigh_fix2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei,unsigned& nei_row,unsigned& nei_col) const
{
  switch((nei << 4) | neinei)
    {
    case 0: case 20: case 35: case 50: case 65: case 88: case 103: case 118: case 133: goto nei0;
    case 1: case 16: case 38: case 53: case 83: case 98: goto nei1;
    case 2: case 23: case 32: case 69: case 84: case 113: goto nei2;
    case 3: case 24: case 48: case 70: case 100: case 129: goto nei3;
    case 4: case 40: case 55: case 64: case 115: case 130: goto nei4;
    case 5: case 18: case 33: case 80: goto nei5;
    case 6: case 19: case 49: case 96: goto nei6;
    case 7: case 36: case 66: case 112: goto nei7;
    case 8: case 52: case 67: case 128: goto nei8;
    case 17: case 86: case 101: goto nei9;
    case 34: case 87: case 117: goto nei10;
    case 51: case 104: case 134: goto nei11;
    case 68: case 120: case 135: goto nei12;
    case 21: case 81: goto nei13;
    case 22: case 97: goto nei14;
    case 37: case 82: goto nei15;
    case 54: case 99: goto nei16;
    case 39: case 114: goto nei17;
    case 56: case 131: goto nei18;
    case 71: case 116: goto nei19;
    case 72: case 132: goto nei20;
    case 85: goto nei21;
    case 102: goto nei22;
    case 119: goto nei23;
    case 136: goto nei24;
    default:
      std::cerr << "CA2D::xy_neigh_fix2() Error, got an invalid argument value(s), nei= " << nei << " neinei="<< neinei << std::endl;
      exit(-1);
    }

 nei0:/* C */
  nei_row=row;
  nei_col=col;
  return;
 nei1:/* N */
  nei_row=row-1;
  nei_col=col;
  return;
 nei2:/* W */
  nei_row=row;
  nei_col=col-1;
  return;
 nei3:/* E */
  nei_row=row;
  nei_col=col+1;
  return;
 nei4:/* S */
  nei_row=row+1;
  nei_col=col;
  return;
 nei5:/* NW */
  nei_row=row-1;
  nei_col=col-1;
  return;
 nei6:/* NE */
  nei_row=row-1;
  nei_col=col+1;
  return;
 nei7:/* SW */
  nei_row=row+1;
  nei_col=col-1;
  return;
 nei8:/* SE */
  nei_row=row+1;
  nei_col=col+1;
  return;
 nei9:/* Nn */
  nei_row=((row>2)?(row-2):0);
  nei_col=col;
  return;
 nei10:/* Ww */
  nei_row=row;
  nei_col=((col>2)?(col-2):0);
  return;
 nei11:/* Ee */
  nei_row=row;
  nei_col=((col!=ncol)?(col+2):(ncol+1));
  return;
 nei12:/* Ss */
  nei_row=((row!=nrow)?(row+2):(nrow+1));
  nei_col=col;
  return;
 nei13:/* NnW */
  nei_row=((row>2)?(row-2):0);
  nei_col=col-1;
  return;
 nei14:/* NnE */
  nei_row=((row>2)?(row-2):0);
  nei_col=col+1;
  return;
 nei15:/* NWw */
  nei_row=row-1;
  nei_col=((col>2)?(col-2):0);
  return;
 nei16:/* NEe */
  nei_row=row-1;
  nei_col=((col!=ncol)?(col+2):(ncol+1));
  return;
 nei17:/* SWw */
  nei_row=row+1;
  nei_col=((col>2)?(col-2):0);
  return;
 nei18:/* SEe */
  nei_row=row+1;
  nei_col=((col!=ncol)?(col+2):(ncol+1));
  return;
 nei19:/* Ssw */
  nei_row=((row!=nrow)?(row+2):(nrow+1));
  nei_col=col-1;
  return;
 nei20:/* SsE */
  nei_row=((row!=nrow)?(row+2):(nrow+1));
  nei_col=col+1;
  return;
 nei21:/* NnWw */
  nei_row=((row>2)?(row-2):0);
  nei_col=((col>2)?(col-2):0);
  return;
 nei22:/* NnEe */
  nei_row=((row>2)?(row-2):0);
  nei_col=((col!=ncol)?(col+2):(ncol+1));
  return;
 nei23:/* SsWw */
  nei_row=((row!=nrow)?(row+2):(nrow+1));
  nei_col=((col>2)?(col-2):0);
  return;
 nei24:/* SsEe */
  nei_row=((row!=nrow)?(row+2):(nrow+1));
  nei_col=((col!=ncol)?(col+2):(ncol+1));
  return;
}

template <class T> unsigned CA2D<T>::index_neigh_wrap2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const
{
  switch((nei << 4) | neinei)
    {
    case 0: case 20: case 35: case 50: case 65: case 88: case 103: case 118: case 133: goto nei0;
    case 1: case 16: case 38: case 53: case 83: case 98: goto nei1;
    case 2: case 23: case 32: case 69: case 84: case 113: goto nei2;
    case 3: case 24: case 48: case 70: case 100: case 129: goto nei3;
    case 4: case 40: case 55: case 64: case 115: case 130: goto nei4;
    case 5: case 18: case 33: case 80: goto nei5;
    case 6: case 19: case 49: case 96: goto nei6;
    case 7: case 36: case 66: case 112: goto nei7;
    case 8: case 52: case 67: case 128: goto nei8;
    case 17: case 86: case 101: goto nei9;
    case 34: case 87: case 117: goto nei10;
    case 51: case 104: case 134: goto nei11;
    case 68: case 120: case 135: goto nei12;
    case 21: case 81: goto nei13;
    case 22: case 97: goto nei14;
    case 37: case 82: goto nei15;
    case 54: case 99: goto nei16;
    case 39: case 114: goto nei17;
    case 56: case 131: goto nei18;
    case 71: case 116: goto nei19;
    case 72: case 132: goto nei20;
    case 85: goto nei21;
    case 102: goto nei22;
    case 119: goto nei23;
    case 136: goto nei24;
    default:
      std::cerr << "CA2D::index_neigh_wrap2() Error, got an invalid argument value(s), nei= " << nei << " neinei="<< neinei << std::endl;
      exit(-1);
    }
  
 nei0:/* C */
  return index(                                       row,                                        col);
 nei1:/* N */
  return index(                  ((row!=1)? (row-1):nrow),                                        col);
 nei2:/* W */
  return index(                                       row,                   ((col!=1)? (col-1):ncol));
 nei3:/* E */
  return index(                                       row,                   ((col!=ncol)? (col+1):1));
 nei4:/* S */
  return index(                  ((row!=nrow)? (row+1):1),                                        col);
 nei5:/* NW */
  return index(                  ((row!=1)? (row-1):nrow),                   ((col!=1)? (col-1):ncol));
 nei6:/* NE */
  return index(                  ((row!=1)? (row-1):nrow),                   ((col!=ncol)? (col+1):1));
 nei7:/* SW */
  return index(                  ((row!=nrow)? (row+1):1),                   ((col!=1)? (col-1):ncol));
 nei8:/* SE */
  return index(                  ((row!=nrow)? (row+1):1),                   ((col!=ncol)? (col+1):1));
 nei9:/* Nn */
  return index(((row>2)?(row-2):((row==1)?(nrow-1):nrow)),                                        col);
 nei10:/* Ww */
  return index(                                       row, ((col>2)?(col-2):((col==1)?(ncol-1):ncol)));
 nei11:/* Ee */
  return index(                                       row,   ((col<ncol-1)?(col+2):((col==ncol)?2:1)));
 nei12:/* Ss */
  return index(  ((row<nrow-1)?(row+2):((row==nrow)?2:1)),                                        col);
 nei13:/* NnW */
  return index(((row>2)?(row-2):((row==1)?(nrow-1):nrow)),                   ((col!=1)? (col-1):ncol));
 nei14:/* NnE */
  return index(((row>2)?(row-2):((row==1)?(nrow-1):nrow)),                   ((col!=ncol)? (col+1):1));
 nei15:/* NWw */
  return index(                  ((row!=1)? (row-1):nrow), ((col>2)?(col-2):((col==1)?(ncol-1):ncol)));
 nei16:/* NEe */
  return index(                  ((row!=1)? (row-1):nrow),   ((col<ncol-1)?(col+2):((col==ncol)?2:1)));
 nei17:/* SWw */
  return index(                  ((row!=nrow)? (row+1):1), ((col>2)?(col-2):((col==1)?(ncol-1):ncol)));
 nei18:/* SEe */
  return index(                  ((row!=nrow)? (row+1):1),   ((col<ncol-1)?(col+2):((col==ncol)?2:1)));
 nei19:/* Ssw */
  return index(  ((row<nrow-1)?(row+2):((row==nrow)?2:1)),                   ((col!=1)? (col-1):ncol));
 nei20:/* SsE */
  return index(  ((row<nrow-1)?(row+2):((row==nrow)?2:1)),                   ((col!=ncol)? (col+1):1));
 nei21:/* NnWw */
  return index(((row>2)?(row-2):((row==1)?(nrow-1):nrow)), ((col>2)?(col-2):((col==1)?(ncol-1):ncol)));
 nei22:/* NnEe */
  return index(((row>2)?(row-2):((row==1)?(nrow-1):nrow)),   ((col<ncol-1)?(col+2):((col==ncol)?2:1)));
 nei23:/* SsWw */
  return index(  ((row<nrow-1)?(row+2):((row==nrow)?2:1)), ((col>2)?(col-2):((col==1)?(ncol-1):ncol)));
 nei24:/* SsEe */
  return index(  ((row<nrow-1)?(row+2):((row==nrow)?2:1)),   ((col<ncol-1)?(col+2):((col==ncol)?2:1)));
}

template <class T> T& CA2D<T>::neigh_wrap2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei)
{
  return *(cells + index_neigh_wrap2(row,col,nei,neinei));
}

template <class T> const T& CA2D<T>::neigh_wrap2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei) const
{
  return *(cells + index_neigh_wrap2(row,col,nei,neinei));
}

template <class T> void CA2D<T>::xy_neigh_wrap2(const unsigned row,const unsigned col,const unsigned nei,const unsigned neinei,unsigned& nei_row,unsigned& nei_col) const
{
  switch((nei << 4) | neinei)
    {
    case 0: case 20: case 35: case 50: case 65: case 88: case 103: case 118: case 133: goto nei0;
    case 1: case 16: case 38: case 53: case 83: case 98: goto nei1;
    case 2: case 23: case 32: case 69: case 84: case 113: goto nei2;
    case 3: case 24: case 48: case 70: case 100: case 129: goto nei3;
    case 4: case 40: case 55: case 64: case 115: case 130: goto nei4;
    case 5: case 18: case 33: case 80: goto nei5;
    case 6: case 19: case 49: case 96: goto nei6;
    case 7: case 36: case 66: case 112: goto nei7;
    case 8: case 52: case 67: case 128: goto nei8;
    case 17: case 86: case 101: goto nei9;
    case 34: case 87: case 117: goto nei10;
    case 51: case 104: case 134: goto nei11;
    case 68: case 120: case 135: goto nei12;
    case 21: case 81: goto nei13;
    case 22: case 97: goto nei14;
    case 37: case 82: goto nei15;
    case 54: case 99: goto nei16;
    case 39: case 114: goto nei17;
    case 56: case 131: goto nei18;
    case 71: case 116: goto nei19;
    case 72: case 132: goto nei20;
    case 85: goto nei21;
    case 102: goto nei22;
    case 119: goto nei23;
    case 136: goto nei24;
    default:
      std::cerr << "CA2D::xy_neigh_wrap2() Error, got an invalid argument value(s), nei= " << nei << " neinei="<< neinei << std::endl;
      exit(-1);
    }
  
 nei0:/* C */
  nei_row=row;
  nei_col=col;
  return;
 nei1:/* N */
  nei_row=((row!=1)? (row-1):nrow);
  nei_col=col;
  return;
 nei2:/* W */
  nei_row=row;
  nei_col=((col!=1)? (col-1):ncol);
  return;
 nei3:/* E */
  nei_row=row;
  nei_col=((col!=ncol)? (col+1):1);
  return;
 nei4:/* S */
  nei_row=((row!=nrow)? (row+1):1);
  nei_col=col;
  return;
 nei5:/* NW */
  nei_row=((row!=1)? (row-1):nrow);
  nei_col=((col!=1)? (col-1):ncol);
  return;
 nei6:/* NE */
  nei_row=((row!=1)? (row-1):nrow);
  nei_col=((col!=ncol)? (col+1):1);
  return;
 nei7:/* SW */
  nei_row=((row!=nrow)? (row+1):1);
  nei_col=((col!=1)? (col-1):ncol);
  return;
 nei8:/* SE */
  nei_row=((row!=nrow)? (row+1):1);
  nei_col=((col!=ncol)? (col+1):1);
  return;
 nei9:/* Nn */
  nei_row=((row>2)?(row-2):((row==1)?(nrow-1):nrow));
  nei_col=col;
  return;
 nei10:/* Ww */
  nei_row=row;
  nei_col=((col>2)?(col-2):((col==1)?(ncol-1):ncol));
  return;
 nei11:/* Ee */
  nei_row=row;
  nei_col=((col<ncol-1)?(col+2):((col==ncol)?2:1));
  return;
 nei12:/* Ss */
  nei_row=((row<nrow-1)?(row+2):((row==nrow)?2:1));
  nei_col=col;
  return;
 nei13:/* NnW */
  nei_row=((row>2)?(row-2):((row==1)?(nrow-1):nrow));
  nei_col=((col!=1)? (col-1):ncol);
  return;
 nei14:/* NnE */
  nei_row=((row>2)?(row-2):((row==1)?(nrow-1):nrow));
  nei_col=((col!=ncol)? (col+1):1);
  return;
 nei15:/* NWw */
  nei_row=((row!=1)? (row-1):nrow);
  nei_col=((col>2)?(col-2):((col==1)?(ncol-1):ncol));
  return;
 nei16:/* NEe */
  nei_row=((row!=1)? (row-1):nrow);
  nei_col=((col<ncol-1)?(col+2):((col==ncol)?2:1));
  return;
 nei17:/* SWw */
  nei_row=((row!=nrow)? (row+1):1);
  nei_col=((col>2)?(col-2):((col==1)?(ncol-1):ncol));
  return;
 nei18:/* SEe */
  nei_row=((row!=nrow)? (row+1):1);
  nei_col=((col<ncol-1)?(col+2):((col==ncol)?2:1));
  return;
 nei19:/* Ssw */
  nei_row=((row<nrow-1)?(row+2):((row==nrow)?2:1));
  nei_col=((col!=1)? (col-1):ncol);
  return;
 nei20:/* SsE */
  nei_row=((row<nrow-1)?(row+2):((row==nrow)?2:1));
  nei_col=((col!=ncol)? (col+1):1);
  return;
 nei21:/* NnWw */
  nei_row=((row>2)?(row-2):((row==1)?(nrow-1):nrow));
  nei_col=((col>2)?(col-2):((col==1)?(ncol-1):ncol));
  return;
 nei22:/* NnEe */
  nei_row=((row>2)?(row-2):((row==1)?(nrow-1):nrow));
  nei_col=((col<ncol-1)?(col+2):((col==ncol)?2:1));
  return;
 nei23:/* SsWw */
  nei_row=((row<nrow-1)?(row+2):((row==nrow)?2:1));
  nei_col=((col>2)?(col-2):((col==1)?(ncol-1):ncol));
  return;
 nei24:/* SsEe */
  nei_row=((row<nrow-1)?(row+2):((row==nrow)?2:1));
  nei_col=((col<ncol-1)?(col+2):((col==ncol)?2:1));
  return;
}

template <class T> unsigned CA2D<T>::neighbor_converter(const unsigned nei,const unsigned neinei)
{
  switch((nei << 4) | neinei)
    {
    case 0: case 20: case 35: case 50: case 65: case 88: case 103: case 118: case 133: return 0;
    case 1: case 16: case 38: case 53: case 83: case 98: return 1;
    case 2: case 23: case 32: case 69: case 84: case 113: return 2;
    case 3: case 24: case 48: case 70: case 100: case 129: return 3;
    case 4: case 40: case 55: case 64: case 115: case 130: return 4;
    case 5: case 18: case 33: case 80: return 5;
    case 6: case 19: case 49: case 96: return 6;
    case 7: case 36: case 66: case 112: return 7;
    case 8: case 52: case 67: case 128: return 8;
    case 17: case 86: case 101: return 9;
    case 34: case 87: case 117: return 10;
    case 51: case 104: case 134: return 11;
    case 68: case 120: case 135: return 12;
    case 21: case 81: return 13;
    case 22: case 97: return 14;
    case 37: case 82: return 15;
    case 54: case 99: return 16;
    case 39: case 114: return 17;
    case 56: case 131: return 18;
    case 71: case 116: return 19;
    case 72: case 132: return 20;
    case 85: return 21;
    case 102: return 22;
    case 119: return 23;
    case 136: return 24;
    default:
      std::cerr << "CA2D::neighbor_converter(nei,neinei) Error, got an invalid argument value(s), nei= " << nei << " neinei="<< neinei << std::endl;
      exit(-1);
    }
}

#endif

