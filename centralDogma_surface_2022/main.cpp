#include "molecule.h"
#include "newca.h"
#include "para.h"
#include "random.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string_view>

char **getCmd(char **begin, char **end, const std::string_view opt) {
  char **itr{std::find(begin, end, opt)};
  if (itr != end && ++itr != end) { // note that the last element of argv (i.e.
                                    // end - 1) is a nullptr
    return itr;
  }
  return nullptr;
}

void singleRun(std::string_view cType) {
  long init_time{0};
  newCA randomizedGrid(Para::sys_nrow, Para::sys_ncol);

  for (long t{init_time}; t <= 2000; t++) {
    /* if (t % Para::display_interval == 0) { */
    /* randomizedGrid.visualize(); */
    /* } */

    randomizedGrid.writeFile(t, "qp", &newCA::testComplex);
    /* randomizedGrid.visualize(); */
    randomizedGrid.update_squares();
    /* for (int i{0}; i < Para::grid_size; i++) */
    /*   randomizedGrid.reallycba(); */
    /* randomizedGrid.diffuse(row, col); */

    /* randomizedGrid.update_squares(); */
  }
}

int main(int argc, char *argv[]) {

  /* if (argc <= 2) { */
  /*   std::cout << "Enter valid arguments to proceed '\n'"; */
  /*   return 1; */
  /* } */

  /* char **compout{getCmd(argv, argv + argc, "-co")}; */
  /* auto getCompout{[](std::string_view coOpt) { */
  /*   return (coOpt == "pp" || coOpt == "pq" || coOpt == "qq" || coOpt ==
   * "qp"); */
  /* }}; */
  /* char **simpout{getCmd(argv, argv + argc, "-so")}; */

  /* if (compout) { */
  /*   for (int i{static_cast<int>(compout - &argv[0])}; i < argc; ++i) { */
  /*     std::cout << *std::find_if(&argv[i], argv + argc, getCompout); */
  /*   } */
  /* } */

  singleRun("pp");
  return 0;
}
