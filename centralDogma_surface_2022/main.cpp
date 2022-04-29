#include "molecule.h"
#include "newca.h"
#include "para.h"
#include "random.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string_view>

char **getCmd(char **begin, char **end, const std::string_view opt) {
  char **itr{std::find(begin, end, opt)};
  if (itr != end &&
      ++itr != end) { // note that the last element of argv (i.e.
                      // end - 1) is a nullptr, therefore we need itr to not be
                      // the last element, i.e. ++itr != end
    --itr;
    return itr;
  }
  return nullptr;
}

void useBoolCmd(char **myOpt) {
  std::map<std::string_view, unsigned *> myBoolCmd{{"-v", &Para::visualization},
                                                   {"-m", &Para::movie}};
  using namespace std::literals;
  auto optInd{myBoolCmd.find(static_cast<std::string_view>(*myOpt))};
  std::cout << optInd->first;
  std::cout << *(myOpt + 1);
  if (*(myOpt + 1) == "t"sv) {
    std::cout << "true";
    *optInd->second = 1;
    std::cout << ' ' << *optInd->second << Para::visualization;
  } else if (*(myOpt + 1) != "f"sv) // if arg for bool option is neither t or f
    std::cout << "Usage: bool (-v & -m) take t or f as arguments\n";
}

void singleRun(std::string_view cType) {
  long init_time{0};
  newCA randomizedGrid(Para::sys_nrow, Para::sys_ncol);

  for (long t{init_time}; t <= 2000; t++) {
    /* if (t % Para::display_interval == 0) { */
    /* randomizedGrid.myOptize(); */
    /* } */
    randomizedGrid.writeFile(t, &newCA::testComplex);
    randomizedGrid.update_squares();
  }
}

int main(int argc, char *argv[]) {

  if (argc <= 2) {
    std::cout << "Enter valid arguments to proceed \n";
    return 1;
  }

  if (argc > 2) {
    std::string_view cmdList{"-v-m-d-t-S-De-Di-M-A-B-G-D"};
    for (int i{0}; i < cmdList.length(); i++) {
      if (i % 2 == 0) {
        /* std::cout << cmd << ' '; */
        char **cmdItr{getCmd(argv, argv + argc, cmdList.substr(i, 2))};
        std::cout << *cmdItr;
        if (cmdItr && i <= 2) // so far, only the frist 2 of cmdList are bool
          useBoolCmd(cmdItr);
      }
      /* singleRun("pp"); */
    }
  }
  return 0;
}
