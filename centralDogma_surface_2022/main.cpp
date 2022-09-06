#include "molecule.h"
#include "newca.h"
#include "para.h"
#include "random.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <string_view>

// Returns the iterator pointing to opt in argv[], where opt is the shorthand of
// the option (e.g. -v for visualization), nullptr if opt is not found. Begin &
// end should be begin() and end() of argv[] (see function call).
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

// Sets the parameter that corresponds to the myOpt iterator of argv[] using
// std::map; takes an optType parameter to make the correct string-to-number
// conversion and emit the correct warning messege.
void useCmd(char **myOpt,
            const char optType) { // optType: i = int, d = double
  const static std::map<std::string_view, unsigned *> myIntCmd{
      {"-C", &Para::sys_ncol},      {"-R", &Para::sys_nrow},
      {"-v", &Para::visualization}, {"-d", &Para::display_interval},
      {"-s", &Para::mt_seed},       {"-t", &Para::max_time},
      {"-m", &Para::movie},         {"-r", &Para::read}};

  const static std::map<std::string_view, double *> myDoubleCmd{
      {"-D", &Para::decay_probability},
      {"-I", &Para::diffusion_probability},
      {"-M", &Para::mutation_probability},
      {"-A", &Para::alpha},
      {"-B", &Para::beta},
      {"-G", &Para::mutation_interval}};

  using namespace std::literals;
  if (optType == 'i') {
    auto optInd{myIntCmd.find(static_cast<std::string_view>(*myOpt))};
    /* std::cout << optInd->first; */
    /* std::cout << *(myOpt + 1); */
    if (*(myOpt + 1) == "t"sv) {
      *optInd->second = 1;
      /* std::cout << ' ' << *optInd->second << Para::visualization; */
    } else if (*(myOpt + 1) ==
               "f"sv) { // if arg for bool option is neither t or f
      *optInd->second = 0;
    } else {
      try {
        *optInd->second = static_cast<unsigned>(std::stoi(*(myOpt + 1)));
        if (static_cast<std::string_view>(*myOpt) == "-C" ||
            static_cast<std::string_view>(*myOpt) == "-R") {
          Para::grid_size = Para::sys_nrow * Para::sys_ncol;
          DiceRoller::randomRowOrCol = std::uniform_int_distribution<>{
              1, static_cast<int>(Para::sys_nrow)};
        }
        /* if (Para::mt_seed) */
        /*   DiceRoller::twister.seed(Para::mt_seed); */
        /* int seeds[8]; */
        /* DiceRoller::ss.param(seeds); */
        /* for (int i{0}; i <= 7; i++) */
        /*   std::cout << seeds[i] << '\n'; */
      } catch (std::invalid_argument) {
        std::cerr << "Usage:\n    " << (*myOpt)[1]
                  << " takes a positive integer argument (negatives/fractions "
                     "cause undefined behaviour and no warnings will be "
                     "emitted). "
                     "Alternatively, use "
                     "'t'/'f' for 1 and 0 (applicable to -v and -m)\n";
        exit(1);
      }
    }
  } else if (optType == 'd') {
    auto optDouble{myDoubleCmd.find(static_cast<std::string_view>(*myOpt))};
    try {
      *optDouble->second = std::stod(*(myOpt + 1));
    } catch (std::invalid_argument) {
      std::cerr
          << "Usage:\n    " << (*myOpt)[1]
          << " takes a positive integer/fraction argument (negatives cause "
             "undefined behaviour and no warnings will be emitted)\n";
      exit(1);
    }
  }
}

void singleRun() {
  long init_time{0};
  newCA randomizedGrid(Para::sys_nrow, Para::sys_ncol);
  /* randomizedGrid.testDummy(1); */

  for (long t{init_time}; t <= Para::max_time; t++) {
    if (Para::visualization == 1) {
      if (t % Para::display_interval == 0)
        randomizedGrid.visualize();
    }

    /* std::cout << t << '\n'; */

    if (t % 10000 == 0) {
      randomizedGrid.writeField(t);
    }

    if (t % 1000 == 0) {
      randomizedGrid.writeDensity(t, &newCA::testComplex);
      randomizedGrid.writeDistribution(t);
      for (int i{0}; i <= 7; i++)
        randomizedGrid.writeAverageK(i, t);
    }

    randomizedGrid.update_squares();
  }
  /* randomizedGrid.testDummy(1); */
}

int main(int argc, char *argv[]) {

  if (argc <= 2) {
    std::cout << "Enter valid arguments to proceed; for help, use \"-h 1\"\n";
    return 1;
  }

  if (argc > 2) {
    std::string_view cmdList{"-h-v-m-s-d-t-r-C-R-D-I-M-A-B-G"};
    for (int i{0}; i < cmdList.length(); i++) {
      if (cmdList[i] == '-') {
        char **cmdItr = getCmd(argv, argv + argc, cmdList.substr(i, 2));
        /* std::cout << *cmdItr; */
        if (cmdItr && i == 0) {
          std::cout << "CLI options are as follows:"
                    << "\n-h show this messege"
                    << "\n-C ncol (CA)"
                    << "\n-R nrow (CA)"
                    << "\n-v visualization (png)"
                    << "\n-d display_interval"
                    << "\n-s seed"
                    << "\n-t max_time"
                    << "\n-m movie (view in real time)"
                    << "\n-D decay_probability"
                    << "\n-I diffusion_probability"
                    << "\n-M mutation_probability"
                    << "\n-A alpha"
                    << "\n-B beta"
                    << "\n-G mutation_interval";
          exit(0);
        } else if (cmdItr && i <= 16) // 12 is the index of -R, which is the
                                      // last int parameter
          useCmd(cmdItr, 'i');
        else if (cmdItr && i <= 28) // 24 is theindex of -G, which is the last
                                    // double parameter
          useCmd(cmdItr, 'd');
      }
    }

    if (Para::sys_nrow != Para::sys_ncol) {
      /* std::cout << Para::sys_nrow << ' ' << Para::sys_ncol << '\n'; */
      std::cerr << "CA can only be square for now. Enter the same "
                   "number for rows and columns of CA.";
      exit(1);
    }

    std::cout << "mutation probability = " << Para::mutation_probability
              << '\n';
    std::cout << "mutation interval = " << Para::mutation_interval << '\n';
    std::cout << Para::sys_nrow << 'x' << Para::sys_ncol << '\n';
    singleRun();
    return 0;
  }
}
