#include "newca.h"
#include "cash-display.hpp"
#include "cash.h"
#include "molecule.h"

newCA::newCA(const unsigned a_nrow, const unsigned a_ncol)
    : nrow{a_nrow}, ncol{a_ncol}, plane(a_nrow, a_ncol) {
  for (unsigned row = 0; row < nrow; row++) {
    for (unsigned col = 0; col < ncol; col++) {

      plane.cell(row, col).setTypeRep(
