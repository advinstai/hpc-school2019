/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/cut/demo5,PairLJCutDemo5)

#else

#ifndef LMP_PAIR_LJ_CUT_DEMO5_H
#define LMP_PAIR_LJ_CUT_DEMO5_H

#include "pair_lj_cut.h"

namespace LAMMPS_NS {

class PairLJCutDemo5 : public PairLJCut {
 public:
  PairLJCutDemo5(class LAMMPS *);
  virtual ~PairLJCutDemo5() {};
  virtual void compute(int, int);

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
