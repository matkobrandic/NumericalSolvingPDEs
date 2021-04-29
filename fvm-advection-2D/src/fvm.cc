#include "config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/uggrid.hh>

// Za čitanje UUGrida iz DGF datoteke trebamo dvije datoteke
#include <dune/grid/io/file/dgfparser/gridptr.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>

// Za čitanje .input datoteke
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>

// VTK
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <stdexcept>
#include "data.hh"
#include "initialize.hh"
#include "evolve.hh"

// GLANI PROGRAM
int main (int argc , char ** argv)
{
    // Vaš kod dolazi ovdje.
    return 0;
}
