#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

#include "bc_extension.hh"
#include "driver.hh"




int main(int argc, char **argv) {

  Dune::MPIHelper::instance(argc, argv);

  // ovdje zadajemo sve parametre
  typedef double RF; //range field type
  RF mu             = 1;
  RF rho            = 1;
  RF dt             = 0.001;
  RF T_time         = 0.2;

  // Učitavanje mreže
  const int dim = 2;
  typedef Dune::UGGrid<dim> GridType;
  std::unique_ptr<GridType> pgrid{Dune::GmshReader<GridType>::read("prepreka1.msh", true, false)};
  //pgrid->globalRefine(1); //ako želimo profiniti mrežu

  using GV =  GridType::LeafGridView;
  const GV &gv = pgrid->leafGridView();

  // konstrukcija klase paramatara, rubni i inicijalni uvjeti dani u bc_extension
  using BdryPressure = ZeroFunction<GV, RF, 1>;
  using BdryVelocity = Velocity<GV, RF, 2>;
  using BdrySolution =  Dune::PDELab::CompositeGridFunction<BdryVelocity, BdryPressure>;

  BdryVelocity bdry_velocity(gv);
  BdryPressure bdry_pressure(gv);
  BdrySolution bdry_solution(bdry_velocity, bdry_pressure);

  BCTypeParam bdry_type;

  using SourceFunction = ZeroFunction<GV, RF, dim>; //općenito: f, u nasem slucaju f=0
  using NeumannFlux    = ZeroFunction<GV, RF, dim>;

  NeumannFlux    neumann_flux(gv);
  SourceFunction source_function(gv);

  const bool navier = true; // odreduje imamo li N-S ili S
  const bool tensor = false; // Treba li raditi sa simetriziranim gradijentom ili ne
                            // (utječe na interpretaciju rubnih uvjeta)

  using Parameters = Dune::PDELab::NavierStokesDefaultParameters<
      GV, RF, SourceFunction, BCTypeParam, BdrySolution, NeumannFlux, navier, tensor>;

  Parameters parameters(mu, rho, source_function, bdry_type, bdry_solution,
                           neumann_flux);


  driver<GV, BdrySolution, Parameters>(gv, parameters, bdry_solution, dt, T_time);

  return 0;
}
