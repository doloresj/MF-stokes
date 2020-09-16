/*
 * File:   driver.hh
 *
 */

#ifndef DRIVER_HH
#define	DRIVER_HH

#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>
//#include <dune/istl/superlu.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/seqistlsolverbackend.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/localoperator/taylorhoodnavierstokes.hh>
#include <dune/pdelab/newton/newton.hh>


template <int dim>
using Point = Dune::FieldVector<double, dim>;

template <int dim>
Point<dim> fun(Point<dim> const &x) {
  Point<dim> y(0.0);
  double a = 2.0;
  assert(dim <= 3);

  y[0] = 1.0;   //x[0] * x[0] * x[1] * x[1] * x[1] * std::cos(a * (x[0] + x[1]));
  if (dim > 1)
    y[1] = 1.0;     //(x[0] * x[0] + x[1] * x[1]) * std::exp(x[0] * x[1]);
  if (dim > 2)
    y[2] = 1.0;     //std::sin(a * x[2]) * std::exp(x[1]) * x[2] * x[2] * x[2] * x[2];
  return y;
}//ovo nam je samo testna fja na kojoj smo testirali radi li nam bdry_integration()... tu bi trebalo ici pod integral -p+grad{u}.. više manje jel i te vrijednost se spremaju i racunaju u driver fji varijabls je x_0


/* Integracija po rubu domene vektorskog polja f.n.
 * p = red točnosti integracijske formule
 */
template <typename Grid>
double bdry_integration(Grid &gridView, int p) {

  const int dim = Grid::dimension;
 // auto gridView = grid.leafGridView();
  double integral = 0.0;
  // Petlja po svim elementima
  for (auto const &element : elements(gridView)) {
    double elem_integral = 0.0;
    // petlja po svim stranicama elementa
    for (auto const &side : intersections(gridView, element)) {

      if (side.boundary()) // Jesmo li na granici domene?
      {
        const auto sidegeo = side.geometry();
        auto outerNormal = side.centerUnitOuterNormal();
        // Zatražimo kvadraturnu formulu na stranici (dimenzija dim -1)
        const auto &rule =
            Dune::QuadratureRules<double, dim - 1>::rule(sidegeo.type(), p);
        if (rule.order() < p)
             std::cerr << "Integracijska formula reda " << p << " nije dostupna.\n";
        double result = 0.0;
        // Petlja po svim integracijskim točkama
        for (auto const &qpoint : rule) {
          auto fval = fun(sidegeo.global(qpoint.position()));
          double weight = qpoint.weight();
          // | det (grad g) | daje Geometry objekt
          double detjac = sidegeo.integrationElement(qpoint.position());
          result += (fval * outerNormal) * weight * detjac;
        }

        elem_integral += result;
      }
    } // kraj petlje po svim stranicama
    integral += elem_integral;
  } // kraj petlje po svim elementima

  return integral;
}//doslovno samo copy paste od juraks

template<typename GV, typename IF, typename PARAMS>
void driver(const GV& gv,  PARAMS & parameters, IF & bdry_solution, double dt, double T_time)
{   
    using namespace Dune::PDELab;

    static const unsigned int dim = GV::dimension;


    typedef double RF;
    Dune::Timer timer;
    std::cout << "=== Initialize:" << timer.elapsed() << std::endl;
    timer.reset();
    double time = 0.0;

    // Konstrukcija prostora konačnih elemenata
    typedef typename GV::Grid::ctype DF;
    const int k = 2;
    const int q = 2 * k;  // preciznost integracijske formule lokalnog operatora

    // Taylor-Hoodovi elementi -- P2 za brzinu P1 za tlak
    typedef PkLocalFiniteElementMap<GV, double, double, k>      V_FEM;  // komponenta brzine
    typedef PkLocalFiniteElementMap<GV, double, double, k - 1 > P_FEM;  // tlak
    V_FEM vFem(gv);
    P_FEM pFem(gv);

    using CDC = ConformingDirichletConstraints;
    using VB = ISTL::VectorBackend<>;

    // Ova klasa direktno konstruira vektorske elemente u R^dim.
    // Prostor mrežnih funkcija za brzinu (vektorski):  V_h
    using V_GFS = VectorGridFunctionSpace<GV, V_FEM, dim, VB, VB, CDC>;
    V_GFS VGfs(gv, vFem);
    VGfs.name("velocity");
    // Prostor mrežnih funkcija za tlak (skalarni): W_h
    using P_GFS = GridFunctionSpace<GV, P_FEM, CDC, VB>;
    P_GFS pGfs(gv, pFem);
    pGfs.name("pressure");
    // Prostor V_h x W_h
    // LexicographicOrderingTag daje poredak varijabli: v1,v2,p
    using GFS = CompositeGridFunctionSpace<VB, LexicographicOrderingTag, V_GFS, P_GFS> ;
    GFS gfs(VGfs, pGfs);

    // Primjena Dirichletovih ograničenja
    using C = typename GFS::template ConstraintsContainer<double>::Type;
    C cg;
    cg.clear();

    // Određivanje Dirichletove granice. Ovdje se koriste pomoćne klase koje određuju
    // Dirichletovu granicu za svaku komponentu vektorske funkcije (v1,v2,p).
    using ScalarVelConstraints = StokesVelocityDirichletConstraints<PARAMS>;
    using VelocityConstraints = PowerConstraintsParameters<ScalarVelConstraints, dim>;
    using PressureConstraints = StokesPressureDirichletConstraints<PARAMS>;
    using Constraints = CompositeConstraintsParameters<VelocityConstraints, PressureConstraints>;

    ScalarVelConstraints scalarvelocity_constraints(parameters);
    VelocityConstraints  velocity_constraints(scalarvelocity_constraints);
    PressureConstraints  pressure_constraints(parameters);
    Constraints          bconst(velocity_constraints, pressure_constraints);

    // Odredi Dirichletova ograničenja
    constraints(bconst, gfs, cg);

    // Prostorni lokalni operator - definiran u Dune::PDELab-u.
    using LOP = TaylorHoodNavierStokes<PARAMS>;
    LOP lop(parameters, q);
    using TLOP = NavierStokesMass<PARAMS>;
    TLOP tlop(parameters);

    // Mrežni operator
    using MBE = ISTL::BCRSMatrixBackend<>;
    MBE mbe(5); // maksimalan broj ne-nul elemenat u retku (samo pretpostavka)
    using GOV = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, RF, RF, RF, C, C>;
    GOV gov(gfs, cg, gfs, cg, lop, mbe);
    using GOT = Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, RF, RF, RF, C, C>;
    GOT got(gfs, cg, gfs, cg, tlop, mbe);
    using GO = Dune::PDELab::OneStepGridOperator<GOV,GOT>;
    GO go (gov,got);

//samo dodano di poziva bdry_integration i ovaj dio bi trebao (barem tako mislimo) ici dolje u while petlju di se računaju nase vrijednosti

for (int l = 1; l < 5; ++l) {
    //gv.globalRefine(1);
    double bdry_integral = bdry_integration(gv, 5);
    std::cout << "elements=" << std::setw(8) << std::right << gv.size(0)
              << " bdry_integral = " << std::scientific << std::setprecision(12)
              << bdry_integral << "\n";
  }


    // Vektor koeficijenata i interpolacija rubnog uvjeta
    using U = typename GO::Traits::Domain;
    U x0(gfs);
    x0 = 0.0;
    bdry_solution.setTime(time);
    interpolate(bdry_solution, gfs, x0);
    //std::cout << "=== Finished interpolation:" << timer.elapsed() << std::endl;
    timer.reset();

    // Postavi sve stupnjeve slobode koji nisu Dirichletovi na nulu.
    set_shifted_dofs(cg, 0.0, x0);

    // Linear solver
    //using LS = Dune::PDELab::ISTLBackend_SEQ_BCGS_AMG_SSOR<GO>;
    using LS = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LS ls(false);

    typedef Dune::PDELab::Newton<GO,LS,U> PDESOLVER;
    PDESOLVER newton(go, x0, ls);
    newton.setReassembleThreshold(0.0);
    newton.setVerbosityLevel(2);
    newton.setMaxIterations(25);
    newton.setLineSearchMaxIterations(30);
    //  Izbor metode konačnih diferencija
     Dune::PDELab::OneStepThetaParameter<RF> method(1);  // implicit, theta = 1
     Dune::PDELab::OneStepMethod<RF,GO,PDESOLVER,U,U> osm(method,go,newton);
     osm.setVerbosityLevel(2);

    timer.reset();

    using VTKW =Dune::SubsamplingVTKWriter<GV>;
        VTKW vtkwriter(gv, 2);
        Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x0);  // new
        Dune::VTKSequenceWriter<GV> writer(std::make_shared<VTKW>(vtkwriter), "out");
        writer.write(time);
        U x1(gfs,0.0);         // sljedeći vremenski sloj -jednokoračna metoda x1=x^{n+1}, x0=x^{n}

        while (time < T_time-1e-8) {
            // postavi novo vrijeme u BC klasu
            bdry_solution.setTime(time+dt);
            cg.clear();
            Dune::PDELab::constraints(bconst,gfs,cg);

            osm.apply(time, dt, x0, bdry_solution, x1);                           // riješi sustav
            
           // std::cout << x0 << std::endl;            
            
      // e sada, tu bi po meni trebalo ici taj dio di integral računamo, sads valjda imamo tu neku novu varijablu tipa U x2 i nekako iz toga izvuci brzinu i tlak, dakle ako je dobro razumijemo U je "V_{h brzine}×V_{h tlaka} " taj dio me razumijemo pa smo stali
            U r(gfs);
            r = 0.;
            go.residual(x0, r);
            std::cout << "Final Residual: " << r.two_norm() << std::endl;
            // graphics

            x0 = x1;                                              // pripremi sljedeći vremenski korak
            time += dt;
            writer.write(time);
          }

}

#endif	/* DRIVER_HH */
