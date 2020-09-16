#ifndef CG_STOKES_INITIAL_HH
#define CG_STOKES_INITIAL_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/localoperator/stokesparameter.hh>
#include <dune/pdelab/common/function.hh>


// Klasa koja određuje tip granice
class BCTypeParam
{
private:
  double time;
public:
    // Ova klasa daje indekse: DoNothing i VelocityDirichlet i StressNeumann
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
    using BoundaryCondition = Dune::PDELab::StokesBoundaryCondition;
  };


  template<typename I>
  inline void evaluate (const I & intersection,
                        const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord,
                        BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::coorddimension>
        xg = intersection.geometry().global( coord );
    if( xg[0] > 10-1e-6 )
      y = BC::DoNothing; //homogen Neumannov uvjet
    else
      y = BC::VelocityDirichlet;
  }
  template <typename T>
  void setTime(T t){
  time = t;
  }
};


// Inicijalni i Dirichletov rubni uvjet za brzinu
template<typename GV, typename RF, int dim>
class Velocity :
  public Dune::PDELab::AnalyticGridFunctionBase<
                             Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
                             Velocity<GV,RF,dim>
                                               >
{
private:
  RF time;
  //RF intensity; //stac

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Velocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  /*Velocity(const GV & gv, double intensity_ = 1.0) : BaseT(gv), intensity(intensity_) {
        time=0.0;
    }*/ //stac

  Velocity(const GV & gv) : BaseT(gv){}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  { 
    y[1]=0;

    if(x[0] < 1e-6){
    //  y[0] =-x[1]*x[1]*x[1]*(x[1]-4)/20;
     y[0]= (x[1])*(4-x[1])*(150/4);
     //  y[0] = (x[1] + 2)*(3 - x[1]) * (400/5);
      
    }
    else{
      y[0]=0;
     
    }

//krug
/*    if ((x[0]-5)*(x[0]-5)+(x[1]-2)*(x[1]-2)<1+1e-3){
      double a = x[0]-5;
      double b = x[1]-2;
      double norm_ab=std::sqrt(a*a+b*b);
      y[0]=b/norm_ab;
      y[1]=-a/norm_ab;
    }
*/
//kvadrat
/*
  if (x[0]<4+1e-6&&x[0]>4-1e-6&&x[1]<3+1e-6&&x[1]>1+1e-6){
    y[0]=3*(x[1]-1)*(x[1]-3);
    y[1]=0;
  }
*/
  //krugovi2
/*
    if ((x[0]-3)*(x[0]-3)+(x[1]-3)*(x[1]-3)<1+1e-3){
        double a = x[0]-3;
        double b = x[1]-3;
        double norm_ab=std::sqrt(a*a+b*b);
        y[0]=b/norm_ab;
        y[1]=-a/norm_ab;
      }
    if ((x[0]-3)*(x[0]-3)+(x[1]-1)*(x[1]-1)<1+1e-3){
        double a = x[0]-3;
        double b = x[1]-1;
        double norm_ab=std::sqrt(a*a+b*b);
        y[0]=b/norm_ab;
        y[1]=-a/norm_ab;
      }
    if ((x[0]-7)*(x[0]-7)+(x[1]-1)*(x[1]-1)<1+1e-3){
        double a = x[0]-7;
        double b = x[1]-1;
        double norm_ab=std::sqrt(a*a+b*b);
        y[0]=b/norm_ab;
        y[1]=-a/norm_ab;
        }
    if ((x[0]-5)*(x[0]-5)+(x[1]-2)*(x[1]-2)<1+1e-3){
        double a = x[0]-5;
        double b = x[1]-2;
        double norm_ab=std::sqrt(a*a+b*b);
        y[0]=b/norm_ab;
        y[1]=-a/norm_ab;
        }
    if ((x[0]-7)*(x[0]-7)+(x[1]-3)*(x[1]-3)<1+1e-3){
        double a = x[0]-7;
        double b = x[1]-3;
        double norm_ab=std::sqrt(a*a+b*b);
        y[0]=b/norm_ab;
        y[1]=-a/norm_ab;
        }
*/
//    if (time == 0) y[0] = x[1]*(4-x[1])/4; //inicijalni uvjet
  }


  template <typename T>
  void setTime(T t){
    time = t;
  }

};




// Vektorska funkcija jednaka nuli

template<typename GV, typename RF, std::size_t dim_range>
class ZeroFunction :
  public Dune::PDELab::AnalyticGridFunctionBase<
                            Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim_range>,
                            ZeroFunction<GV,RF,dim_range>
                                               >,
  public Dune::PDELab::InstationaryFunctionDefaults
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim_range> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ZeroFunction> BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  ZeroFunction(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    y=0;
  }
};





#endif
