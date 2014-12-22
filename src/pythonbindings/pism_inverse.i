
/* Does not seem like this is used anywhere, but if I don't compile
   this, it will rot. */
%include "stressbalance/ssa/SNESProblem.hh"
%template(SNESScalarProblem) pism::SNESProblem<1,double>;
%template(SNESVectorProblem) pism::SNESProblem<2,pism::Vector2>;
%{
#include "SNESProblem.hh"
%}

/* Inverse model classes */
%{
#include "inverse/IP_SSATaucForwardProblem.hh"
#include "inverse/IP_SSAHardavForwardProblem.hh"
#include "inverse/IPDesignVariableParameterization.hh"
#include "inverse/functional/IPFunctional.hh"
#include "inverse/functional/IP_L2NormFunctional.hh"
#include "inverse/functional/IP_H1NormFunctional.hh"
#include "inverse/functional/IPGroundedIceH1NormFunctional.hh"
#include "inverse/functional/IPTotalVariationFunctional.hh"
#include "inverse/functional/IPMeanSquareFunctional.hh"
#include "inverse/functional/IPLogRelativeFunctional.hh"
#include "inverse/functional/IPLogRatioFunctional.hh"
#include "inverse/IP_SSATaucTikhonovGNSolver.hh"
#ifdef PISM_USE_TAO
#include "inverse/TaoUtil.hh"
#include "inverse/IP_SSATaucTaoTikhonovProblem.hh"
#include "inverse/IP_SSATaucTaoTikhonovProblemLCL.hh"
#include "inverse/IP_SSAHardavTaoTikhonovProblem.hh"
#endif
%}

%{
#include "TerminationReason.hh"
%}

%typemap(in, numinputs=0, noblock=1) pism::TerminationReason::Ptr & OUTPUT(pism::TerminationReason::Ptr temp) {
  $1 = &temp;
}

%typemap(argout,noblock=1) pism::TerminationReason::Ptr & OUTPUT
{
  {
    pism::TerminationReason::Ptr *smartresult = new pism::TerminationReason::Ptr(*$1);
    %append_output(SWIG_NewPointerObj(%as_voidptr(smartresult), $descriptor, SWIG_POINTER_OWN));
  }
};

%apply pism::TerminationReason::Ptr & OUTPUT { pism::TerminationReason::Ptr &reason };

%shared_ptr(pism::TerminationReason)
%shared_ptr(pism::KSPTerminationReason)
%shared_ptr(pism::SNESTerminationReason)
%shared_ptr(pism::GenericTerminationReason)
%include "TerminationReason.hh"


%include "inverse/functional/IPFunctional.hh"
%template(IPFunctional2S) pism::IPFunctional< pism::IceModelVec2S >;
%template(IPFunctional2V) pism::IPFunctional< pism::IceModelVec2V >;
%template(IPInnerProductFunctional2S) pism::IPInnerProductFunctional< pism::IceModelVec2S >;
%template(IPInnerProductFunctional2V) pism::IPInnerProductFunctional< pism::IceModelVec2V >;
%include "inverse/functional/IP_L2NormFunctional.hh"
%include "inverse/functional/IP_H1NormFunctional.hh"
%include "inverse/functional/IPGroundedIceH1NormFunctional.hh"
%include "inverse/functional/IPTotalVariationFunctional.hh"
%include "inverse/functional/IPMeanSquareFunctional.hh"
%include "inverse/functional/IPLogRatioFunctional.hh"
%include "inverse/functional/IPLogRelativeFunctional.hh"
%include "inverse/IPDesignVariableParameterization.hh"
%include "inverse/IP_SSATaucForwardProblem.hh"
%include "inverse/IP_SSAHardavForwardProblem.hh"
%include "inverse/IP_SSATaucTikhonovGNSolver.hh"

#ifdef PISM_USE_TAO
%ignore TaoConvergedReasons;
%shared_ptr(pism::TAOTerminationReason)
%include "inverse/TaoUtil.hh"

%include "inverse/IPTaoTikhonovProblem.hh"

//################### IP_SSATauc... #############################

// Instantiate the base class for IP_SSATaucTaoTikhonovProblem
// so that SWIG will implement the base class methods.
%template(IP_SSATaucTaoTikhonovProblemBaseClass) pism::IPTaoTikhonovProblem<pism::IP_SSATaucForwardProblem>;

%shared_ptr(pism::IPTaoTikhonovProblemListener<pism::IP_SSATaucForwardProblem>)

%feature("director") pism::IPTaoTikhonovProblemListener<pism::IP_SSATaucForwardProblem>;

%template(IP_SSATaucTaoTikhonovProblemListener)  pism::IPTaoTikhonovProblemListener<pism::IP_SSATaucForwardProblem>;

%include "inverse/IP_SSATaucTaoTikhonovProblem.hh"

%template(IP_SSATaucTaoTikhonovSolver) pism::TaoBasicSolver<pism::IP_SSATaucTaoTikhonovProblem>;


%shared_ptr(pism::IP_SSATaucTaoTikhonovProblemLCLListener)
%feature("director") pism::IP_SSATaucTaoTikhonovProblemLCLListener;
%include "inverse/IP_SSATaucTaoTikhonovProblemLCL.hh"
%template(IP_SSATaucTaoTikhonovProblemLCLSolver) pism::TaoBasicSolver< pism::IP_SSATaucTaoTikhonovProblemLCL >;


//################### IP_SSAHardav... #############################

%template(IP_SSAHardavTaoTikhonovProblemBaseClass) pism::IPTaoTikhonovProblem<pism::IP_SSAHardavForwardProblem>;

%shared_ptr(pism::IPTaoTikhonovProblemListener<pism::IP_SSAHardavForwardProblem>)

%feature("director") pism::IPTaoTikhonovProblemListener<pism::IP_SSAHardavForwardProblem>;

%template(IP_SSAHardavTaoTikhonovProblemListener) pism::IPTaoTikhonovProblemListener<pism::IP_SSAHardavForwardProblem>;

%include "inverse/IP_SSAHardavTaoTikhonovProblem.hh"

%template(IP_SSAHardavTaoTikhonovSolver) pism::TaoBasicSolver<pism::IP_SSAHardavTaoTikhonovProblem>;

#endif  /* end of ifdef PISM_USE_TAO */
