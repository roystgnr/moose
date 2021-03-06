#include "NodalAreaVarAction.h"

#include "Factory.h"
#include "FEProblem.h"
#include "Parser.h"
#include "MooseApp.h"
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<NodalAreaVarAction>()
{
  MooseEnum orders("CONSTANT FIRST SECOND THIRD FOURTH", "FIRST");

  InputParameters params = validParams<Action>();
  params.addParam<MooseEnum>("order", orders, "The finite element order: " + orders.getRawNames());
  return params;
}

NodalAreaVarAction::NodalAreaVarAction(const std::string & name, InputParameters params) :
  Action(name, params)
{
}

void
NodalAreaVarAction::act()
{
  std::string short_name(_name);
  // Chop off "Contact/"
  short_name.erase(0, 8);

  _problem->addAuxVariable("nodal_area_"+ short_name,
                           FEType(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
                                  Utility::string_to_enum<FEFamily>("LAGRANGE")));

}
