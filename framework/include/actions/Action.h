/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef ACTION_H
#define ACTION_H

#include "InputParameters.h"
#include "ConsoleStreamInterface.h"

#include <string>
#include <ostream>

class Action;
class ActionWarehouse;
class ActionFactory;
class MooseMesh;
class FEProblem;
class Executioner;
class MooseApp;
class Factory;

template<>
InputParameters validParams<Action>();

/**
 * Base class for actions.
 */
class Action : public ConsoleStreamInterface
{
public:
  Action(const std::string & name, InputParameters params);
  virtual ~Action() {}                  // empty virtual destructor for proper memory release

  virtual void act() = 0;

  const std::string & name() const { return _name; }

  const std::string & type() const { return _action_type; }

  InputParameters & parameters() { return _pars; }

  const std::string & specificTaskName() const { return _specific_task_name; }

  const std::set<std::string> & getAllTasks() const { return _all_tasks; }

  ///@{
  /**
   * Retrieve a parameter for the object
   * @param name The name of the parameter
   * @return The value of the parameter
   */
  template <typename T>
  const T & getParam(const std::string & name);

  template <typename T>
  const T & getParam(const std::string & name) const;
  ///@}

  inline bool isParamValid(const std::string &name) const { return _pars.isParamValid(name); }

  inline InputParameters & getParams() { return _pars; }

  /**
   * Returns the short name which is the final string after the last delimiter for the
   * current ParserBlock
   */
  std::string getShortName() const;

  /**
   * Returns the base name which is the string before the last delimiter for the
   * current ParserBlock
   * Note:
   *  1) If there are multiple slashes, like ./foo/bar/baz, this will return
   *     ./foo/bar while getShortName() will return baz.
   *  2) If there are no slashes, this will return an empty string and
   *     getShortName() will return the action name.
   */
  std::string getBaseName() const;

  void appendTask(const std::string & task) { _all_tasks.insert(task); }

protected:
  /// The name of the action
  std::string _name;

  /// Input parameters for the action
  InputParameters _pars;

  // The registered syntax for this block if any
  std::string _registered_identifier;

  // The type name of this Action instance
  std::string _action_type;

  /// The MOOSE application this is associated with
  MooseApp & _app;

  /// The Factory associated with the MooseApp
  Factory & _factory;

  /// Builds Actions
  ActionFactory & _action_factory;

  /**
   * This member will only be populated if this Action instance is only designed to
   * handle one task.  This happens when an Action is registered with several pieces
   * of syntax in which case separate instances are built to handle the different
   * incoming parameter values.
   */
  std::string _specific_task_name;

  /**
   * A list of all the tasks that this Action will satisfy.
   * Note: That this is _not_ populated at construction time.  However, all tasks will be
   *       added prior to act().
   */
  std::set<std::string> _all_tasks;

  /// Reference to ActionWarehouse where we store object build by actions
  ActionWarehouse & _awh;

  /// The current action (even though we have seperate instances for each action)
  const std::string & _current_task;

  MooseSharedPointer<MooseMesh> & _mesh;
  MooseSharedPointer<MooseMesh> & _displaced_mesh;
  /// Convenience reference to a problem this action works on

public:
  MooseSharedPointer<FEProblem> & _problem;
  /// Convenience reference to an executioner
  MooseSharedPointer<Executioner> & _executioner;
};

template <typename T>
const T &
Action::getParam(const std::string & name)
{
  return InputParameters::getParamHelper(name, _pars, static_cast<T *>(0));
}

template <typename T>
const T &
Action::getParam(const std::string & name) const
{
  return InputParameters::getParamHelper(name, _pars, static_cast<T *>(0));
}

#endif // ACTION_H
