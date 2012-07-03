// Copyright (C) 2012  David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef TERMINATIONREASON_HH_JW17MC8V
#define TERMINATIONREASON_HH_JW17MC8V

#include <petscsnes.h>
#include <petscksp.h>
#include <tr1/memory>
#include <string>
#include <sstream>

class TerminationReason {
public:
  TerminationReason() :m_reason(0) {};
  
  typedef std::tr1::shared_ptr<TerminationReason> Ptr;
  
  virtual int reason() {
    return m_reason;
  };

  virtual std::string description() {
    std::stringstream sdesc;
    this->get_description(sdesc);
    return sdesc.str();
  }
  virtual void get_description( std::ostream &desc,int indent_level=0) = 0;

  virtual std::string nested_description() {
    std::stringstream sdesc;
    this->get_nested_description(sdesc);
    return sdesc.str();
  }
  virtual void get_nested_description( std::ostream &desc,int indent_level=0) {
    this->get_description(desc,indent_level);
    if(this->has_root_cause()) {
      indent_level++;
      desc << std::endl;
      this->root_cause()->get_nested_description(desc,indent_level);
    }
  }

  virtual bool has_root_cause() {
    return m_root_cause;
  };

  TerminationReason::Ptr root_cause() {
    return m_root_cause;
  };

  void set_root_cause(TerminationReason::Ptr cause) {
    m_root_cause = cause;
  }
  
  bool succeeded() {
    return (this->reason())>0;
  }
  
  bool failed() {
    return (this->reason())<0;
  }
  
protected:
  int m_reason;
  TerminationReason::Ptr m_root_cause;
  static const char *sm_indent;

private:
  TerminationReason(TerminationReason const &reason);
  TerminationReason &operator =(TerminationReason const &reason);
};

class KSPTerminationReason: public TerminationReason {
public:
  KSPTerminationReason( KSPConvergedReason r);
  virtual void get_description( std::ostream &desc,int indent_level=0);
};

class SNESTerminationReason: public TerminationReason {
public:
  SNESTerminationReason( SNESConvergedReason r);
  virtual void get_description( std::ostream &desc,int indent_level=0);
};

class GenericTerminationReason: public TerminationReason {
public:
  GenericTerminationReason(int code, std::string &desc) :
  m_description(desc) {
    m_reason = code;
  };

  GenericTerminationReason(int code, const char*desc) :
  m_description(desc) {
    m_reason = code;
  };
  
  virtual void get_description( std::ostream &desc, int indent_level=0); 
protected:
  std::string m_description;
};

#endif /* end of include guard: TERMINATIONREASON_HH_JW17MC8V */
