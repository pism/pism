#include "SimpleConfig.hh"

namespace pism {

extern std::ostream &operator<<(std::ostream &out, SimpleConfig const &config)
{
  for (std::map<std::string, double>::const_iterator ii = config.doubles.begin(); ii != config.doubles.end(); ++ii) {
    out << ii->first << " = " << ii->second << std::endl;
  }

  for (std::map<std::string, std::string>::const_iterator ii = config.strings.begin(); ii != config.strings.end(); ++ii) {
    out << ii->first << " = " << ii->second << std::endl;
  }

  for (std::map<std::string, bool>::const_iterator ii = config.bools.begin(); ii != config.bools.end(); ++ii) {
    out << ii->first << " = " << ii->second << std::endl;
  }
  return out;
}

}
