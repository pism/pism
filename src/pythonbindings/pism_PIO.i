%{
#include "util/io/PIO.hh"
#include "util/io/io_helpers.hh"
%}

%include "util/io/IO_Flags.hh"
%include "util/io/PIO.hh"               // include before VariableMetadata
%include "util/io/io_helpers.hh"

%ignore pism::PIO::get_vara_double(const std::string &, const std::vector<unsigned int> &, const std::vector<unsigned int> &, double *) const;
%ignore pism::PIO::put_vara_double(const std::string &, const std::vector<unsigned int> &, const std::vector<unsigned int> &, const double *) const;

%extend pism::PIO
{
  void put_vara_double(const std::string &variable_name,
                       const std::vector<unsigned int> &start,
                       const std::vector<unsigned int> &count,
                       const std::vector<double> &data) const {
    $self->put_vara_double(variable_name, start, count, data.data());
  }

  std::vector<double> get_vara_double(const std::string &variable_name,
                                      const std::vector<unsigned int> &start,
                                      const std::vector<unsigned int> &count) {
    size_t size = 1;
    for (auto c : count) {
      size *= c;
    }

    std::vector<double> result(size);
    $self->get_vara_double(variable_name, start, count, result.data());

    return result;
  }
};
