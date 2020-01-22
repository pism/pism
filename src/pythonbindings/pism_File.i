%{
#include "util/io/File.hh"
#include "util/io/io_helpers.hh"
%}

%ignore pism::File::read_variable(const std::string &, const std::vector<unsigned int> &, const std::vector<unsigned int> &, double *) const;
%ignore pism::File::write_variable(const std::string &, const std::vector<unsigned int> &, const std::vector<unsigned int> &, const double *) const;

%include "util/io/IO_Flags.hh"
%include "util/io/File.hh"
%include "util/io/io_helpers.hh"

%extend pism::File
{
  void write_variable(const std::string &variable_name,
                      const std::vector<unsigned int> &start,
                      const std::vector<unsigned int> &count,
                      const std::vector<double> &data) const {
    $self->write_variable(variable_name, start, count, data.data());
  }

  std::vector<double> read_variable(const std::string &variable_name,
                                    const std::vector<unsigned int> &start,
                                    const std::vector<unsigned int> &count) {
    size_t size = 1;
    for (auto c : count) {
      size *= c;
    }

    std::vector<double> result(size);
    $self->read_variable(variable_name, start, count, result.data());

    return result;
  }
};
