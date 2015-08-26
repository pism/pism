/* Copyright (C) 2015 Florian Ziemen and the PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>

#include <netcdf>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace netCDF;

void write_xml_header();
void write_xml_reference(const std::string filename, double step = -9e9);
void write_xml_footer();
vector<double> interp_field(vector<double> & tx , vector<double> & ty , const string name , NcFile & nc);
vector<double> read_var(string variable, NcFile & file);
int parse_options(int ac, char * argv[], string & tracer_file, string & coordinates_file, vector <string> & input_files );

int main(int argc, char *argv[])
{
  string tracer_file_name, coordinates_file_name;
  vector <string> input_files;
  int retval = parse_options (argc, argv, tracer_file_name, coordinates_file_name, input_files);
  if (retval != 0 ) return retval;
  NcFile tracer_file (tracer_file_name, NcFile::read);
  NcFile coords_file (coordinates_file_name, NcFile::read);
  vector<string> variables;
  variables.push_back("tracer_x");
  variables.push_back("tracer_y");
  variables.push_back("tracer_z");
  variables.push_back("time");
  std::cerr<<"Reading tracers from " << tracer_file_name << std::endl;
  vector<vector<double> > creations;
  for (vector<string>::iterator v = variables.begin() ; v != variables.end() ; v++){

    creations.push_back(read_var(*v, tracer_file));}

  write_xml_header();

  for (vector<string>::iterator filename = input_files.begin() ; filename != input_files.end() ; filename++)
    {
      string outfile = *filename;
      outfile.replace(outfile.end()-3,outfile.end(), ".vtu");
			std::cerr<<"Reading from " << *filename << std::endl;
      NcFile extra_file (*filename, NcFile::read);

      vector<double>
				tx = read_var("tracer_x", extra_file),
				ty = read_var("tracer_y", extra_file),
				tz = read_var("tracer_z", extra_file),
				t_id = read_var("tracer_id", extra_file);

      vector<double>::iterator
	iid = t_id.begin(),
	ix  = tx.begin(),
	iy = ty.begin(),
	iz = tz.begin();

      for ( ; ix != tx.end() ; ){
	if ( *iz < 0 ) {
	  iid = t_id.erase(iid);
	  ix = tx.erase(ix);
	  iy = ty.erase(iy);
	  iz = tz.erase(iz);
	} else{
	  iid++;
	  ix++;
	  iy++;
	  iz++;
	}
      }

	std::cerr<<"Reading coordinates from " << "coordinates.nc" << std::endl;
	vector<double>
		lon = interp_field(tx, ty, "lon", coords_file),
		lat = interp_field(tx, ty, "lat", coords_file);

      vector<double> topgz (tx.size(), 0.);

      if(not extra_file.getVar("topg").isNull())
				topgz = interp_field(tx, ty, "topg", extra_file);
      else if(not extra_file.getVar("tracer_topg").isNull())
				topgz = read_var("tracer_topg", extra_file);



      vtkSmartPointer<vtkPoints> points =
				vtkSmartPointer<vtkPoints>::New();

      ix  = tx.begin();
      iy = ty.begin();
      iz = tz.begin();
      vector<double>::iterator
				itg = topgz.begin(),
				ilon = lon.begin(),
				ilat = lat.begin();

      for ( ;ilon != lon.end(); )	{
				const double r = 6.371e6 + (*iz + *itg)*100;
				*ilon *= M_PI/180.;
				*ilat *= M_PI/180.;
				points->InsertNextPoint(cos(*ilat)*cos(*ilon)*r, cos(*ilat)*sin(*ilon)*r, sin(*ilat)*r);
				ilon++; ilat++; iz++; itg++;
      }


      vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
				vtkSmartPointer<vtkUnstructuredGrid>::New();
      unstructuredGrid->SetPoints(points);

      size_t npoints = t_id.size();

      vtkSmartPointer<vtkFloatArray> floatArray =
	vtkSmartPointer<vtkFloatArray>::New();
      floatArray->SetNumberOfValues(npoints);
      floatArray->SetNumberOfComponents(1);
      floatArray->SetName("id");
      vector<double>::iterator i_id = t_id.begin();
      for(vtkIdType i = 0; i < npoints; i++)
	  floatArray->SetValue(i,* i_id++);
      unstructuredGrid->GetPointData()->AddArray(floatArray);

      vector<vector<double> >::iterator ic = creations.begin();
      for (vector<string>::iterator v = variables.begin() ; v != variables.end() ; v++)
	{
	  vtkSmartPointer<vtkFloatArray> fArray =
	    vtkSmartPointer<vtkFloatArray>::New();
	  fArray->SetNumberOfValues(npoints);
	  fArray->SetNumberOfComponents(1);
	  fArray->SetName((*v).c_str());
	  vector<double>::iterator i_id = t_id.begin();
	  for(vtkIdType i = 0; i < npoints; i++)
	    fArray->SetValue(i,(*ic)[*i_id++]);
	  unstructuredGrid->GetPointData()->AddArray(fArray);
	  ic++;
	}




  // Write file
			vector<double> time = read_var("time", extra_file);
			write_xml_reference(outfile, time.back());
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(outfile.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(unstructuredGrid);
#else
  writer->SetInputData(unstructuredGrid);
#endif
  writer->SetEncodeAppendedData(false);
  writer->Write();
    }
  write_xml_footer();
	std::cerr<< "\n\nPIPE STDOUT INTO A .pvd FILE TO FEED INTO PARAVIEW. ... > test.pvd\n";
  return EXIT_SUCCESS;
	    }

void write_xml_header(){
  std::cout<< "<?xml version=\"1.0\" ?>\n"
	   <<"    <VTKFile byte_order=\"LittleEndian\" type=\"Collection\" version=\"0.1\"> \n"
	   << " <Collection>\n" ;
}

void write_xml_reference(const std::string filename, double step){
  static double default_step = 0 ;
  if (step == -9e9) {
    step = default_step;
    default_step += 1 ;
  }
  std::cout<< "<DataSet file=\""<< filename << "\" group=\"\" part=\"0\" timestep=\"" << step << "\"/> \n";
}

void write_xml_footer(){
  std::cout<<"</Collection> \n</VTKFile>\n";
}

vector<double> read_var(string variable, NcFile & file){
    NcVar var = file.getVar(variable);

    int dc = var.getDimCount();
    size_t npoints =  var.getDim(dc-1).getSize();
    vector <double> data (npoints);
    var.getVar(&data[0] );
    return data;
}


  vector<double> interp_field(vector<double> & tx , vector<double> & ty , const string name , NcFile & nc){
  NcVar var = nc.getVar(name);
  vector<double> interpolated(tx.size());
  int dc = var.getDimCount();
  size_t nx = var.getDim(dc-2).getSize(),
         ny = var.getDim(dc-1).getSize();
  double data [nx][ny];
  var.getVar(&data[0][0]);

  vector<double> ix = read_var ("x", nc),
                 iy = read_var ("y", nc);
  const double x0 = ix[0],
               dx = ix[1]-ix[0],
               y0 = iy[0],
               dy = iy[1]-iy[0];

  vector<double>::iterator x = tx.begin(),
                           y = ty.begin();
  vector<double>::iterator d = interpolated.begin();
  for (;x != tx.end(); ){
    const double i = (*x - x0) / dx,
                 j = (*y - y0) / dy;
    const double ir = i - floor(i),
                 jr = j - floor(j);
    const int fi = floor(i)
      ,ci = ceil(i),
      fj = floor(j),
      cj = ceil(j);
      *d =
	(1-ir)*(1-jr) * data[fi][fj] +
	(  ir)*(1-jr) * data[ci][fj] +
	(  ir)*(  jr) * data[ci][cj] +
	(1-ir)*(  jr) * data[fi][cj];
    x++; y++; d++;
  }

  return interpolated;
}

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
	copy(v.begin(), v.end(), ostream_iterator<T>(os, "\n"));
	return os;
}

int  parse_options(int ac, char * av[], string & creation_file, string & coordinates_file, vector <string> & input_files)
	{

		try {
			po::options_description desc("Allowed options");
			desc.add_options()
				("help", "produce help message")
				("creation_file,c", po::value<string>(&creation_file),
				 "Tracer creation log")
				("coordinates_file,l", po::value<string>(&coordinates_file),
				 "file with coordinates")
				("input-file", po::value< vector<string> >(&input_files), "input file, can just come as normal argument without --input_file")
				;

			po::positional_options_description p;
			p.add("input-file", -1);

			po::variables_map vm;
			po::store(po::command_line_parser(ac, av).
								options(desc).positional(p).run(), vm);
			po::notify(vm);

			if (vm.count("help") || (coordinates_file == "") || (creation_file == "") || (!vm.count("input-file"))) {
				cerr << "Usage: options_description [options] input_files\n";
				cerr << desc;
				return -7;
			}

			cerr << "Input files are:\n"
					 << vm["input-file"].as< vector<string> >() << "\n";


			cerr << "Creation file is " << creation_file << "\n";

			cerr << "Coordinates file is " << coordinates_file << "\n";
		}
		catch(std::exception& e)
			{
				cout << e.what() << "\n";
				return 1;
			}
		return 0;
	}
