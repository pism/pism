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

using namespace std;
using namespace netCDF;

void write_xml_header();
void write_xml_reference(const std::string filename, double step = -9e9);
void write_xml_footer();
vector<double> interp_field(vector<double> & tx , vector<double> & ty , const string name , NcFile & nc);

vector<double> read_var(string variable, NcFile & file);


  
int main(int argc, char *argv[])
{
  NcFile tracer_file ("tracers_created.nc", NcFile::read);
  NcFile coords_file ("coordinates.nc", NcFile::read);
  vector<string> variables;
  variables.push_back("x");
  variables.push_back("y");
  variables.push_back("z");
  variables.push_back("time");
  vector<vector<double> > creations;
  for (vector<string>::iterator v = variables.begin() ; v != variables.end() ; v++){
    
    creations.push_back(read_var(*v, tracer_file));}
  
 
  // Parse command line arguments
  if(argc < 2)
    {
      std::cerr << "Required arguments: OutputFilename" << std::endl;
      return EXIT_FAILURE;
    }
  write_xml_header();
  
  
  for (int file_i = 1 ; file_i < argc ; file_i++)
    {
      string filename = argv[file_i];
      string outfile = filename;
      outfile.replace(outfile.end()-3,outfile.end(), ".vtu");
      NcFile extra_file (filename, NcFile::read);
      
      vector<double>
	tx = read_var("tracer_x", extra_file),
	ty = read_var("tracer_y", extra_file),
	tz = read_var("tracer_z", extra_file),
	t_id = read_var("tracer_id", extra_file),
	lon = interp_field(tx, ty, "lon", coords_file),
	lat = interp_field(tx, ty, "lat", coords_file);

      vector<double> topgz (tx.size(), 0.);

      if(not extra_file.getVar("topg").isNull())
	topgz = interp_field(tx, ty, "topg", extra_file);
	 


      vtkSmartPointer<vtkPoints> points =
	vtkSmartPointer<vtkPoints>::New();

      vector<double>::iterator ix  = tx.begin(),
	iy = ty.begin(),
	iz = tz.begin(),
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
  write_xml_reference(outfile, file_i);
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
    
    size_t npoints =  var.getDim(0).getSize();
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


  // #include <boost/program_options.hpp>
	    
 // namespace {
 //   const size_t ERROR_IN_COMMAND_LINE = 1;
 //     const size_t SUCCESS = 0;
 //     const size_t ERROR_UNHANDLED_EXCEPTION = 2;
     
 //     } // namespace
	      
	 
 // vector<string> parse_options(int argc, char* argv){
 //   vector<string> files;
 //     namespace po = boost::program_options;
 //       po::options_description desc("Options");
 // 	 desc.add_options()
 // 	   ("help,h", "Print help messages")
 // 	   ("tracer_file,t", po::value<string>(&tracer_file)->required(),  "Tracer file")
 // 	   ("coordinates,c", po::value<string>(&coord_file)->required(),
 // 	    "Coordinate file (if different from files parsed)");
 // 	   po::variables_map vm;
 // 	     po::store(po::parse_command_line(ac, av, desc), vm);
 // 	       po::notify(vm);
 // 		 if (vm.count("help")) {
 // 		   cout << desc << "\n";
 // 		   return files;
 // 		   }
		   
 // };
