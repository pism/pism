#include <netcdf>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>



using namespace std;
using namespace netCDF;

const bool debug = false;

typedef struct {
  double id, x, y, z, topg;
}          Dataset;

void write_xml_header();
void write_xml_reference(const std::string filename, double step = -9e9);
void write_xml_footer();


vector<double> read_var(string variable, NcFile & file);
vector<double> interp_field(vector<double> & tx, vector<double> & ty, const string name, NcFile & nc);
std::vector<Dataset> find_tracers(std::vector<Dataset> tracers, std::vector<double> ids);

void print_tracer(const Dataset & p);
bool  lt (const Dataset a , const Dataset b);

std::vector<Dataset> read_file( string filename) ;
void write_data(const std::vector<Dataset> & data, double time, NcFile & file);

  
int main(int argc, char *argv[])
{
  std::vector<double> ids ;
  double id;
  while (cin >> id )
    ids.push_back(id);
  
  std::sort(ids.begin(), ids.end());


  NcFile tracer_file ("tracers_traced.nc", NcFile::replace);
  tracer_file.addDim("time");
  tracer_file.addDim("tracer_id", ids.size());
  NcVar idv = tracer_file.addVar("tracer_id", "double", "tracer_id");
  idv.putVar(&ids[0]);
  tracer_file.addVar("time", "double", "time");
  std::vector<std::string> dimensions;
  dimensions.push_back("time");
  dimensions.push_back("tracer_id");
  vector<string> variables;
  variables.push_back("tracer_x");
  variables.push_back("tracer_y");
  variables.push_back("tracer_z");
  variables.push_back("tracer_topg");
  vector<vector<double> > creations;
  for (vector<string>::iterator v = variables.begin() ; v != variables.end() ; v++){
    NcVar cv = tracer_file.addVar(*v, "float", dimensions);
    cv.putAtt("_FillValue", ncFloat, -9e9);
  }
  
  // Parse command line arguments
  if(argc < 2)
    {
      std::cerr << "Required arguments: OutputFilename" << std::endl;
      return EXIT_FAILURE;
    }
  
  for (int file_i = 1 ; file_i < argc ; file_i++)
    {
      string filename = argv[file_i];
      vector<double> time ; 
      {
	NcFile extra_file (filename, NcFile::read);
	time = read_var("time", extra_file );
	if (file_i == 1){
	  std::string timeunit, calendar;
	  extra_file.getVar("time").getAtt("units").getValues(timeunit);
	  extra_file.getVar("time").getAtt("calendar").getValues(calendar);
	  tracer_file.getVar("time").putAtt("units", timeunit);
	  tracer_file.getVar("time").putAtt("calendar", calendar);
	}
      }
      
      vector<Dataset> particles = read_file(filename);
      std::sort(particles.begin(), particles.end(), lt);
      vector<Dataset> particles_found = find_tracers(particles, ids);
      
      write_data(particles_found, *time.begin(), tracer_file);
    }
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
  if (debug)    std::cerr<< "reading " << variable << " ... ";
    NcVar var = file.getVar(variable);
    if (var.isNull())
      {
	vector<double> data ;
	return data;
      }
    size_t npoints =  var.getDim(0).getSize();
    vector <double> data (npoints);
    var.getVar(&data[0] );
    if (debug) std::cerr<< "done. \n" ;
    return data;
}

  
  vector<double> interp_field(vector<double> & tx , vector<double> & ty , const string name , NcFile & nc){
    if (debug)    std::cerr<< "interpolating " << name << "\n";
    if (tx.size() == 0 )
    {
      vector<double> data ;
      return data;
    }
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


bool  lt (const Dataset a , const Dataset b){
  return a.id < b.id;
}

void print_tracer(const Dataset & p){
  std::cout <<" Dataset " << p.id << "\t x="<< p.x << "\t y="<< p.y << "\t z="<< p.z<<"\t topg" << p.topg <<  "\n";
}
 
 std::vector<Dataset> read_file( string filename) {
   NcFile extra_file (filename, NcFile::read);
   if (debug) std::cerr<< "reading tracer information from " << filename << "\n";
   vector<double>
     tx = read_var("tracer_x", extra_file),
     ty = read_var("tracer_y", extra_file),
     tz = read_var("tracer_z", extra_file),
     t_id = read_var("tracer_id", extra_file);
   std::vector<double>     topgz;
   if (not extra_file.getVar("topg").isNull())
     topgz = interp_field(tx, ty, "topg", extra_file);
   else
     topgz.resize(tx.size(), -9e9);
   
   vector<Dataset> result(tx.size());
   
   vector<double>::iterator
     ix  = tx.begin(),
     iy = ty.begin(),
     iz = tz.begin(),
     itg = topgz.begin(),
     iid = t_id.begin();
   for (   std::vector<Dataset>::iterator ip = result.begin(); ip != result.end(); ip++)
     {
       ip->id = *iid++;
       ip->x = *ix++;
       ip->y = *iy++;
       ip->z = *iz++;
       ip->topg = *itg++;
     }
   return result;
 }

std::vector<Dataset> find_tracers(std::vector<Dataset> tracers, std::vector<double> ids){
  Dataset nothing={-9e9, -9e9, -9.e9, -9e9};
  std::vector<Dataset> c(ids.size(), nothing);
  std::vector<Dataset>::iterator it = tracers.begin(),
    ic = c.begin();
  std::vector<double>::iterator ii = ids.begin();
  for ( ; it != tracers.end() && ii != ids.end() &&  ic != c.end(); ){
    if (it->id < *ii){
      //      std::cout<< it-> id << "\t< " << *ii << std::endl;
      it++;
      continue;
    }
    else if (it->id > *ii){
      ic->id = *ii;
      //      std::cout<< it-> id << "\t> " << *ii << std::endl;
      ic++;
      ii++;
      continue;
    }
    *ic = *it;
    //    std::cout<< it-> id << "\t= " << *ii << std::endl;
    it++; ii++; ic++;
  }
  for ( ; ic != c.end() && ii != ids.end() ; ){
    ic->id = *ii ;
    ic++; ii++;
  }
  return c;
}
 
void write_data(const std::vector<Dataset> & data, double time, NcFile & file) {
     if (debug) std::cerr<< "Writing tracer information\n";

   vector<float>
     tx    (data.size()),
     ty    (data.size()),
     tz    (data.size()),
     t_id  (data.size()),
     topgz (data.size());
   
   vector<float>::iterator
     ix  = tx.begin(),
     iy  = ty.begin(),
     iz  = tz.begin(),
     itg = topgz.begin(),
     iid = t_id.begin();
   for (   std::vector<Dataset>::const_iterator ip = data.begin(); ip != data.end(); ip++){
     *iid++ = ip->id;
     *ix++ = ip->x ;
     *iy++ = ip->y; 
     *iz++ = ip->z;
     *itg++ = ip->topg;
   }
   
   std::vector<double*> pointers ;
   std::vector<std::string> names ;
   vector<size_t> offsets, counts;
   if (debug) std::cerr<< "Writing time\n";
   offsets.push_back(file.getDim("time").getSize());
   file.getVar("time").putVar(offsets, time);
   offsets.push_back(0);
   counts.push_back(1);
   counts.push_back(tx.size());
   if (debug) std::cerr<< "tracer data\n";
   file.getVar("tracer_x").putVar(offsets, counts,  &tx[0]);
   file.getVar("tracer_y").putVar(offsets, counts, &ty[0]);
   file.getVar("tracer_z").putVar(offsets, counts, &tz[0]);
   file.getVar("tracer_topg").putVar(offsets, counts, &topgz[0]);
   
}
