#define CL_HPP_TARGET_OPENCL_VERSION 300

#ifdef __APPLE__
#include "cl.hpp" 
#else
#include <CL/opencl.hpp>
#endif

class CLManager {
  public:

    cl::CommandQueue queue;
    cl::Context context;

    cl::Buffer cpm;
    cl::Buffer numberofedges;

    int pde_AB;
    cl::Buffer pdeA;
    cl::Buffer pdeB;
    cl::Buffer diffco;

    cl::Program make_program(std::string filename, std::string head = "");
  
  private:
    cl::Device device;
    int make_context();
    bool context_prepared = false;
};

extern CLManager clm;
