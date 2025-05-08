#include <cstdio>
#include <main.h>
#include <gperftools/profiler.h>
int main(const int argc, const char **argv) {

  if(!ProfilerStart("/tmp/nbody.prof")){
    fprintf(stderr, "could not start profiler");
    return 1;
  }
  int er =  libMain(argc, argv);
  ProfilerStop();
  return er;
}
