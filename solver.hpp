#ifndef SERIAL_HPP
#define SERIAL_HPP

void init(double *C0, double length, double width, int nx_, int ny_,
          double dt_, double D, int rank, int num_procs);
void step();
void transfer(double *result);
void free_memory();

#endif