# rk-sparse-ode-solver

This project provides C utilities and a Runge--Kutta based solver for sparse linear ODE systems derived from TASEP models. The tools generate matrices and vectors, convert between sparse formats and run explicit or embedded Runge--Kutta schemes. The solver is parallelized with OpenMP and offers several optimisation features.

## Tools

The following helper programs are included (see the header of each source file for details):

- **genmat** – generate a sparse matrix in COO text format. Example:
  ```bash
  ./genmat -n 4 -i zahlen.txt -o test_n_4.mat
  ```
- **gencrs** – convert a matrix from COO to CRS format:
  ```bash
  ./gencrs -i test_n_4.mat -o test_n_4.crs
  ```
- **gencrs_new** – convert a matrix from COO to an alternative CRS format:
  ```bash
  ./gencrs_new -i test_n_4.mat -o test_n_4_new.crs
  ```
- **gendia** – convert a matrix from COO to diagonal storage:
  ```bash
  ./gendia -i test_n_4.mat -o test_n_4_dia.mat
  ```
- **geninit** – create an initial vector file with a normal distribution:
  ```bash
  ./geninit -d 4 -n test_n_4.init
  ```
- **vis** – produce a text visualization of a matrix in COO format:
  ```bash
  ./vis -i test_n_4.mat -o test_n_4_vis_new.mat
  ```

## `rk_new` solver

`rk_new` solves sparse linear ODE systems using several Runge--Kutta variants. Supported features are summarised below (see `rk_new.c` for full usage information):

- Standard RK4 without an external tableau
- Explicit methods using a MAT or CRS Butcher tableau (e.g. RK4 or E5)
- Richardson extrapolation on explicit methods
- Embedded methods (RKF45, BS23) with user-specified tolerance
- Support for CRS, new CRS and diagonal matrix formats
- Multithreading with configurable CPU affinity
- Loop tiling optimisations (two versions with adjustable tile size)

Typical call for RK4 with a MAT tableau:
```bash
./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat
```

For more examples see the comments in `rk_new.c`.

### Building

Run `make` inside the `src` directory to build all utilities. The resulting binaries will be created in the same directory.

```bash
cd src
make
```

## License

This project is provided for research and educational purposes. See individual source files for copyright
information.
