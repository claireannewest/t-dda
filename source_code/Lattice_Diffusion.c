#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h> 
#include "lattice_prototypes.h"

//	This code implements the T-DDA. Author: Chris Baldwin, enhanced by Claire A. West.
//	This program computes steady-state temperature distributions for an externally-illuminated target
//  (represented by a finite list of points on a cubic lattice) immersed in a homogeneous background.
//	The syntax for running this program is as follows: 
//  'Lattice_Diffusion <Green's function value list> <parameter list> <input file name> <output file name>'
//	What should be in each input is as follows:
//		-<Green's funion value list>: A list of the lattice Green's function values for a cubic lattice, 
//									  as described in the paper "Application of the Lattice Green's Function 
//									  for Calculating the Resistance of an Infinite Network of Resistors" by 
// 									  Jozsef Cserti, published in AJP, volume 68, pg. 896, in 2000.
//									  The list should be formatted as follows: "%d %d %d %lf\n", where the 
// 									  first three integers (which should be non-negative) give the position of 
// 									  each point relative to a source at the origin, and the fourth number gives 
//									  the value of the lattice Green's function at that point.
//		-<parameter list>: 	A list of the parameter values to be used in the calculation. A description of the 
//							parameters is given below.
//		-<input file name>: The location of the file containing the target geometry and heat source.
//							This file should be formatted in one of two ways, depending on the value of 'input_mode' 
//							in the parameter list.
//							If input_mode == 1, then each line of this file should be of the form: 
//							"%lf %lf %lf %d %lf %lf %lf %lf %lf %lf\n". The first three floats are the positions of each point in the target (in units of the lattice spacing), the next int is an index for the composition of that point, and the next six floats give the real and imaginary parts of each component of the electric field at that point (which should include the scattered electric field).
//							If input_mode != 1, then each line of this file should be of the form: "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n". The first three floats are the positions of each point in the target, and the next six floats are the real and imaginary parts of each component of the electric field at that point (which should include the scattered electric field).
//		-<output file name>: The location of the file to write the output into (existing files will be overwritten).
//							 Each line of the output file is of the form: "%d %d %d %lf %lf %lf\n", where the first three integers are the positions of each point in the target (in units of the lattice spacing), the next float is the temperature at that point, then the original heat source at that point, and finally the "effective heat source" at that point, i.e., the heat source including the thermal equivalent of "bound charges".

int main(int argc, char *argv[]) {
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	printf("Start time:\n%02d:%02d:%02d\n\n", tm.tm_hour, tm.tm_min, tm.tm_sec);


    if (argc != 5) {
        fprintf(stderr,
            "Invalid syntax.\nUsage:\n  %s <greens.dat> <params.txt> <input.txt> <output.txt>\n",
            argv[0]);
        return 1;
    }


	const char *dat_name = argv[1];
	const char *par_name = argv[2];
	const char *inp_name = argv[3];
	const char *out_name = argv[4];

    FILE *dat_id = fopen(dat_name, "r"); 
    if (!dat_id) { perror("Open Green's file"); return 1; }

    FILE *par_id = fopen(par_name, "r");
    if (!par_id) { perror("Open parameter file"); fclose(dat_id); return 1; }

    FILE *inp_id = fopen(inp_name, "r");  
    if (!inp_id) { perror("Open input file"); fclose(dat_id); fclose(par_id); return 1; }

    FILE *out_id = fopen(out_name, "w"); 
    if (!out_id) { perror("Open output file"); fclose(dat_id); fclose(par_id); fclose(inp_id); return 1; }



    int i = 0, j = 0, k = 0; // Related to G


	//	Reading in parameters.
	int num_k = 0; // Number of materials that make up the target, e.g., for a target composed of a homogeneous material set num_k = 1 in the "parameter list" file
	double *kappa = NULL; // kappa is the list of conductivities (the first element is the background conductivity)
	double n_m = 0.0; // background medium's index of refraction
	double lambda = 0.0; // the illuminating light's wavelength
	double I_0 = 0.0; // the intensity of the illuminating beam
	double unit = 0.0; // dipole spacing

	double x_plane = 0.0; 
	int d = 0; // the lattice spacing in units of "unit"
	int x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0; //the size of the grid on which calculations are done (the entire target should be included in the region defined by these values). 
	int input_mode = 0; // = 1 read in E's & P's from DDA; = 2 read in heat_pow per particle 
	int total_points = 0; // total number of points where temperature will be calculated


	//	Read in parameter file. 
	fscanf(par_id, "num_k: %d\n", &num_k);
	kappa = (double *)malloc((num_k + 2)*sizeof(double));
	fscanf(par_id, "k_out: %lf\n", kappa);
	for (i = 0; i < num_k; i++) {
	  fscanf(par_id, "k_in: %lf\n", kappa + i + 1);}
	fscanf(par_id, "k_sub: %lf\n", kappa + num_k + 1);
	fscanf(par_id, "lambda: %lf\nn_m: %lf\n", &lambda, &n_m);
	fscanf(par_id, "I_0: %lf\nunit: %lf\n\n", &I_0, &unit);
	fscanf(par_id, "d: %d\nx_min: %d\nx_max: %d\ny_min: %d\ny_max: %d\nz_min: %d\nz_max: %d\n\n", &d, &x_min, &x_max, &y_min, &y_max, &z_min, &z_max);
	fscanf(par_id, "x_plane: %lf\n", &x_plane);
	fscanf(par_id, "input_mode: %d\n", &input_mode);
	fscanf(par_id, "total_points: %d\n", &total_points);

    // Sanity checks that prevent later divide-by-zero or sign bugs 
    if (d <= 0)              { fprintf(stderr, "d must be > 0\n"); return 1; }
    if (unit <= 0.0)         { fprintf(stderr, "unit must be > 0\n"); return 1; }
    if (*kappa <= 0.0)       { fprintf(stderr, "k_out must be > 0\n"); return 1; }
    if (*(kappa + num_k + 1) <= 0.0) { fprintf(stderr, "k_sub must be > 0\n"); return 1; }
    if (input_mode != 1 && input_mode != 2) { fprintf(stderr, "input_mode must be 1 or 2\n"); return 1; }

	//  Calculating Substrate Ratios
	double upper = (kappa[num_k + 1] - kappa[0]) / (kappa[num_k + 1] + kappa[0]);
	double lower = (2.0*kappa[0])/(kappa[num_k + 1] + kappa[0]);
    printf("Upper = %.6f\n", upper);
    printf("Lower = %.6f\n", lower);

	//	Approximate numbers of coordinates to run calculation on
	int num_x = (x_max - x_min)/d + 1;
	int num_y = (y_max - y_min)/d + 1;
	int num_z = (z_max - z_min)/d + 1;
		
	//	Creating arrays of all possible x coordinates, y coordinates, z coordinates
	int *x_coords = (int *)malloc(num_x*sizeof(int));
	int *y_coords = (int *)malloc(num_y*sizeof(int));
	int *z_coords = (int *)malloc(num_z*sizeof(int));
	for (i = 0; i < num_x; i++) {
		x_coords[i] = x_min + i*d;
	}
	for (j = 0; j < num_y; j++) {
		y_coords[j] = y_min + j*d;
	}
	for (k = 0; k < num_z; k++) {
		z_coords[k] = z_min + k*d;
	}
	
	//	Reading in lattice Green's function values
	int x_lim, y_lim, z_lim;
	fscanf(dat_id, "%d %d %d", &x_lim, &y_lim, &z_lim);
	printf("Reading Green's grid of size %d x %d x %d\n", x_lim, y_lim, z_lim);

	// Allocate memory after reading header
	double *G_r = calloc((size_t)x_lim * y_lim * z_lim, sizeof(double));

	int ii, jj, kk;
	double dat_val;

	// Now read the rest
	while (fscanf(dat_id, "%d %d %d %lf", &ii, &jj, &kk, &dat_val) == 4) {
	    if (ii >= x_lim || jj >= y_lim || kk >= z_lim) {
	        printf("Index out of bounds in Green's function: i=%d j=%d k=%d (limits %d %d %d)\n",
	               ii, jj, kk, x_lim, y_lim, z_lim);
	        exit(1);
	    }
	    *(G_r + ii + x_lim*jj + x_lim*y_lim*kk) = dat_val;
	}

	time_t tG = time(NULL); struct tm tmG = *localtime(&tG);
	printf("Completed reading lattice Green's: %02d:%02d:%02d\n\n", tmG.tm_hour, tmG.tm_min, tmG.tm_sec);


	// Is the Green grid large enough?
	if (num_x > x_lim || num_y > y_lim || num_z > z_lim) {
	    fprintf(stderr, "ERROR: shape exceeds Green's grid "
	                    "(need %d×%d×%d, have %d×%d×%d)\n",
	            num_x, num_y, num_z, x_lim, y_lim, z_lim);
	    free(G_r); return 1;
	}


	//	Declaring useful variables
	int    N_init = 0; // Counts the number of target points in the input file
	int    N_init_guess = total_points; // Initial guess as to how many target points are in the input file
	int   *initial_reading = malloc((size_t)N_init_guess * sizeof *initial_reading);
	int   *material_ir     = malloc((size_t)N_init_guess * sizeof *material_ir);
	double*Q_ir            = malloc((size_t)N_init_guess * sizeof *Q_ir);
	

	// Local vars for parsing 
	int x, y, z, m_val = 0;
	double x_f, y_f, z_f;
	double E_xr, E_xi, E_yr, E_yi, E_zr, E_zi;
	double P_xr, P_xi, P_yr, P_yi, P_zr, P_zi;
	double Q_byhand;

	// Tight bounds on target region
	int x_tar_min = INT_MAX, x_tar_max = INT_MIN; 
	int y_tar_min = INT_MAX, y_tar_max = INT_MIN; 
	int z_tar_min = INT_MAX, z_tar_max = INT_MIN; 
	

	// Precompute constants used for Q 
	const double k0    = 2.0 * M_PI * n_m / lambda;      // 1/length
	const double voxel = unit * d;                       // length
	const double voxel3 = voxel * voxel * voxel;
	const double voxel6 = voxel3 * voxel3;


	// Read input lines
	for (;;) {  
	    if (input_mode == 1) {
	        int got = fscanf(inp_id,
	            "%lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	            &x_f, &y_f, &z_f, &m_val,
	            &E_xr, &E_xi, &E_yr, &E_yi, &E_zr, &E_zi,
	            &P_xr, &P_xi, &P_yr, &P_yi, &P_zr, &P_zi);
	        if (got != 16) break;  
	    } else if (input_mode == 2) {
	        int got = fscanf(inp_id, "%lf %lf %lf %lf", &x_f, &y_f, &z_f, &Q_byhand);
	        if (got != 4) break;   
	        m_val = 1;
	    } else {
	        fprintf(stderr, "ERROR: unsupported input_mode = %d\n", input_mode);
        	return 1;
	    }

		x = (x_f >= 0 ? (int)(x_f + 0.5) : (int)(x_f - 0.5));
		y = (y_f >= 0 ? (int)(y_f + 0.5) : (int)(y_f - 0.5));
		z = (z_f >= 0 ? (int)(z_f + 0.5) : (int)(z_f - 0.5));

		if (x < x_tar_min) x_tar_min = x; //We've found a point outside the existing bound on the region the target is in, so extend the region to contain this new point.
		if (x > x_tar_max) x_tar_max = x; //'' '' '' '' ''...
		if (y < y_tar_min) y_tar_min = y;
		if (y > y_tar_max) y_tar_max = y;
		if (z < z_tar_min) z_tar_min = z;
		if (z > z_tar_max) z_tar_max = z;
		
		int i = (x - x_min)/d;
		int j = (y - y_min)/d;
		int k = (z - z_min)/d;

		// Ensure big enough before writing
	    if (N_init >= N_init_guess) {
	        int new_guess = N_init_guess*2;
	        int   *ir  = realloc(initial_reading, (size_t)new_guess * sizeof *ir);
	        int   *mi  = realloc(material_ir,     (size_t)new_guess * sizeof *mi);
	        double*qi  = realloc(Q_ir,            (size_t)new_guess * sizeof *qi);
	        if (!ir || !mi || !qi) {
	            fprintf(stderr, "ERROR: realloc failed while reading input\n");
	            exit(1);
	        }
	        initial_reading = ir;
	        material_ir     = mi;
	        Q_ir            = qi;
	        N_init_guess    = new_guess;
	    }

		initial_reading[N_init] = k + (num_z + 1)*(j + (num_y + 1)*i); //Each triplet giving a target point (i.e., (x, y, z)) is compressed into a single integer index, then stored
		material_ir[N_init] = m_val;

		if (input_mode == 1) {
        	const double c_ext_prefactor = 4.0 * M_PI * k0 * I_0 * voxel3 * 1e-18;
        	const double k0_4 = k0 * k0 * k0 * k0;
        	const double c_sca_prefactor = (8.0 * M_PI / 3.0) * k0_4 * I_0 * voxel6 * 1e-45;
	        double ImE_dot_P = (E_xr*P_xi - E_xi*P_xr)
	                         + (E_yr*P_yi - E_yi*P_yr)
	                         + (E_zr*P_zi - E_zi*P_zr);
        	double PdotP_im  = (P_xr*P_xi) + (P_yr*P_yi) + (P_zr*P_zi);
	        Q_ir[N_init]  = c_ext_prefactor * ImE_dot_P;
	        Q_ir[N_init] -= c_sca_prefactor * PdotP_im;
        }

		if (input_mode == 2) {
        	Q_ir[N_init] = Q_byhand;
		}

		N_init++;

	}
	
	time_t t1 = time(NULL);
	struct tm tm1 = *localtime(&t1);
	printf("Completed reading input values:\n%02d:%02d:%02d\n\n", tm1.tm_hour, tm1.tm_min, tm1.tm_sec);

	int l, index;
	size_t n, p;

	// 	Sorting initial_reading (and material_ir and Q_ir) in order of increasing target point index, 
	//	so that initial_reading contains entries in increasing order.
	// 	This makes it easier to search for specific values in the array later.
	double q_val;
	int node_num; // Node #, i.e., location in array, of minimum index

	for (n = 0; n < N_init; n++) {
		node_num = n;
		index = *(initial_reading + node_num);

		for (l = n + 1; l < N_init; l++) {
			if (*(initial_reading + l) < index) { //Found new lowest index, this becomes the new minimum
				node_num = l;
				index = *(initial_reading + node_num);
			}
		}

		m_val = *(material_ir + node_num);
		q_val = *(Q_ir + node_num);
		*(initial_reading + node_num) = *(initial_reading + n);
		*(material_ir + node_num) = *(material_ir + n);
		*(Q_ir + node_num) = *(Q_ir + n);
		*(initial_reading + n) = index; //Put lowest target point index in the next spot in array
		*(material_ir + n) = m_val;
		*(Q_ir + n) = q_val;
	}
	



	//	Constructing a list of all points that need to have a modified thermal conductivity with at least one of their neighbors, 
	//	i.e., points that are either inside the target or directly adjacent to points that are inside the target. 
	//	This constitutes the "extended target": the original target plus a surrounding layer of background points.

	//	Also constructing a list of points in the background medium at which we'd like to know the steady-state temperature.
	
	int N = N_init, N_out = 0;
	int N_guess = N_init_guess, N_out_guess = 5000;
	int *target_points = (int *)malloc(N_guess*sizeof(int)); // List of indices for points that are in the extended target, i.e., either part of the target or have at least one neighbor that's part of the target
	int *outside_points = (int *)malloc(N_out_guess*sizeof(int));
	



	for (n = 0; n < N_init; n++) {
		*(target_points + n) = *(initial_reading + n); // All points from the initial reading of the input file are target points, so they should be in the list of extended target points
	}





	for (i = 0; i < num_x; i++) {
		x = x_coords[i];

		for (j = 0; j < num_y; j++) {
			y = y_coords[j];

			for (k = 0; k < num_z; k++) {
				z = z_coords[k];
				
				// We should check to see if the current point (given by x, y, z) is even close to the target. If not, there's no chance that this point is part of the extended target, so we shouldn't waste time checking for it.
				if (x_tar_min - 1 <= x && x <= x_tar_max + 1 && y_tar_min - 1 <= y && y <= y_tar_max + 1 && z_tar_min - 1 <= z && z <= z_tar_max + 1) { 
					// If the current point (given by x, y, z) isn't in the list of target points from the input file...
					if (find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*i)) < 0) { 
						if (find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*(i + 1))) >= 0 || // Adjacent point at (x + d, y, z)
							find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*(i - 1))) >= 0 || // Adjacent point at (x - d, y, z)
							find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + 1 + (num_y + 1)*i)) >= 0 ||   // Adjacent point at (x, y + d, z)
							find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j - 1 + (num_y + 1)*i)) >= 0 ||   // Adjacent point at (x, y - d, z)
							find_neighbor(initial_reading, N_init, k + 1 + (num_z + 1)*(j + (num_y + 1)*i)) >= 0 ||
							find_neighbor(initial_reading, N_init, k - 1 + (num_z + 1)*(j + (num_y + 1)*i)) >= 0) {
								// And if at least one of the adjacent points is in the list of target points from the input file... 
								*(target_points + N) = k + (num_z + 1)*(j + (num_y + 1)*i); // The current point is also part of the extended target
								*(material_ir + N) = 0; // The current point is technically still part of the background medium, so it has material index 0
								*(Q_ir + N) = 0.0;
								N++;

								if (N >= N_guess) { // Arrays have run out of memory, reallocate more for them
									target_points = (int *)realloc(target_points, 2*N_guess*sizeof(int));
									material_ir = (int *)realloc(material_ir, 2*N_guess*sizeof(int));
									Q_ir = (double *)realloc(Q_ir, 2*N_guess*sizeof(double));
									N_guess *= 2;
								}
						// See if the current point, although far away from the target, is still one that we'd like to know the temperature at. 
						// This should be manually changed to get the temperature at different outside points.		
						} else if (x >= x_min && x <= x_max) { 
							*(outside_points + N_out) = k + (num_z + 1)*(j + (num_y + 1)*i);
							N_out++;
							if (N_out >= N_out_guess) { // Array has run out of memory, reallocate more
								outside_points = (int *)realloc(outside_points, 2*N_out_guess*sizeof(int));
								N_out_guess *= 2;
							}
						}

					}
				} else if (x >= x_min && x <= x_max) { //See if the current point, although far away from the target, is still one that we'd like to know the temperature at. This should be manually changed to get the temperature at different outside points.
					*(outside_points + N_out) = k + (num_z + 1)*(j + (num_y + 1)*i);
					N_out++;
					if (N_out >= N_out_guess) { //Array has run out of memory, reallocate more
						outside_points = (int *)realloc(outside_points, 2*N_out_guess*sizeof(int));
						N_out_guess *= 2;
					}
				}
			}
		}
	}
	
//	Determining the composition of each neighbor to each node 
	int *neighbor_material_ir = (int *)calloc(6*N, sizeof(int)); //6xN array listing the material composition of each neighbor to each node, i.e., *(neighbor_material_ir + n + i*N) gives the material index of the (i+1)th neighbor of the nth extended target point
	int ne_index;
	for (n = 0; n < N; n++) {
		index = *(target_points + n);
		i = index/((num_z + 1)*(num_y + 1));
		j = (index % ((num_z + 1)*(num_y + 1)))/(num_z + 1);
		k = index % (num_z + 1);
		if ((ne_index = find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*(i + 1)))) >= 0) *(neighbor_material_ir + n + 0*N) = *(material_ir + ne_index); //The neighbor at (x + d, y, z) is part of the target, so copy the neighbor's material index
		if ((ne_index = find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + (num_y + 1)*(i - 1)))) >= 0) *(neighbor_material_ir + n + 1*N) = *(material_ir + ne_index); //The neighbor at (x - d, y, z) is '' '' '' ''...
		if ((ne_index = find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j + 1 + (num_y + 1)*i))) >= 0) *(neighbor_material_ir + n + 2*N) = *(material_ir + ne_index); //The neighbor at (x, y + d, z) is '' '' '' ''...
		if ((ne_index = find_neighbor(initial_reading, N_init, k + (num_z + 1)*(j - 1 + (num_y + 1)*i))) >= 0) *(neighbor_material_ir + n + 3*N) = *(material_ir + ne_index); //The neighbor at (x, y - d, z) is '' '' '' ''...
		if ((ne_index = find_neighbor(initial_reading, N_init, k + 1 + (num_z + 1)*(j + (num_y + 1)*i))) >= 0) *(neighbor_material_ir + n + 4*N) = *(material_ir + ne_index);
		if ((ne_index = find_neighbor(initial_reading, N_init, k - 1 + (num_z + 1)*(j + (num_y + 1)*i))) >= 0) *(neighbor_material_ir + n + 5*N) = *(material_ir + ne_index);
	}
	
//	Splitting the target points into boundary points and interior points: boundary points have at least one neighbor that has a different material composition, interior points have all neighbors being of the same composition as they are
	size_t N_in = 0, N_bound = 0; 
	int *interior_points = (int *)malloc(N*sizeof(int)); //List of interior points
	int *boundary_points = (int *)malloc(N*sizeof(int)); //List of boundary points
	int *material = (int *)malloc(N*sizeof(int)); //Final list of material indices (this will now be sorted so that all interior points come first)
	int *neighbor_material = (int *)malloc(6*N*sizeof(int)); //Final list of the material indices of each point's neighbors (also sorted so that all interior points come first)
	double *Q = (double *)malloc(N*sizeof(double)); //Final list of the heat source at each point (also sorted so that all interior points come first)
	for (n = 0; n < N; n++) {
		m_val = *(material_ir + n);
		if (m_val == *(neighbor_material_ir + n + 0*N) && m_val == *(neighbor_material_ir + n + 1*N) && m_val == *(neighbor_material_ir + n + 2*N) && m_val == *(neighbor_material_ir + n + 3*N) && m_val == *(neighbor_material_ir + n + 4*N) && m_val == *(neighbor_material_ir + n + 5*N)) { //If this point has the same composition as all its neighbors...
			*(interior_points + N_in) = *(target_points + n); //This point is an interior point
			*(material + N_in) = *(material_ir + n);
			*(neighbor_material + N_in + 0*N) = *(neighbor_material_ir + n + 0*N);
			*(neighbor_material + N_in + 1*N) = *(neighbor_material_ir + n + 1*N);
			*(neighbor_material + N_in + 2*N) = *(neighbor_material_ir + n + 2*N);
			*(neighbor_material + N_in + 3*N) = *(neighbor_material_ir + n + 3*N);
			*(neighbor_material + N_in + 4*N) = *(neighbor_material_ir + n + 4*N);
			*(neighbor_material + N_in + 5*N) = *(neighbor_material_ir + n + 5*N);
			*(Q + N_in) = *(Q_ir + n);
			N_in++;
		}
	}
	for (n = 0; n < N; n++) {
		m_val = *(material_ir + n);
		if (m_val != *(neighbor_material_ir + n + 0*N) || m_val != *(neighbor_material_ir + n + 1*N) || m_val != *(neighbor_material_ir + n + 2*N) || m_val != *(neighbor_material_ir + n + 3*N) || m_val != *(neighbor_material_ir + n + 4*N) || m_val != *(neighbor_material_ir + n + 5*N)) { //If this point's composition differs from at least one of its neighbors'...
			*(boundary_points + N_bound) = *(target_points + n); //This point is a boundary point
			*(material + N_in + N_bound) = *(material_ir + n);
			*(neighbor_material + N_in + N_bound + 0*N) = *(neighbor_material_ir + n + 0*N);
			*(neighbor_material + N_in + N_bound + 1*N) = *(neighbor_material_ir + n + 1*N);
			*(neighbor_material + N_in + N_bound + 2*N) = *(neighbor_material_ir + n + 2*N);
			*(neighbor_material + N_in + N_bound + 3*N) = *(neighbor_material_ir + n + 3*N);
			*(neighbor_material + N_in + N_bound + 4*N) = *(neighbor_material_ir + n + 4*N);
			*(neighbor_material + N_in + N_bound + 5*N) = *(neighbor_material_ir + n + 5*N);
			*(Q + N_in + N_bound) = *(Q_ir + n);
			N_bound++;
		}
	}
	


//	Constructing surface denominator matrix and volume "heat source" (see a description of the T-DDA method for an explanation of this part)
	double *D = (double *)calloc(N_bound*N_bound, sizeof(double)); // Denominator matrix
	double *Q_bound = (double *)calloc(N_bound, sizeof(double)); // "Heat source"
	double factor;
	int d_index, d_x, d_y, d_z; // Final total coordinates in lattice units, CAW
	int s_index, s_x, s_y, s_z; // Heat source points 



	for (n = 0; n < N_bound; n++) {
		d_index = *(boundary_points + n);
		d_x = x_coords[d_index/((num_z + 1)*(num_y + 1))];
		d_y = y_coords[(d_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
		d_z = z_coords[d_index % (num_z + 1)];
		
		*(D + n + N_bound*n) = 1.0;
		for (p = 0; p < N_bound; p++) {
			s_index = *(boundary_points + p);
			s_x = x_coords[s_index/((num_z + 1)*(num_y + 1))];
			s_y = y_coords[(s_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
			s_z = z_coords[s_index % (num_z + 1)];
			

			if ( d_x - d < 0 ){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x -d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			}else if ( d_x < 0){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								 (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));

			}else if( d_x + d < 0 ){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));

			}else{

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(D + n + N_bound*p) -= factor/(*kappa)*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								 (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))));

			}
		}
		
		*(Q_bound + n) = *(Q + N_in + n);
		for (l = 0; l < N_in; l++) {
			s_index = *(interior_points + l);
			s_x = x_coords[s_index/((num_z + 1)*(num_y + 1))];
			s_y = y_coords[(s_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
			s_z = z_coords[s_index % (num_z + 1)];
			
			if ( d_x - d < 0 ){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) 
			               + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x -d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) +
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) 
				       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) 
                                       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z + d - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) +
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z - d - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));

			}else if (d_x < 0){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) 
			               + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
				       (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) 
				       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) 
                                       + x_lim*y_lim*abs(d_z - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z + d - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) +
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z - d - s_z)))) - (*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + 
                                         x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));

			
			}else if ( d_x + d < 0){

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((*(G_r + abs(d_x + d - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)) - upper*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) 
                                       + x_lim*y_lim*abs(d_z - s_z)))) - (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
                                                                               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));

			}else{

			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 0*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 0*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 1*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 1*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x - d + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
                                                                               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 2*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 2*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y + d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 3*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 3*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - d - s_y) + x_lim*y_lim*abs(d_z - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 4*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 4*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z + d - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));
			
			factor = 2*(*(kappa + *(material + N_in + n)))*(*(kappa + *(neighbor_material + N_in + n + 5*N)))/(*(kappa + *(material + N_in + n)) + *(kappa + *(neighbor_material + N_in + n + 5*N))) - *kappa;
			*(Q_bound + n) += factor/(*(kappa + *(material + l)))*((lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - d - s_z)))) - 
								               (lower*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z)))))*(*(Q + l));

			}
		}
	}
	
  time_t t2 = time(NULL);
  struct tm tm2 = *localtime(&t2);
  printf("Completed making denom. matrix:\n%02d:%02d:%02d\n\n", tm2.tm_hour, tm2.tm_min, tm2.tm_sec);

double *Q_eff = NULL;


/* ---- Calculating "effective heat", i.e., bound charges ---- */
{
    /* LAPACK uses 32-bit ints */
    if (N_bound > INT_MAX) {
        fprintf(stderr, "System too large for LAPACK int: N_bound=%zu\n", N_bound);
        exit(1);
    }

    int n    = (int)N_bound;
    int nrhs = 1;
    int lda  = n;
    int ldb  = n;
    int info = 0;

    Q_eff = (double *)calloc(N, sizeof *Q_eff);
    if (!Q_eff) { perror("calloc Q_eff"); exit(1); }

    int *ipiv = (int *)malloc(n * sizeof *ipiv);
    if (!ipiv) { perror("malloc ipiv"); exit(1); }

    /* Solve A X = B with A=D, B=Q_bound (in-place overwrite B with solution) */
    dgesv_(&n, &nrhs, D, &lda, ipiv, Q_bound, &ldb, &info);
    free(ipiv);

    if (info != 0) {
        fprintf(stderr, "dgesv_ failed, info = %d\n", info);
        exit(1);
    }

    /* Build Q_eff: interior points scaled by kappa ratio; boundary from solution */
    for (size_t idx = 0; idx < N; idx++) {
        if (idx < N_in) {
            Q_eff[idx] = (*kappa) / (*(kappa + material[idx])) * (Q[idx]);
        } else {
            Q_eff[idx] = Q_bound[idx - N_in];
        }
    }

    time_t t3 = time(NULL);
    struct tm tm3 = *localtime(&t3);
    printf("Completed inverting matrix:\n%02d:%02d:%02d\n\n", tm3.tm_hour, tm3.tm_min, tm3.tm_sec);
}









//	Calculate and write out temperatures
//	Iterate through every point at which we want to determine the temperature

	double *T = (double *)calloc(N + N_out, sizeof(double));

	for (n = 0; n < N + N_out; n++) { 


		if (n < N_in) {
			d_index = *(interior_points + n);

		} else if (n < N) {
			d_index = *(boundary_points + n - N_in);
		} else {
			d_index = *(outside_points + n - N);
		}



		d_x = x_coords[d_index/((num_z + 1)*(num_y + 1))];
		d_y = y_coords[(d_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
		d_z = z_coords[d_index % (num_z + 1)];
		

		 // Iterate through every point at which there is a heat source
		for (l = 0; l < N; l++) {
			

			if (l < N_in) {
				s_index = *(interior_points + l);

			} else {
				s_index = *(boundary_points + l - N_in);
			}



			s_x = x_coords[s_index/((num_z + 1)*(num_y + 1))];
			s_y = y_coords[(s_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
			s_z = z_coords[s_index % (num_z + 1)];
			

		//	This source point's contribution to the temperature at (d_x, d_y, d_z) is its effective heat multiplied by the Green's function value that takes (s_x, s_y, s_z) to (d_x, d_y, d_z)

		// Particle space
			if ( d_x < 0) {
			  *(T + n) += 1/((*kappa)*unit*d)*( 
					  	(*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))) -
					  	 upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))))
			  			* (Q_eff[l]); 
			}

		// In substrate
			else{			  
			  *(T + n) += 1/((*kappa)*unit*d)*(
					  	lower*(*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))))
					  	* (Q_eff[l]); 
			}

		}






		if (n < N) {
			fprintf(out_id, "%d %d %d %lf %lf %lf\n", d_x, d_y, d_z, *(T + n), *(Q + n), *(Q_eff + n));
		} 
		else {
			fprintf(out_id, "%d %d %d %lf %lf %lf\n", d_x, d_y, d_z, *(T + n), 0.0, 0.0);
		}

	}













	fclose(dat_id);
	fclose(par_id);
	fclose(inp_id);	
	fclose(out_id);
	
	time_t t4 = time(NULL); struct tm tm4 = *localtime(&t4);
	printf("Completed temperature calculation:\n%02d:%02d:%02d\n\n", tm4.tm_hour, tm4.tm_min, tm4.tm_sec);

	return 0;
}
