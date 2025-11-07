#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[]) {
	char *dat_name = argv[1];
	FILE *dat_id; 
	dat_id = fopen(dat_name, "r");

	int i, j, k;
  	int x_lim = 0, y_lim = 0, z_lim = 0;



	fscanf(dat_id, "%d %d %d\n", &x_lim, &y_lim, &z_lim);

	float *G_r = (float *)malloc(x_lim*y_lim*z_lim*sizeof(float));
	float dat_val = 0.0;


	while (!feof(dat_id)) {
		fscanf(dat_id, "%d %d %d %f\n", &i, &j, &k, &dat_val);
		*(G_r + i + x_lim*j + x_lim*y_lim*k) = dat_val;
	}



	int N, N_out, N_in; 
	N_in = 8;
	N = 10;
	N_out = 5;

	float *T = (float *)calloc(N + N_out, sizeof(float));
	int d_index, d_x, d_y, d_z; // Final total coordinates in lattice units, CAW
	int n; 

	for (n = 0; n < N + N_out; n++) { 


		if (n < N_in) {
			d_index = 1; //*(interior_points + n);

		} else if (n < N) {
			d_index = 1; //*(boundary_points + n - N_in);
		} else {
			d_index = 1; //*(outside_points + n - N);
		}

		// d_x = x_coords[d_index/((num_z + 1)*(num_y + 1))];
		// d_y = y_coords[(d_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
		// d_z = z_coords[d_index % (num_z + 1)];
		

		//  // Iterate through every point at which there is a heat source
		// for (l = 0; l < N; l++) {
			

		// 	if (l < N_in) {
		// 		s_index = *(interior_points + l);

		// 	} else {
		// 		s_index = *(boundary_points + l - N_in);
		// 	}



		// 	s_x = x_coords[s_index/((num_z + 1)*(num_y + 1))];
		// 	s_y = y_coords[(s_index % ((num_z + 1)*(num_y + 1)))/(num_z + 1)];
		// 	s_z = z_coords[s_index % (num_z + 1)];
			

		// //	This source point's contribution to the temperature at (d_x, d_y, d_z) is its effective heat multiplied by the Green's function value that takes (s_x, s_y, s_z) to (d_x, d_y, d_z)

		// // Particle space
		// 	if ( d_x < 0) {
		// 	  *(T + n) += 1/((*kappa)*unit*d)*( 
		// 			  	(*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))) -
		// 			  	 upper*(*(G_r + abs(d_x + s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))))*
		// 	  			(*(Q_eff + l)); 
		// 	}

		// // In substrate
		// 	else{			  
		// 	  *(T + n) += 1/((*kappa)*unit*d)*(
		// 			  	lower*(*(G_r + abs(d_x - s_x) + x_lim*abs(d_y - s_y) + x_lim*y_lim*abs(d_z - s_z))))*
		// 			  	(*(Q_eff + l)); 
		// 	}

		// }






		// if (n < N) {
		// 	fprintf(out_id, "%d %d %d %f %f %f\n", d_x, d_y, d_z, *(T + n), *(Q + n), *(Q_eff + n));
		// } 
		// else {
		// 	fprintf(out_id, "%d %d %d %f %f %f\n", d_x, d_y, d_z, *(T + n), 0.0, 0.0);
		// }



	}
















	return 0;
}
