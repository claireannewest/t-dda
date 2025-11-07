import numpy as np
import time
starting_time = time.time()


class SteadyState_Temp:
    def __init__(self,
                 green_tab,
                 var_par,
                 optical_data):
        """Defines the different system parameters.

        Keyword arguments:
        """
        self.green_tab = green_tab
        self.var_par = var_par
        self.optical_data = optical_data

    def read_varpar(self):
        varpar = open(self.var_par, "r")
        file = varpar.readlines()
        num_k = file[0].split(' ')[1]
        kap = np.zeros(len(num_k) + 2)
        for idx, kappai in enumerate(kap):
            kap[idx] = file[idx + 1].split(' ')[1]
        wave = float(file[idx + 2].split(' ')[1])
        n_back = float(file[idx + 3].split(' ')[1])
        I0 = float(file[idx + 4].split(' ')[1])
        unit = float(file[idx + 5].split(' ')[1])
        d = float(file[idx + 7].split(' ')[1])
        x_min = float(file[idx + 8].split(' ')[1])
        x_max = float(file[idx + 9].split(' ')[1])
        y_min = float(file[idx + 10].split(' ')[1])
        y_max = float(file[idx + 11].split(' ')[1])
        z_min = float(file[idx + 12].split(' ')[1])
        z_max = float(file[idx + 13].split(' ')[1])
        varcoords = [unit, d, x_min, x_max, y_min, y_max, z_min, z_max]
        plane = file[idx + 15].split(' ')[1]
        input_mode = file[idx + 16].split(' ')[1]
        tot_pnts = file[idx + 17].split(' ')[1]
        return kap, wave, n_back, I0, varcoords, plane, input_mode, tot_pnts

    def define_coords(self):
        #  Create arrays of all possible x, y, and z coordinates
        #  Will narrow down later
        _, _, _, _, varcoords, _, _, _ = self.read_varpar()
        _, d, x_min, x_max, y_min, y_max, z_min, z_max = varcoords
        num_x = int((x_max - x_min) / d + 1)
        num_y = int((y_max - y_min) / d + 1)
        num_z = int((z_max - z_min) / d + 1)
        x_coords = np.arange(x_min, x_max + d, d)
        y_coords = np.arange(y_min, y_max + d, d)
        z_coords = np.arange(z_min, z_max + d, d)
        return num_x, num_y, num_z, x_coords, y_coords, z_coords

    def readingreen(self):
        num_x, num_y, num_z, _, _, _ = self.define_coords()
        file = open(self.green_tab, "r")
        green = file.readlines()
        _, xlim, ylim, zlim = green[0].rstrip().split(' ')
        if num_x > int(xlim) or num_y > int(ylim) or num_z > int(zlim):
            print('ERROR: Shape extends greater than Lattice Green grid.\n')
        G_coords = np.zeros((len(green), 3))
        G_val = np.zeros((len(green), 1))
        idx = 0
        with open(self.green_tab) as f:
            for line in f:
                if idx == 0:
                    idx = 1
                    continue
                Gx, Gy, Gz, G_val[idx] = line.rstrip().split(' ')
                G_coords[idx, 0] = Gx
                G_coords[idx, 1] = Gy
                G_coords[idx, 2] = Gz
                idx += 1

    def readinoptics(self):
        # x_o, y_o, z_o are the coordinates from the optical calculation
        # N_o is number of coordinates in optical calculation
        _, wave, n_back, I0, varcoords, plane, input_mode, _ = self.read_varpar()
        unit, d, _, _, _, _, _, _ = varcoords
        if input_mode.rstrip() == str(1):
            file = open(self.optical_data, "r")
            opt_data = file.readlines()
            N_init = len(opt_data)
            Q_ir = np.zeros(N_init)
            material_ir = np.zeros(N_init)
            coords_init = np.zeros((N_init, 3))
            with open(self.optical_data) as f:
                idx = 0
                for line in f:
                    raw_list = line.strip().split(' ')
                    fin_list = ' '.join(raw_list).split()
                    x_f, y_f, z_f, m_val, E_xr, E_xi, E_yr, E_yi, E_zr, E_zi, P_xr, P_xi, P_yr, P_yi, P_zr, P_zi = fin_list
                    E = [float(E_xr), float(E_xi),
                         float(E_yr), float(E_yi),
                         float(E_zr), float(E_zi), ]
                    P = [float(P_xr), float(P_xi),
                         float(P_yr), float(P_yi),
                         float(P_zr), float(P_zi), ]
                    k = 2 * np.pi * n_back / wave
                    extinc = 4 * np.pi * k * I0 * (unit * d)**3 * (E[0]*P[1] - E[1]*P[0] + E[2]*P[3] - E[3]*P[2] + E[4]*P[5] - E[5]*P[4]) * (1e-18)
                    scatt = 8 / 3 * np.pi * k**4 * I0*(unit*d)**6*(P[0]*P[1] + P[2]*P[3] + P[4]*P[5])*(1e-45)
                    Q_ir[idx] = extinc - scatt
                    material_ir[idx] = m_val
                    coords_init[idx, 0] = x_f
                    coords_init[idx, 1] = y_f
                    coords_init[idx, 2] = z_f
                    idx += 1
        # if input_mode.rstrip() == str(2):
        #     with open(self.optical_data) as f:
        #         for line in f:
        #             raw_list = line.strip().split(' ')
        #             fin_list = ' '.join(raw_list).split()
        #             x_f, y_f, z_f, Q_byhand = fin_list
                    # material_ir[idx] = 1
                    # coords_init[idx, 0] = x_f
                    # coords_init[idx, 1] = y_f
                    # coords_init[idx, 2] = z_f
                    # Q_ir[idx] = Q_byhand
        return N_init, Q_ir, material_ir, coords_init

    def find_neighbor_ij(self, coordi, neigh, mat_ir):
        coord_neigh = coordi + neigh
        idx_neigh = np.where((coordi[0] == coordi[0] + neigh[0])
                             & (coordi[1] == coordi[1] + neigh[1])
                             & (coordi[2] == coordi[2] + neigh[2]))[0]
        if len(idx_neigh) == 1:
            mat_neigh = mat_ir[idx_neigh]
        else:
            mat_neigh = 0
        return coord_neigh, mat_neigh

    def find_neighbors_all_together(self):
        N_init, Q_ir, material_ir, coords_init = self.readinoptics()
        dicts_for_N_init = [{}] * N_init
        for i in range(N_init):
            coordi = coords_init[i, :]
            pnti = dicts_for_N_init[i]
            pnti["coord"] = coordi
            pnti["mat"] = material_ir[i]
            coord_ipx, mat_ipx = self.find_neighbor_ij(coordi=coordi,
                                                       neigh=[1, 0, 0],
                                                       mat_ir=material_ir,
                                                       )
            pnti["coord_px"] = coord_ipx
            pnti["mat_px"] = mat_ipx
            coord_imx, mat_imx = self.find_neighbor_ij(coordi=coordi,
                                                       neigh=[-1, 0, 0],
                                                       mat_ir=material_ir,
                                                       )
            pnti["coord_mx"] = coord_imx
            pnti["mat_mx"] = mat_imx
            coord_ipy, mat_ipy = self.find_neighbor_ij(coordi=coordi,
                                                       neigh=[0, 1, 0],
                                                       mat_ir=material_ir,
                                                       )
            pnti["coord_py"] = coord_ipy
            pnti["mat_py"] = mat_ipy
            coord_imy, mat_imy = self.find_neighbor_ij(coordi=coordi,
                                                       neigh=[0, -1, 0],
                                                       mat_ir=material_ir,
                                                       )
            pnti["coord_my"] = coord_imy
            pnti["mat_my"] = mat_imy
            coord_ipz, mat_ipz = self.find_neighbor_ij(coordi=coordi,
                                                       neigh=[0, 0, 1],
                                                       mat_ir=material_ir,
                                                       )
            pnti["coord_pz"] = coord_ipz
            pnti["mat_pz"] = mat_ipz
            coord_imz, mat_imz = self.find_neighbor_ij(coordi=coordi,
                                                       neigh=[0, 0, -1],
                                                       mat_ir=material_ir,
                                                       )
            pnti["coord_mz"] = coord_imz
            pnti["mat_mz"] = mat_imz
        return N_init, dicts_for_N_init

    def define_int_and_boundaries(self):
        # If all of your neighbors are in dda, you're an interior point.
        N_init, dicts_for_N_init = self.find_neighbors_all_together()
        for i in range(N_init):
            pnti = dicts_for_N_init[i]
            print(pnti)
            if pnti['mat_px'] == 0 :
                print('here')
# 			if sum([x_pd, x_md, y_pd, y_md, z_pd, z_md]) == 6:
# 				interior[idx_int, :] = val
# 				idx_int = idx_int + 1

# 			# If any of your neighbors are NOT in the optical list, you're a boundary point.
# 			# And that neighbor that isn't in optical list needs to be added to the extended target.
# 			if sum([x_pd, x_md, y_pd, y_md, z_pd, z_md]) < 6:
# 				boundary[idx_boun, :] = val
# 				idx_boun = idx_boun + 1
# 				# Whichever point was False, add it to extended_target_pnts.
# 				if x_pd == False: extended.append(val + [self.d, 0, 0])
# 				if x_md == False: extended.append(val + [-self.d, 0, 0])
# 				if y_pd == False: extended.append(val + [0, self.d, 0])
# 				if y_md == False: extended.append(val + [0, -self.d, 0])
# 				if z_pd == False: extended.append(val + [0, 0, self.d])
# 				if z_md == False: extended.append(val + [0, 0, -self.d])

# 			## Need to make a loop about the extra temperature points!


# 		# Trim up the arrays
# 		target_pnts_interior = interior[:idx_int,:]
# 		target_pnts_boundary = boundary[:idx_boun,:]
# 		target_pnts_extended = np.vstack(extended)


#
#
#
#
#
#
#
#
#
#
#
#
#
#

#
#
#
#
#
#
    def findD(self,
               ii,  # which side, 0 -> 5
               term1, # lower or upper
               term2, # lower or upper
               ):
        kap_at_mat_Nin_plus_n = 1
        kap_at_neighmat_Nin_plus_n_iiN = 1
        numer = 2 * kap_at_mat_Nin_plus_n * kap_at_neighmat_Nin_plus_n_iiN
        denom = kap_at_mat_Nin_plus_n + kap_at_neighmat_Nin_plus_n_iiN
        factor = numer / denom - self.kappa[0]

        Gr = 1  # Gr(dx-sx, dy-sy, dz-sz)
        Grp = 1  # Gr(dx+sx, dy-sy, dz-sz)
        Gr_1 = 1  # Gr(dx+sx+d, dy-sy, dz-sz)

        Gr_px = 1  # Gr(dx-sx+d, dy-sy, dz-sz)
        Gr_mx = 1  # Gr(dx-sx-d, dy-sy, dz-sz)
        Gr_py = 1  # Gr(dx-sx, dy-sy+d, dz-sz)
        Gr_my = 1  # Gr(dx-sx, dy-sy-d, dz-sz)
        Gr_pz = 1  # Gr(dx-sx, dy-sy, dz-sz+d)
        Gr_mz = 1  # Gr(dx-sx, dy-sy, dz-sz-d)
        if ii == 0:
            greens = (Gr_px - term1 * Gr_1) - (Gr - term2 * Grp)
        if ii == 1:
            greens = (Gr_mx - term1 * Gr_1) - (Gr - term2 * Grp)
        if ii == 2:
            greens = (Gr_py - term1 * Gr_1) - (Gr - term2 * Grp)
        if ii == 3:
            greens = (Gr_my - term1 * Gr_1) - (Gr - term2 * Grp)
        if ii == 4:
            greens = (Gr_pz - term1 * Gr_1) - (Gr - term2 * Grp)
        if ii == 5:
            greens = (Gr_mz - term1 * Gr_1) - (Gr - term2 * Grp)
        D = factor / self.kappa[0] * greens
        return D

    def findQbound(self,
                   ii,  # which side, 0 -> 5
                   term1,  # lower or upper
                   term2,  # lower or upper
                   ):
        kap_at_mat_Nin_plus_n = 1
        kap_at_neighmat_Nin_plus_n_iiN = 1
        numer = 2 * kap_at_mat_Nin_plus_n * kap_at_neighmat_Nin_plus_n_iiN
        denom = kap_at_mat_Nin_plus_n + kap_at_neighmat_Nin_plus_n_iiN
        factor = numer / denom - self.kappa[0]
        if ii == 0:
            greens = (Gr_px - term1 * Gr_1) - (Gr - term2 * Grp)

        Qbound = factor / kappa_at_ll * greens * Q_atll
        return Qbound

    def construct_D(self,
                    N_bound,
                    boundary_pnts,
                    x_coords,
                    y_coords,
                    z_coords,
                    N_in,
                    Q,
                    interior_points,
                    ):
        # Constructing surface denominator matrix and volume "heat source"
        D = np.zeros(N_bound * N_bound)  # denominator matrix
        Q_bound = np.zeros(N_bound * N_bound)  # "Heat source"
        num_k = len(self.kappa)
        upper = ((self.kappa[num_k + 1] - self.kappa[0])
                 / (self.kappa[num_k + 1] + self.kappa[0]))
        lower = 2 * self.kappa[0] / (self.kappa[num_k + 1] + self.kappa[0])
        # d_x, d_y, d_z; Final total coordinates in lattice units, CAW
        # s_x, s_y, s_z; Heat source points

        for n in range(0, N_bound):
            d_idx = boundary_pnts[n]
            d_x = x_coords[d_idx]
            # d_y = y_coords[d_idx]
            # d_z = z_coords[d_idx]
            D[n + N_bound] = 1.0
            for p in range(0, N_bound):
                s_idx = boundary_pnts[p]
                s_x = x_coords[s_idx]
                s_y = y_coords[s_idx]
                s_z = z_coords[s_idx]
                if d_x - self.d < 0:
                    D[n + N_bound * p] -= self.findD(0, upper, upper)
                    D[n + N_bound * p] -= self.findD(1, upper, upper)
                    D[n + N_bound * p] -= self.findD(2, upper, upper)
                    D[n + N_bound * p] -= self.findD(3, upper, upper)
                    D[n + N_bound * p] -= self.findD(4, upper, upper)
                    D[n + N_bound * p] -= self.findD(5, upper, upper)
                elif d_x < 0:
                    D[n + N_bound * p] -= self.findD(0, upper, upper)
                    D[n + N_bound * p] -= self.findD(1, upper, upper)
                    D[n + N_bound * p] -= self.findD(2, upper, upper)
                    D[n + N_bound * p] -= self.findD(3, upper, upper)
                    D[n + N_bound * p] -= self.findD(4, upper, upper)
                    D[n + N_bound * p] -= self.findD(5, upper, upper)
                elif d_x + self.d < 0:
                    D[n + N_bound * p] -= self.findD(0, upper, lower)
                    D[n + N_bound * p] -= self.findD(1, upper, lower)
                    D[n + N_bound * p] -= self.findD(2, upper, lower)
                    D[n + N_bound * p] -= self.findD(3, upper, lower)
                    D[n + N_bound * p] -= self.findD(4, upper, lower)
                    D[n + N_bound * p] -= self.findD(5, upper, lower)
                else:
                    D[n + N_bound * p] -= self.findD(0, lower, lower)
                    D[n + N_bound * p] -= self.findD(1, lower, lower)
                    D[n + N_bound * p] -= self.findD(2, lower, lower)
                    D[n + N_bound * p] -= self.findD(3, lower, lower)
                    D[n + N_bound * p] -= self.findD(4, lower, lower)
                    D[n + N_bound * p] -= self.findD(5, lower, lower)
            Q_bound[n] = Q[N_in + n]
            for ll in range(0, N_in):
                s_idx = interior_points[ll]
                # s_x = x_coords[s_idx]
                # s_y = y_coords[s_idx]
                # s_z = z_coords[s_idx]
                if d_x - self.d < 0:
                    Q_bound[n] += self.findQbound(0, upper, upper)
                    Q_bound[n] += self.findQbound(1, upper, upper)
                    Q_bound[n] += self.findQbound(2, upper, upper)
                    Q_bound[n] += self.findQbound(3, upper, upper)
                    Q_bound[n] += self.findQbound(4, upper, upper)
                    Q_bound[n] += self.findQbound(5, upper, upper)
                elif d_x < 0:
                    Q_bound[n] += self.findQbound(0, upper, upper)
                    Q_bound[n] += self.findQbound(1, upper, upper)
                    Q_bound[n] += self.findQbound(2, upper, upper)
                    Q_bound[n] += self.findQbound(3, upper, upper)
                    Q_bound[n] += self.findQbound(4, upper, upper)
                    Q_bound[n] += self.findQbound(5, upper, upper)
                elif d_x + self.d < 0:
                    Q_bound[n] += self.findQbound(0, upper, lower)
                    Q_bound[n] += self.findQbound(1, upper, lower)
                    Q_bound[n] += self.findQbound(2, upper, lower)
                    Q_bound[n] += self.findQbound(3, upper, lower)
                    Q_bound[n] += self.findQbound(4, upper, lower)
                    Q_bound[n] += self.findQbound(5, upper, lower)
                else:
                    Q_bound[n] += self.findQbound(0, lower, lower)
                    Q_bound[n] += self.findQbound(1, lower, lower)
                    Q_bound[n] += self.findQbound(2, lower, lower)
                    Q_bound[n] += self.findQbound(3, lower, lower)
                    Q_bound[n] += self.findQbound(4, lower, lower)
                    Q_bound[n] += self.findQbound(5, lower, lower)

    def calc_Q(self, D, Q, Q_bound):
        # Calculating "effective heat"
        # AX = B ....> A = D, B = Q_bound.
        # sgesv_(&N_bound, &n_rhs, D, &N_bound, ipiv, Q_bound, &N_bound, &info);

        for n in range(0, N):
            # If the current point is an interior point
            if n < N_in:
                # Heat source at this pnt reduced by factor of target conduct
                kap_at_n = 1
                Q_eff[n] = self.kappa[0] / (kap_at_n) * Q[n]
            else:
                # Heat source is D*Q_eff = Q_bound
                Q_eff[n] = Q_bound[n - N_in]
        return Q_eff

    def calcd_idx_temp(self,
                       N,
                       N_out,
                       N_in,
                       x_coords,
                       y_coords,
                       z_coords,
                       interior_pnts,
                       boundary_pnts,
                       outside_pnts,
                       Q_eff
                       ):
        T = np.zeros(N + N_out)
        # Iterate through every point at which we want the temperature
        for n in range(0, N + N_out):
            # Find d_x, d_y, d_z coords.
            if n < N_in:
                d_idx = interior_pnts[n]
            elif n < N:
                d_idx = boundary_pnts[n - N_in]
            else:
                d_idx = outside_pnts[n - N]
            d_x = x_coords[d_idx]
            # d_y = y_coords[d_idx]
            # d_z = z_coords[d_idx]

            # Iterate through every point at which there is a heat source
            for ll in range(0, N):
                if ll < N_in:
                    s_idx = interior_pnts[ll]
                else:
                    s_idx = boundary_pnts[ll - N_in]
                s_x = x_coords[s_idx]
                s_y = y_coords[s_idx]
                s_z = z_coords[s_idx]

                # This source point's contribution to the temperature at
                # (d_x, d_y, d_z) is its effective heat multiplied by the
                # Green's function value that takes (s_x, s_y, s_z) to
                # (d_x, d_y, d_z)
                G_d_s = 1  # G_r(dx-sx, dy-sy, dz-sz)
                G_d_px_s = 1  # G_r(dx+sx, dy-sy, dz-sz)

                # In particle space
                if d_x < 0:
                    T[n] += (1 / (kap[0] * self.unit * self.d)
                             * (G_d_s - self.upper * G_d_px_s) * Q_eff[ll])

                # In substrate space
                else:
                    T[n] += (1 / (kap[0] * self.unit * self.d)
                             * (self.lower * G_d_s) * Q_eff[ll])

    def write_temp(self, N, N_out, d, T, Q, Q_eff):
        d_x, d_y, d_z = d
        fout = open('temp.out', 'w')
        for n in range(0, N + N_out):
            if n < N:
                fout.write(d_x, d_y, d_z, T, Q, Q_eff)
            else:
                fout.write(d_x, d_y, d_z, T, 0, 0)
        fout.close()


path2green = '../../lattice_greenfunction/Green_grid100.txt'
sst = SteadyState_Temp(green_tab=path2green,
                       var_par='var.par',
                       optical_data='tdda_input_test',
                       )
sst.define_int_and_boundaries()
# sst.find_target_and_extended_target_points()

# ending_time = time.time()
# print('\nCompute time:', ending_time - starting_time)




