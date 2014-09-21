/**
 * Copyright (C) 2014~2014 by Yichao Yu
 * yyc1992@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

typedef __global const float *gcfloat_p;
typedef __global const unsigned *gcuint_p;

typedef struct {
    unsigned dim_x;
    unsigned dim_y;
    unsigned dim_z;
    unsigned total_dim;

    // direction averaged recoil coupling matrix in x, y, z directions
    gcfloat_p gamma_x; // dim_x * dim_x
    gcfloat_p gamma_y; // dim_y * dim_y
    gcfloat_p gamma_z; // dim_z * dim_z

    // index range in which gamma_x/y/z is not zero
    // minimum index
    gcuint_p gidxmin_x; // dim_x
    gcuint_p gidxmin_y; // dim_y
    gcuint_p gidxmin_z; // dim_z

    // maximum index
    gcuint_p gidxmax_x; // dim_x
    gcuint_p gidxmax_y; // dim_y
    gcuint_p gidxmax_z; // dim_z

    gcfloat_p pump_branch;

    // Rabi frequencies of the Raman transition for unit power
    gcfloat_p omegas;
} ODT;

typedef struct {
    ODT odt;
    float gamma_branch[3][3];
    float gamma_total[3];

    unsigned delta_x;
    unsigned delta_y;
    unsigned delta_z;

    gcfloat_p omega_x;
    gcfloat_p omega_y;
    gcfloat_p omega_z;
} CoolingStep;

typedef struct {
    float h;
    gcfloat_p gamma_total;

    gcuint_p delta_x;
    gcuint_p delta_y;
    gcuint_p delta_z;

    gcuint_p omega_x_offset;
    gcuint_p omega_y_offset;
    gcuint_p omega_z_offset;
} CoolingSequence;

typedef struct {
    gcfloat_p pas;
    gcfloat_p pbs;
    gcfloat_p pcs;
    gcfloat_p qs;
} DensityMatrix;

static inline unsigned
calc_idx_3d(const ODT *odt, unsigned i_x, unsigned i_y, unsigned i_z)
{
    return (i_x * odt->dim_y + i_x) * odt->dim_z + i_z;
}

static float
diff_pump(const CoolingStep *step, const DensityMatrix *mat, float cur_val,
          unsigned i_x, unsigned i_y, unsigned i_z, unsigned i_3d,
          unsigned branch)
{
    const ODT *odt = &step->odt;
    float res = 0;
    float *gamma_branch = &step->gamma_branch[branch][0];
    // Also try using the other dimension
    gcfloat_p gamma_xs = &odt->gamma_x[i_x];
    gcfloat_p gamma_ys = &odt->gamma_y[i_y];
    gcfloat_p gamma_zs = &odt->gamma_z[i_z];
    for (unsigned i = odt->gidxmin_x;i < odt->gidxmax_x;i++) {
        for (unsigned j = odt->gidxmin_y;j < odt->gidxmax_y;j++) {
            for (unsigned k = odt->gidxmin_z;k < odt->gidxmax_z;k++) {
                unsigned ijk = calc_idx_3d(odt, i, j, k);
                unsigned idx = ijk * odt->total_dim + i_3d;
                float pa = mat->pas[idx];
                float pb = mat->pbs[idx];
                float pc = mat->pcs[idx];

                float gamma_x = gamma_xs[i * odt->dim_x];
                float gamma_y = gamma_ys[j * odt->dim_y];
                float gamma_z = gamma_zs[k * odt->dim_z];

                res += (gamma_x * gamma_y * gamma_z
                        * (gamma_branch[0] * pa + gamma_branch[1] * pb
                           + gamma_branch[2] * pc));
            }
        }
    }
    return res - cur_val * step->gamma_total[branch];
}

static float
diff_pa(const CoolingStep *step, const DensityMatrix *mat, float cur_val,
        unsigned i_x, unsigned i_y, unsigned i_z, unsigned i_3d)
{
    const ODT *odt = &step->odt;
    float res = diff_pump(step, mat, cur_val, i_x, i_y, i_z, i_3d, 0);
    return res + 2 * omega_x[i_x] * omega_y[i_y] * omega_z[i_z] * mat->qs[i_3d];
}

static float
diff_pb(const CoolingStep *step, const DensityMatrix *mat, float cur_val,
        unsigned i_x, unsigned i_y, unsigned i_z, unsigned i_3d)
{
    float res = diff_pump(step, mat, cur_val, i_x, i_y, i_z, i_3d, 1);
    if (i_x < step->delta_x || i_y < step->delta_y || i_z < step->delta_z) {
        return res;
    }
    const ODT *odt = &step->odt;
    unsigned i_x_new = i_x - step->delta_x;
    unsigned i_y_new = i_y - step->delta_y;
    unsigned i_z_new = i_z - step->delta_z;
    return (res + 2 * omega_x[i_x_new] * omega_y[i_y_new] * omega_z[i_z_new] *
            mat->qs[calc_idx_3d(odt, i_x_new, i_y_new, i_z_new)]);
}

static float
diff_pc(const CoolingStep *step, const DensityMatrix *mat, float cur_val,
        unsigned i_x, unsigned i_y, unsigned i_z, unsigned i_3d)
{
    return diff_pump(step, mat, cur_val, i_x, i_y, i_z, i_3d, 2);
}

static float
diff_q(const CoolingStep *step, const DensityMatrix *mat, float cur_val,
       unsigned i_x, unsigned i_y, unsigned i_z, unsigned i_3d)
{
    const ODT *odt = &step->odt;
    unsigned i_x_new = i_x + step->delta_x;
    unsigned i_y_new = i_y + step->delta_y;
    unsigned i_z_new = i_z + step->delta_z;
    if (i_x_new >= odt->dim_x || i_y_new >= odt->dim_y ||
        i_z_new >= odt->dim_z) {
        return 0;
    }
    float res = -cur_val / 2 * (step->gamma_total[0] + step->gamma_total[1]);
    float omega = omega_x[i_x] * omega_y[i_y] * omega_z[i_z];
    return omega * (mat->pbs[calc_idx_3d(odt, i_x_new, i_y_new, i_z_new)]
                    - mat->pas[i_3d]);
}

static float
diff_mat(const CoolingStep *step, const DensityMatrix *mat, float cur_val,
         unsigned glob_idx)
{
    const ODT *odt = &step->odt;
    unsigned cls = glob_idx / odt->total_dim;
    unsigned i_3d = glob_idx % odt->total_dim;
    unsigned dim_xy = odt->dim_x * odt->dim_y;
    unsigned xy = i_3d / dim_xy;
    unsigned i_z = i_3d % dim_xy;
    unsigned i_x = xy / odt->dim_y;
    unsigned i_y = xy % odt->dim_y;
    switch (cls) {
    case 0:
        return diff_pa(step, mat, cur_val, i_x, i_y, i_z, i_3d);
    case 1:
        return diff_pb(step, mat, cur_val, i_x, i_y, i_z, i_3d);
    case 2:
        return diff_pc(step, mat, cur_val, i_x, i_y, i_z, i_3d);
    default:
        return diff_q(step, mat, cur_val, i_x, i_y, i_z, i_3d);
    }
}

static void
fill_step(CoolingStep *step, const CoolingSequence *seq, float t)
{
    unsigned idx_t = t / seq->h;
    gcfloat_p gamma_total = seq->gamma_total + idx_t * 3;
    for (unsigned i = 0;i < 3;i++) {
        float gamma = gamma_total[i];
        step->gamma_total[i] = gamma;
        for (unsigned j = 0;j < 3;j++) {
            step->gamma_branch[i][j] = step->odt.pump_branch[i * 3 + j] * gamma;
        }
    }
    step->delta_x = seq->delta_x[idx_t];
    step->delta_y = seq->delta_y[idx_t];
    step->delta_z = seq->delta_z[idx_t];

    step->omega_x = step->odt.omega + seq->omega_x_offset[idx_t];
    step->omega_y = step->odt.omega + seq->omega_y_offset[idx_t];
    step->omega_z = step->odt.omega + seq->omega_z_offset[idx_t];
}

static float
calc_sbcooling_diff(float t, gcfloat_p mat_in, unsigned glob_idx, float cur_val,
                    unsigned dim_x, unsigned dim_y, unsigned dim_z,
                    gcfloat_p gamma_xyz, gcuint_p gidx_minmax_xyz,
                    gcfloat_p pump_branch, gcfloat_p omegas,
                    float h_t, gcfloat_p gamma_total, unsigned seq_len,
                    gcuint_p delta_xyz, gcuint_p omega_xyz_offset)
{
    CoolingStep step;
    ODT *odt = &step.odt;

    odt->dim_x = dim_x;
    odt->dim_y = dim_y;
    odt->dim_z = dim_z;
    odt->total_dim = dim_x * dim_y * dim_z;

    odt->gamma_x = gamma_xyz;
    odt->gamma_y = gamma_xyz + dim_x * dim_x;
    odt->gamma_z = odt->gamma_y + dim_y * dim_y;

    odt->gidxmin_x = gidx_minmax_xyz;
    odt->gidxmin_y = odt->gidxmin_x + dim_x;
    odt->gidxmin_z = odt->gidxmin_y + dim_y;

    odt->gidxmax_x = odt->gidxmin_z + dim_z;
    odt->gidxmax_y = odt->gidxmax_x + dim_x;
    odt->gidxmax_z = odt->gidxmax_y + dim_y;

    odt->pump_branch = pump_branch;
    odt->omegas = omegas;

    CoolingSequence seq;
    seq.h = h_t;
    seq.gamma_total = gamma_total;

    seq.delta_x = delta_xyz;
    seq.delta_y = delta_xyz + seq_len;
    seq.delta_z = seq.delta_y + seq_len;

    seq.omega_x_offset = omega_xyz_offset;
    seq.omega_y_offset = omega_xyz_offset + seq_len;
    seq.omega_z_offset = seq.omega_y_offset + seq_len;

    fill_step(&step, &seq, t);

    DensityMatrix mat;
    mat.pas = mat_in;
    mat.pbs = mat_in + odt->total_dim;
    mat.pcs = mat.pbs + odt->total_dim;
    mat.qs = mat.pcs + odt->total_dim;

    return diff_mat(&step, &mat, cur_val, glob_idx);
}
