#!/usr/bin/env python

import os.path as _path
import numpy as np
import pyopencl as cl
from cffi import FFI

from pyscical.ocl.ode import ElwiseOdeSolver
from pyscical.ocl.utils import CLArg
from pyscical.utils import cffi_ptr
from pyscical.atomic import sideband_strength, sideband_scatter_strength

_ffi = FFI()
_ffi.cdef('''
void _fill_gidx_minmax(unsigned dim, float *gamma,
                       float *min_out, float *max_out);''')
_lib = _ffi.verify('''
#include <math.h>

static const float gamma_thresh = 1e-10;

static void
_fill_gidx_minmax(unsigned dim, float *gamma, unsigned *min_out,
                  unsigned *max_out)
{
#pragma omp parallel for
    for (unsigned i = 0;i < dim;i++) {
        float *row = gamma + i * dim;
        int j = 0;
        for (;(unsigned)j < dim;j++) {
            if (fabsf(row[j]) > gamma_thresh) {
                break;
            }
        }
        min_out[i] = j;

        j = dim - 1;
        for (;j >= 0;j--) {
            if (fabsf(row[j]) > gamma_thresh) {
                break;
            }
        }
        max_out[i] = j + 1;
    }
}
''', extra_compile_args=['-w', '-fopenmp', '-O2', '-std=gnu99'],
extra_link_args=['-fopenmp'])

def evolve_sideband(ctx, queue, gamma_x, gamma_y, gamma_z, pump_branch,
                    omegas_x, omegas_y, omegas_z, h_t, gamma_totals,
                    delta_xyz, omega_xyz):
    dim_x, d = gamma_x.shape
    if dim_x != d:
        raise ValueError("gamma_x is not a square matrix.")
    if gamma_x.dtype != np.float32:
        raise TypeError("The type of gamma_x should be float32.")

    dim_y, d = gamma_y.shape
    if dim_y != d:
        raise ValueError("gamma_y is not a square matrix.")
    if gamma_y.dtype != np.float32:
        raise TypeError("The type of gamma_y should be float32.")

    dim_z, d = gamma_z.shape
    if dim_z != d:
        raise ValueError("gamma_z is not a square matrix.")
    if gamma_z.dtype != np.float32:
        raise TypeError("The type of gamma_z should be float32.")

    mf = cl.mem_flags
    gamma_xyz = cl.Buffer(ctx, mf.READ_ONLY,
                          (dim_x**2 + dim_y**2 + dim_z**2) * 4)
    events = []
    events.append(cl.enqueue_copy(queue, gamma_xyz, gamma_x,
                                  device_offset=0, is_blocking=False))
    events.append(cl.enqueue_copy(queue, gamma_xyz, gamma_y,
                                  device_offset=dim_x**2 * 4,
                                  is_blocking=False))
    events.append(cl.enqueue_copy(queue, gamma_xyz, gamma_z,
                                  device_offset=(dim_x**2 + dim_y**2) * 4,
                                  is_blocking=False))

    dev = queue.device
    src = """
    #include <sideband.cl>
    """
    extra_args = (CLArg('dim_x', 'unsigned'),
                  CLArg('dim_y', 'unsigned'),
                  CLArg('dim_z', 'unsigned'),
                  CLArg('gamma_xyz', 'gcfloat_p'),
                  CLArg('gidx_minman_xyz', 'gcuint_p'),
                  CLArg('pump_branch', 'gcfloat_p'),
                  CLArg('omegas', 'gcfloat_p'),
                  CLArg('h_t', 'float'),
                  CLArg('seq_len', 'unsigned'),
                  CLArg('gamma_total', 'gcfloat_p'),
                  CLArg('delta_xyz', 'gcuint_p'),
                  CLArg('omega_xyz_offset', 'gcuint_p'))
    solver = ElwiseOdeSolver(ctx, dev, src, "calc_sbcooling_diff",
                             extra_args=extra_args,
                             options=['-I', _path.dirname(__file__)])

    # Probably faster if doing is on GPU
    gidx_minmax_xyz_cpu = np.empty((dim_x + dim_y + dim_z) * 2, np.uint32)

    gidxmin_x = gidx_minmax_xyz_cpu[:dim_x]
    gidxmin_y = gidx_minmax_xyz_cpu[dim_x:dim_x + dim_y]
    gidxmin_z = gidx_minmax_xyz_cpu[dim_x + dim_y:dim_x + dim_y + dim_z]

    gidxmax_x = gidx_minmax_xyz_cpu[dim_x + dim_y + dim_z:
                                    dim_x * 2 + dim_y + dim_z]
    gidxmax_y = gidx_minmax_xyz_cpu[dim_x * 2 + dim_y + dim_z:
                                dim_x * 2 + dim_y * 2 + dim_z]
    gidxmax_z = gidx_minmax_xyz_cpu[dim_x * 2 + dim_y * 2 + dim_z:]
    _lib._fill_gidx_minmax(dim_x, cffi_ptr(gamma_x, _ffi)[0],
                           cffi_ptr(gidxmin_x, _ffi, writable=True)[0],
                           cffi_ptr(gidxmax_x, _ffi, writable=True)[0])
    _lib._fill_gidx_minmax(dim_y, cffi_ptr(gamma_y, _ffi)[0],
                           cffi_ptr(gidxmin_y, _ffi, writable=True)[0],
                           cffi_ptr(gidxmax_y, _ffi, writable=True)[0])
    _lib._fill_gidx_minmax(dim_z, cffi_ptr(gamma_z, _ffi)[0],
                           cffi_ptr(gidxmin_z, _ffi, writable=True)[0],
                           cffi_ptr(gidxmax_z, _ffi, writable=True)[0])

    gidx_minmax_xyz = cl.Buffer(ctx, mf.READ_ONLY, (dim_x + dim_y + dim_z) * 8)
    events.append(cl.enqueue_copy(queue, gidx_minmax_xyz, gidx_minmax_xyz_cpu,
                                  device_offset=0, is_blocking=False))

    pump_branch_gpu = cl.Buffer(ctx, mf.READ_ONLY, 36)
    events.append(cl.enqueue_copy(queue, pump_branch_gpu, pump_branch,
                                  device_offset=0, is_blocking=False))

    return events, pump_branch_gpu
    CLArg('pump_branch', 'gcfloat_p'),
    CLArg('omegas', 'gcfloat_p'),
    CLArg('seq_len', 'unsigned'),
    CLArg('gamma_total', 'gcfloat_p'),
    CLArg('delta_xyz', 'gcuint_p'),
    CLArg('omega_xyz_offset', 'gcuint_p')

    res, evt = solver.run(t0, t1, h, y0, queue,
                          extra_args=(np.float32(h_x), np.int64(len_x)))
    print('queued')
    evt.wait()
    print('finished')
    res_np = [a.get() for a in res]


def main():
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    dim_x = 100
    dim_y = 100
    dim_z = 30
    gamma_x = np.empty([dim_x, dim_x], np.float32)
    gamma_y = np.empty([dim_y, dim_y], np.float32)
    gamma_z = np.empty([dim_z, dim_z], np.float32)

    for i in range(dim_x):
        for j in range(dim_x):
            gamma_x[i, j] = sideband_scatter_strength(i, j, 0.2, 0)

    for i in range(dim_y):
        for j in range(dim_y):
            gamma_y[i, j] = sideband_scatter_strength(i, j, 0.2, np.pi / 2)

    for i in range(dim_z):
        for j in range(dim_z):
            gamma_z[i, j] = sideband_scatter_strength(i, j, 0.6, np.pi / 2)

    pump_branch = np.array([[0.500, 0.333, 0.167],
                            [0.006, 0.171, 0.002],
                            [0.500, 0.333, 0.167]], np.float32)
    omegas_x = None
    omegas_y = None
    omegas_z = None
    h_t = None
    gamma_totals = None
    delta_xyz = None
    omega_xyz = None

    events, pump_branch_gpu = evolve_sideband(ctx, queue, gamma_x, gamma_y,
                                              gamma_z, pump_branch, omegas_x,
                                              omegas_y, omegas_z, h_t,
                                              gamma_totals, delta_xyz,
                                              omega_xyz)

    cl.wait_for_events(events)

    res_np = np.empty(9, np.float32)
    cl.enqueue_copy(queue, res_np, pump_branch_gpu)
    print(res_np)

if __name__ == '__main__':
    main()
