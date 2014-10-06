#!/usr/bin/env python

import os.path as _path
import numpy as np
import pyopencl as cl
from cffi import FFI
from six.moves import range as _range

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

def _get_gidx_minmax_xyz(dim_x, dim_y, dim_z, gamma_x, gamma_y, gamma_z):
    # Probably faster if done on GPU
    gidx_minmax_xyz = np.empty((dim_x + dim_y + dim_z) * 2, np.uint32)

    gidxmin_x = gidx_minmax_xyz[:dim_x]
    gidxmin_y = gidx_minmax_xyz[dim_x:dim_x + dim_y]
    gidxmin_z = gidx_minmax_xyz[dim_x + dim_y:dim_x + dim_y + dim_z]

    gidxmax_x = gidx_minmax_xyz[dim_x + dim_y + dim_z:
                                    dim_x * 2 + dim_y + dim_z]
    gidxmax_y = gidx_minmax_xyz[dim_x * 2 + dim_y + dim_z:
                                dim_x * 2 + dim_y * 2 + dim_z]
    gidxmax_z = gidx_minmax_xyz[dim_x * 2 + dim_y * 2 + dim_z:]
    _lib._fill_gidx_minmax(dim_x, cffi_ptr(gamma_x, _ffi)[0],
                           cffi_ptr(gidxmin_x, _ffi, writable=True)[0],
                           cffi_ptr(gidxmax_x, _ffi, writable=True)[0])
    _lib._fill_gidx_minmax(dim_y, cffi_ptr(gamma_y, _ffi)[0],
                           cffi_ptr(gidxmin_y, _ffi, writable=True)[0],
                           cffi_ptr(gidxmax_y, _ffi, writable=True)[0])
    _lib._fill_gidx_minmax(dim_z, cffi_ptr(gamma_z, _ffi)[0],
                           cffi_ptr(gidxmin_z, _ffi, writable=True)[0],
                           cffi_ptr(gidxmax_z, _ffi, writable=True)[0])
    return gidx_minmax_xyz

def evolve_sideband(ctx, queue, gamma_x, gamma_y, gamma_z, pump_branch,
                    omegas_x, omegas_y, omegas_z, h_t, gamma_total,
                    delta_xyz, omega_xyz, p_b, p_a=None, p_c=None):
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
    events = []

    gamma_xyz = cl.Buffer(ctx, mf.READ_ONLY,
                          (dim_x**2 + dim_y**2 + dim_z**2) * 4)
    events.append(cl.enqueue_copy(queue, gamma_xyz, gamma_x,
                                  device_offset=0, is_blocking=False))
    events.append(cl.enqueue_copy(queue, gamma_xyz, gamma_y,
                                  device_offset=dim_x**2 * 4,
                                  is_blocking=False))
    events.append(cl.enqueue_copy(queue, gamma_xyz, gamma_z,
                                  device_offset=(dim_x**2 + dim_y**2) * 4,
                                  is_blocking=False))

    gidx_minmax_xyz = cl.Buffer(ctx, mf.READ_ONLY, (dim_x + dim_y + dim_z) * 8)
    events.append(cl.enqueue_copy(queue, gidx_minmax_xyz,
                                  _get_gidx_minmax_xyz(dim_x, dim_y, dim_z,
                                                       gamma_x, gamma_y,
                                                       gamma_z),
                                  device_offset=0, is_blocking=False))

    if pump_branch.dtype != np.float32:
        raise TypeError("The type of pump_branch should be float32.")
    pump_branch_gpu = cl.Buffer(ctx, mf.READ_ONLY, 36)
    events.append(cl.enqueue_copy(queue, pump_branch_gpu, pump_branch,
                                  device_offset=0, is_blocking=False))

    num_omg_x, d = omegas_x.shape
    if dim_x != d:
        raise ValueError("The second dimension of omegas_x is not "
                         "the same with dim_x.")
    if omegas_x.dtype != np.float32:
        raise TypeError("The type of omegas_x should be float32.")

    num_omg_y, d = omegas_y.shape
    if dim_y != d:
        raise ValueError("The second dimension of omegas_y is not "
                         "the same with dim_y.")
    if omegas_y.dtype != np.float32:
        raise TypeError("The type of omegas_y should be float32.")

    num_omg_z, d = omegas_z.shape
    if dim_z != d:
        raise ValueError("The second dimension of omegas_z is not "
                         "the same with dim_z.")
    if omegas_z.dtype != np.float32:
        raise TypeError("The type of omegas_z should be float32.")

    omegas_gpu = cl.Buffer(ctx, mf.READ_ONLY,
                           (num_omg_x * dim_x + num_omg_y * dim_y +
                            num_omg_z * dim_z) * 4)
    events.append(cl.enqueue_copy(queue, omegas_gpu, omegas_x,
                                  device_offset=0, is_blocking=False))
    events.append(cl.enqueue_copy(queue, omegas_gpu, omegas_y,
                                  device_offset=num_omg_x * dim_x * 4,
                                  is_blocking=False))
    events.append(cl.enqueue_copy(queue, omegas_gpu, omegas_z,
                                  device_offset=(num_omg_x * dim_x +
                                                 num_omg_y * dim_y) * 4,
                                  is_blocking=False))

    h_t = np.float32(h_t)
    seq_len, d = gamma_total.shape
    t_len = (seq_len - 1) * h_t
    if d != 3:
        raise TypeError("Second dimension of gamma_total should be 3.")
    if gamma_total.dtype != np.float32:
        raise TypeError("The type of gamma_total should be float32.")
    gamma_total_gpu = cl.Buffer(ctx, mf.READ_ONLY, seq_len * 12)
    events.append(cl.enqueue_copy(queue, gamma_total_gpu, gamma_total,
                                  device_offset=0, is_blocking=False))

    d1, d2 = delta_xyz.shape
    if d1 != 3 or d2 != seq_len:
        raise TypeError("Dimensions of delta_xyz should be (3, seq_len).")
    if delta_xyz.dtype != np.uint32:
        raise TypeError("The type of delta_xyz should be uint32.")
    delta_xyz_gpu = cl.Buffer(ctx, mf.READ_ONLY, seq_len * 12)
    events.append(cl.enqueue_copy(queue, delta_xyz_gpu, delta_xyz,
                                  device_offset=0, is_blocking=False))

    d1, d2 = omega_xyz.shape
    if d1 != 3 or d2 != seq_len:
        raise TypeError("Dimensions of omega_xyz should be (3, seq_len).")
    if omega_xyz.dtype != np.uint32:
        raise TypeError("The type of omega_xyz should be uint32.")
    omega_xyz_offset = np.empty(seq_len * 3, np.uint32)
    for i in _range(seq_len):
        _omega_x = omega_xyz[0, i]
        if _omega_x >= num_omg_x:
            raise IndexError("omega_x index too larger")
        omega_xyz_offset[i] = _omega_x * dim_x
        _omega_y = omega_xyz[1, i]
        if _omega_y >= num_omg_y:
            raise IndexError("omega_y index too larger")
        omega_xyz_offset[seq_len + i] = _omega_y * dim_y + num_omg_x * dim_x
        _omega_z = omega_xyz[2, i]
        if _omega_z >= num_omg_z:
            raise IndexError("omega_z index too larger")
        omega_xyz_offset[seq_len * 2 + i] = (_omega_z * dim_z +
                                             num_omg_x * dim_x +
                                             num_omg_y * dim_y)
    omega_xyz_offset_gpu = cl.Buffer(ctx, mf.READ_ONLY, seq_len * 12)
    events.append(cl.enqueue_copy(queue, omega_xyz_offset_gpu, omega_xyz_offset,
                                  device_offset=0, is_blocking=False))


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
    seq_len = np.uint32(seq_len)

    return events, omega_xyz_offset_gpu
    y0

    res, evt = solver.run(0, t_len, h_t, y0, queue,
                          extra_args=(np.uint32(dim_x), np.uint32(dim_y),
                                      np.uint32(dim_z), gamma_xyz,
                                      gidx_minmax_xyz, pump_branch_gpu,
                                      omegas_gpu, h_t, seq_len, gamma_total_gpu,
                                      delta_xyz_gpu, omega_xyz_offset_gpu))
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
    num_omg_x = 10
    num_omg_y = 10
    num_omg_z = 10
    omegas_x = np.ones([num_omg_x, dim_x], np.float32)
    omegas_y = np.zeros([num_omg_y, dim_y], np.float32)
    omegas_z = np.ones([num_omg_z, dim_z], np.float32)
    h_t = 0.1
    seq_len = 1000
    gamma_total = np.ones([seq_len, 3], np.float32)
    delta_xyz = np.ones([3, seq_len], np.uint32)
    omega_xyz = np.ones([3, seq_len], np.uint32)
    p_b = None

    events, omega_xyz_offset_gpu = evolve_sideband(ctx, queue, gamma_x, gamma_y,
                                                   gamma_z, pump_branch, omegas_x,
                                                   omegas_y, omegas_z, h_t,
                                                   gamma_total, delta_xyz,
                                                   omega_xyz, p_b)

    cl.wait_for_events(events)

    res_np = np.empty(seq_len * 3, np.uint32)
    cl.enqueue_copy(queue, res_np, omega_xyz_offset_gpu)
    print(res_np)

if __name__ == '__main__':
    main()
