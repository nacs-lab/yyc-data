#!/usr/bin/env python

import os.path as _path
from pyscical.ocl.ode import ElwiseOdeSolver
from pyscical.ocl.utils import CLArg
import numpy as np
import pyopencl as cl

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

    print(events)
    return

    res, evt = solver.run(t0, t1, h, y0, queue,
                          extra_args=(np.float32(h_x), np.int64(len_x)))
    print('queued')
    evt.wait()
    print('finished')
    res_np = [a.get() for a in res]


def main():
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    gamma_x = np.ones([100, 100]).astype(np.float32)
    gamma_y = np.zeros([100, 100]).astype(np.float32)
    gamma_z = np.ones([30, 30]).astype(np.float32)

    pump_branch = None
    omegas_x = None
    omegas_y = None
    omegas_z = None
    h_t = None
    gamma_totals = None
    delta_xyz = None
    omega_xyz = None

    evolve_sideband(ctx, queue, gamma_x, gamma_y, gamma_z, pump_branch,
                    omegas_x, omegas_y, omegas_z, h_t, gamma_totals,
                    delta_xyz, omega_xyz)

if __name__ == '__main__':
    main()
