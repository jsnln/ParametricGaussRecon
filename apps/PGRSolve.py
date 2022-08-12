# MIT License

# Copyright (c) 2022 Siyou Lin, Dong Xiao, Zuoqiang Shi, Bin Wang

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import numpy as np
from time import time
import argparse

from utils_pgr import *

CHUNK_SIZE = 512
FLT_TYPE = np.float32
R_SQ_STOP_EPS = 1e-20
TARGET_ISO_VALUE = -0.5     # this is compatible with the MC orientation

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--base', type=str, required=True, help='NPY file contatning the surface samples')
parser.add_argument('-s', '--sample', type=str, required=True, help='NPY file contatning the surface samples')
parser.add_argument('-q', '--query', type=str, required=True, help='NPY file contatning the octree grid corners for queries')
parser.add_argument('-o', '--output', type=str, required=True, help='output prefix')
parser.add_argument('-wk', '--width_k', type=int, required=True, help='k in knn for width estimation')
parser.add_argument('-wmin', '--width_min', type=float, required=True, help='minimum width, overrides --width_max')
parser.add_argument('-wmax', '--width_max', type=float, required=True, help='maximum width')
parser.add_argument('-a', '--alpha', type=float, required=True, help='alpha for regularization')
parser.add_argument('--max_iters', type=int, help='maximum iterations for CG')
parser.add_argument('--save_r', action='store_true', help='save the residual list')
parser.add_argument('--cpu', action='store_true', help='use cpu')
args = parser.parse_args()

if args.cpu:
    cp = None
else:
    import cupy as cp

if __name__ == '__main__':
    
    # env setup
    out_prefix = args.output

    y_base_np = load_sample_from_npy(args.base, return_cupy=False, dtype=FLT_TYPE)      # [N_x, 3]
    if args.sample == args.base:
        x_sample_np = y_base_np
    else:
        x_sample_np = load_sample_from_npy(args.sample, return_cupy=False, dtype=FLT_TYPE)  # [N_y, 3]
        
    if args.width_min > args.width_max:
        x_width_np = np.ones(x_sample_np.shape[0], dtype=FLT_TYPE) * args.width_min
        TIME_START_X_WIDTH = 0
        TIME_END_X_WIDTH = 0
    else:    
        TIME_START_X_WIDTH = time()
        x_width_np, base_kdtree = get_width(x_sample_np,
                                         k=args.width_k,
                                         dtype=FLT_TYPE,
                                         width_min=args.width_min,
                                         width_max=args.width_max,
                                         base_set=y_base_np,
                                         return_kdtree=True)
        TIME_END_X_WIDTH = time()
    
    
    x_sample = x_sample_np if args.cpu else cp.array(x_sample_np)
    x_width = x_width_np if args.cpu else cp.array(x_width_np)
    y_base = y_base_np if args.cpu else cp.array(y_base_np)
    
    print(f'[In apps.PGRSolve] x_width range: [{x_width.min().item():.4f}, {x_width.max().item():.4f}], mean: {x_width.mean().item():.4f}')
    print('\033[94m' + f'[Timer] x_width computed in {TIME_END_X_WIDTH-TIME_START_X_WIDTH}' + '\033[0m')
    
    print(f'[In apps.PGRSolve] Starting to solve the system...')
    solved = solve(x_sample,
                   y_base,
                   x_width,
                   chunk_size=CHUNK_SIZE,
                   dtype=FLT_TYPE,
                   iso_value=TARGET_ISO_VALUE,
                   r_sq_stop_eps=R_SQ_STOP_EPS,
                   alpha=args.alpha,
                   max_iters=args.max_iters,
                   save_r=args.save_r)
    if args.save_r:
        lse, r_list = solved
        out_r_list_txt = out_prefix + 'residuals.csv'
        np.savetxt(out_r_list_txt, r_list, fmt="%.16e", delimiter='\n')
    else:
        lse = solved
        
    if args.cpu:
        lse_np = lse
    else:
        lse_np = lse.get()
        cp._default_memory_pool.free_all_blocks()
    
    # saving solution as npy and xyz
    out_lse_array_npy = np.concatenate([y_base_np, -lse_np.reshape(3,-1).T], axis=1)
    out_solve_npy = out_prefix + 'lse'
    np.save(out_solve_npy, out_lse_array_npy)

    # saving solution as xyz
    out_solve_xyz = out_prefix + 'lse.xyz'
    np.savetxt(out_solve_xyz, out_lse_array_npy, fmt="%.8f", delimiter=' ')
    
    # eval on grid
    TIME_START_EVAL = time()
    q_query = load_sample_from_npy(args.query, return_cupy=False, dtype=FLT_TYPE)

    if args.width_min >= args.width_max:
        q_width = np.ones(q_query.shape[0], dtype=FLT_TYPE) * args.width_min
        TIME_START_Q_WIDTH = 0
        TIME_END_Q_WIDTH = 0
    else:
        TIME_START_Q_WIDTH = time()
        q_width = get_width(q_query,
                            k=args.width_k,
                            dtype=FLT_TYPE,
                            width_min=args.width_min,
                            width_max=args.width_max,
                            base_kdtree=base_kdtree,
                            return_kdtree=False)
        TIME_END_Q_WIDTH = time()
    
    print(f'[In apps.PGRSolve] q_width range: [{q_width.min().item():.4f}, {q_width.max().item():.4f}]')
    print('\033[94m' + f'[Timer] q_width computed in {TIME_END_Q_WIDTH-TIME_START_Q_WIDTH}' + '\033[0m')
    print('\033[94m' + f'[Timer] both width computed in {TIME_END_X_WIDTH-TIME_START_X_WIDTH+TIME_END_Q_WIDTH-TIME_START_Q_WIDTH}' + '\033[0m')
    
    sample_vals = get_query_vals(x_sample, x_width, y_base, lse, CHUNK_SIZE)
    iso_val = float(np.median(sample_vals))
    print(f'[In apps.PGRSolve] sample vals range: [{sample_vals.min().item():.4f}, {sample_vals.max().item():.4f}], mean: {sample_vals.mean().item():.4f}, median: {np.median(sample_vals).item():.4f}')
    out_isoval_txt = out_prefix + 'isoval.txt'
    with open(out_isoval_txt, 'w') as isoval_file:
        isoval_file.write(f'{iso_val:.8f}')
        
    CHUNK_SIZE = 1024
    if args.cpu:
        CHUNK_SIZE = 16384
    query_vals = get_query_vals(q_query, q_width, y_base, lse, CHUNK_SIZE)
    
    out_grid_width_npy = out_prefix + 'grid_width'
    print(f'[In apps.PGRSolve] Saving grid widths to {out_grid_width_npy}')
    np.save(out_grid_width_npy, q_width)

    out_eval_grid_npy = out_prefix + 'eval_grid'    
    print(f'[In apps.PGRSolve] Saving grid eval values to {out_eval_grid_npy}')
    np.save(out_eval_grid_npy, query_vals)
    
    TIME_END_EVAL = time()
    print('\033[94m' + f'[Timer] Eval on grid finished in {(TIME_END_EVAL-TIME_START_EVAL)-(TIME_END_Q_WIDTH-TIME_START_Q_WIDTH)}' + '\033[0m')
