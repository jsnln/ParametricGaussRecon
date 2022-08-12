from time import time
import argparse
import os
from os.path import join, basename
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('input', help='xyz format point cloud input, this should be the one to construct the final query octree')
parser.add_argument('-wk', '--width_k', type=int, default=7, help='k in knn for width estimation')
parser.add_argument('-wmax', '--width_max', type=float, default=0.015, help='minimum width, overrides --width_max')
parser.add_argument('-wmin', '--width_min', type=float, default=0.0015, help='maximum width')
parser.add_argument('-a', '--alpha', type=float, default=1.05, help='alpha for regularization')
parser.add_argument('-m', '--max_iters', type=int, help='maximum iterations for CG')
parser.add_argument('-d', '--max_depth', type=int, default=10, help='max depth of octree')
parser.add_argument('-md', '--min_depth', type=int, default=1, help='min depth of octree')
parser.add_argument('--cpu', action='store_true', help='run with cpu')
parser.add_argument('--save_r', action='store_true', help='save the residual list')
args = parser.parse_args()

EXPORT_QUERY_EXE = './apps/PGRExportQuery'
LOAD_QUERY_EXE = './apps/PGRLoadQuery'
SOLVE_APP = 'apps.PGRSolve'
PARAM_MIDFIX = f'_k_{args.width_k}_min_{args.width_min}_max_{args.width_max}_alpha_{args.alpha}_depth_min_{args.min_depth}_depth_max_{args.min_depth}_'

setup_str = f'---------Settings---------\n'\
          + f'min depth:   {args.min_depth}\n'\
          + f'max depth:   {args.max_depth}\n'\
          + f'max iters:   {args.max_iters}\n'\
          + f'width_k:     {args.width_k}\n'\
          + f'width_min:   {args.width_min}\n'\
          + f'width_max:   {args.width_max}\n'\
          + f'alpha:       {args.alpha}\n'\
          + f'---------Settings---------'

print(setup_str)

in_filename = args.input
data_index = in_filename.split('/')[-1][:-4]

results_folder = 'results/'
sample_file_folder = results_folder + data_index + '/samples/'
solve_file_folder = results_folder + data_index + '/solve/'
recon_file_folder = results_folder + data_index + '/recon/'

sample_file_prefix = sample_file_folder + data_index
solve_file_prefix = solve_file_folder + data_index
recon_file_prefix = recon_file_folder + data_index

solve_base_file_prefix = sample_file_prefix
solve_sample_file_prefix = sample_file_prefix

if not os.path.exists(sample_file_folder):
    os.makedirs(f'{sample_file_folder}')
    shutil.copyfile(in_filename, join(sample_file_folder, basename(in_filename)))
    # os.system(f'cp {in_filename} {sample_file_folder}')
if not os.path.exists(solve_file_folder):
    os.makedirs(f'{solve_file_folder}')
if not os.path.exists(recon_file_folder):
    os.makedirs(f'{recon_file_folder}')


# build octree
TIME_START_OCTREE = time()
build_octree_cmd = f'{EXPORT_QUERY_EXE} -i {args.input} -o {sample_file_prefix} -d {args.max_depth} -m {args.min_depth} '
print(f'\n[EXECUTING] {build_octree_cmd}\n')
os.system(build_octree_cmd)
TIME_END_OCTREE = time()

# solve equation
TIME_START_SOLVE = time()
solve_cmd = f"python -m {SOLVE_APP} "\
          + f"--base {solve_base_file_prefix}_normalized.npy "\
          + f"--sample {solve_sample_file_prefix}_normalized.npy "\
          + f"--query {sample_file_prefix}_for_query.npy "\
          + f"--width_k {args.width_k} "\
          + f"--width_max {args.width_max} "\
          + f"--width_min {args.width_min} "\
          + f"--alpha {args.alpha} "\
          + f"-o {solve_file_prefix}{PARAM_MIDFIX} "
if args.max_iters is not None:
    solve_cmd += f"--max_iters {args.max_iters} "
if args.save_r:
    solve_cmd += f"--save_r "
if args.cpu:
    solve_cmd += f"--cpu "

print(f'\n[EXECUTING] {solve_cmd}\n')
os.system(solve_cmd)
TIME_END_SOLVE = time()

# reconstruction
TIME_START_RECON = time()
with open(f'{solve_file_prefix}{PARAM_MIDFIX}isoval.txt', 'r') as isoval_file:
    isoval = isoval_file.read()
    isoval = eval(isoval)

recon_cmd = f"{LOAD_QUERY_EXE} -i {in_filename} -d {args.max_depth} -m {args.min_depth} "\
          + f"--grid_val {solve_file_prefix}{PARAM_MIDFIX}eval_grid.npy "\
          + f"--grid_width {solve_file_prefix}{PARAM_MIDFIX}grid_width.npy "\
          + f"--isov {isoval} "\
          + f"-o {recon_file_prefix}{PARAM_MIDFIX}recon.ply"
print(f'\n[EXECUTING] {recon_cmd}\n')
os.system(recon_cmd)
TIME_END_RECON = time()

print('\033[94m' + f'[Timer] Note: Some preprocessing (width computation) is actually in the Main part.' + '\033[0m')
print('\033[94m' + f'[Timer] Pre:    {TIME_END_OCTREE-TIME_START_OCTREE}' + '\033[0m')
print('\033[94m' + f'[Timer] Main:   {TIME_END_SOLVE-TIME_START_SOLVE}' + '\033[0m')
print('\033[94m' + f'[Timer] Post:   {TIME_END_RECON-TIME_START_RECON}' + '\033[0m')
print('\033[94m' + f'[Timer] Total:  {TIME_END_OCTREE-TIME_START_OCTREE+TIME_END_SOLVE-TIME_START_SOLVE+TIME_END_RECON-TIME_START_RECON}' + '\033[0m')


