% where is your eigen? (Change appropriately)
cg_eigen_dir = '/usr/local/include/eigen3/'; 


% mex solvers
bp = 'solvers/';

fname = [bp 'solver_triang_opt.cpp'];
mex(['-I"' cg_eigen_dir '"'],'-O',fname,'-outdir',bp)

fname = [bp 'solver_triang_null.cpp'];
mex(['-I"' cg_eigen_dir '"'],'-O',fname,'-outdir',bp)

