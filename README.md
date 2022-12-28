# cylinder-triangulation

Code for robust triangulation of 3D-cylinders from silhouette lines and known camera poses

To setup in Matlab add subfolders to path. Possibly you neeed to mex files (see script `mex_solvers.m`)

`triangulate_cylinders_wor_opt.m`  Full least squares 3D triangulation

`triangulate_circle_opt.m` Least squares circle triangulation

`triangulate_circle_minimal.m` Minimal circle triangulation

`triangulate_circle_ransac.m` RANSAC wrapper for robust circle triangulation

For examples run `test_triangulate_circle.m` and `triangulate_kungshuset.m` in the `test` folder.

If you use the code please cite
@article{gummeson2022robust,
  title={Robust and Accurate Cylinder Triangulation},
  author={Gummeson, Anna and Oskarsson, Magnus},
  journal={arXiv preprint arXiv:2212.02319},
  year={2022}
}

