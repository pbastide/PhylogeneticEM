# PhylogeneticEM 1.1.0
* Deprecation Fix:
    * Remove dependency from gglasso (replaced with grplasso). 
    This should not change the results of PhylogeneticEM.
* Bug Fixes:
    * Issue an error message when the tree has some zero-length branches.
* Minor Changes:
    * Added options to plot.equivalent_shifts.

# PhylogeneticEM 1.0.1
* Minor Changes:
    * Impose a maximum value for alpha in find_grid_alpha to respect machine max double.
* New Features
    * added argument `label_font` to `plot` function to control the label font.
    * added argument `axis_las` to `plot` function to control the axis las.
* Bug fixes:
    * Plotting missing values correctly in plot.PhyloEM.
    * Bug fixes in plotting PhyloEM object when p = 1.
    * When p=1 and nbr_alpha is not missing, do not switch to estimated mode for alpha.
* Technical:
    * registration of c++ code to comply with R 3.4 new standards.

# PhylogeneticEM 1.0.0
Initial Release.
