# PhylogeneticEM 1.4.0
* New option:
  * Add an option "show_axis_traits" to control the trait axis in a plot.
* New Feature
  * Function `simul_process` can now handle an "OU/BM" process (i.e. an OU with
  some eigen values equal to zero).
* Bug fix:
  * Documented arguments not in usage in documentation object `equivalent_shifts`.
  * Check with "_R_CHECK_FORCE_SUGGESTS_=FALSE" now runs.

# PhylogeneticEM 1.3.0
* Bug fix:
  * Bug fix in plot function for univariate result.
  * Tree does not need to be ultrametric for a BM
  * Change `tol` to `tol_EM` in `estimateEM` for actual binding
  * Bug fixes in plot function for displaying traits.
* New Option:
  * Add option "name_trait" to plot function for ancestral state reconstruction.
* API Change
  * For a univariate trait, the BM re-scaling is now applied by default.
* Deprecation Fix:
  * Update `ape`, and replace `is.binary.tree` by `is.binary`.
  * Remove dependency to Rphylopars (no longer on CRAN).

# PhylogeneticEM 1.2.1
* CRAN maintenance
  * Added references in Description
  * Reduce checking and vignette building time

# PhylogeneticEM 1.2.0
* Added Options:
    * `option_is.ultrametric` to control the `ape::is.ultrametric` check.
* Deprecation Fix:
    * Go back to updated version of gglasso (faster than grplasso). 
    This should not change the results of PhylogeneticEM.
* Added functions:
    * Functions `find_rotation` and `merge_rotations` to deal with datasets that are
    equal up to a rotation.
    * Functions `get_criterion` and `plot_criterion` to get and plot the criterion used
    for model selection easily.
* New feature:
    * Function PhyloEM can now allow for negative values of alpha, to fit an EB model.
* Bug fix:
    * In function plot, fixed bug that led to the wrong parameters being selected when 
    doing ancestral trait reconstruction.

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
