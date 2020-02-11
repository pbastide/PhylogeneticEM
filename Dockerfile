FROM rocker/verse:latest

RUN R -e 'install.packages(c("lars", "iterators", "shape", "mvtnorm", "elasticnet", "randomForest", "pls", "gtools", "DEoptimR", "ape", "capushe", "foreach", "gglasso", "glmnet", "LINselect", "robustbase", "RcppArmadillo"))'
RUN R -e 'install.packages(c("devtools"))'

COPY . /Phylogenetic-EM

CMD R -e 'Rcpp::compileAttributes("Phylogenetic-EM", verbose = TRUE); devtools::check("Phylogenetic-EM", force_suggests = FALSE, run_dont_test = TRUE, document = FALSE, build_args = c("--no-build-vignettes", "--no-manual"), args = c("--no-vignettes", "--no-manual", "--as-cran"))'