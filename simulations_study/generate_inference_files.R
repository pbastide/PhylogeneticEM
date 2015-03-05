## This is a helper script to distribute inference computations using
## data threads

## Number of parallel data threads
N <- 10
## Inference template file
folder <- "simulations_study"
file <- file.path(folder, "dummy.R")
## Folder to keep the generated files
datestamp_day <- format(Sys.time(), "%Y-%m-%d")
new_folder <- file.path(folder, datestamp_day)
dir.create(new_folder)

## Recover total number of simulations, assumes that the
## first line of the form n <- xxx in file gives the number of simulations
data <- readLines(file)
line2change <- grep("n <-", data)[1]
number.simulations <- data[line2change]
eval(parse(text = number.simulations))

## split n in N chunks of equal size
nchunks <- split(n, ceiling(seq_along(n)/(length(n)/N)))

## replace the range of values of n over which inference is done
## and write resulting R script in a new file
replace_n <- function(i, filename = file) {
    filename <- sub(".R", paste0("_", i, ".R"), file)
    filename <- sub(folder, new_folder, filename)
    range <- nchunks[[i]]
    n.range <- paste0("n.range <- c(",
                      paste0(range, collapse = ", "), ")")
    index <- paste0("inference.index <- ", i)
    data.i <- c(data, n.range, index)
    ## reorder lines of r script
    data.i <- data.i[c(1:line2change, length(data.i) - 0:1,
                       seq_along(data.i)[-c(1:line2change, length(data.i) - 0:1)])]
    write(data.i, filename)
}

for (i in seq_along(nchunks)) {
    replace_n(i)
}
