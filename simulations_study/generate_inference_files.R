## This is a helper script to distribute inference computations using
## data threads

## Number of parallel data threads
N <- 10
## Inference template file
file <- "dummy.R"

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
    range <- nchunks[[i]]
    n.range <- paste0("n.range <- c(",
                      paste0(range, collapse = ", "), ")")
    data.i <- c(data, n.range)
    ## reorder lines of r script
    data.i <- data.i[c(1:line2change, length(data.i),
                       seq_along(data.i)[-c(1:line2change, length(data.i))])]
    write(data.i, filename)
}

for (i in seq_along(nchunks)) {
    replace_n(i)
}
