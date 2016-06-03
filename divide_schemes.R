divide_calls <- function(fname) {

    calls <- read.csv(fname, stringsAsFactors = FALSE,
                      header = TRUE, row.names = 1, check.names = FALSE)


    find_missing_genes <- function(name) {

        perc_missing <- sum(which(calls[,name] == '0')) / length(calls[,name])
        if (perc_missing > 0.01) {
            TRUE
        }
        else {
            FALSE
        }

    }

    missing <- sapply(colnames(calls), find_missing_genes)

    accessory <- calls[,  missing]
    core      <- calls[, !missing]

    write.csv(accessory, file = "accessory_calls.csv", quote = FALSE)
    write.csv(core, file = "core_calls.csv", quote = FALSE)

    accessory_genes <- colnames(calls)[missing]

    accessory_genes
}

divide_markers <- function(accessory_genes, markers_path) {

    markers <- read.table(markers_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    accessory <- subset(markers, Marker.Name %in% accessory_genes)
    core <- subset(markers, !Marker.Name %in% accessory_genes)

    a_path <- paste(dirname(markers_path), "accessory.markers", sep = "/")
    c_path <- paste(dirname(markers_path), "core.markers", sep = "/")

    write.table(accessory, file = a_path, quote = FALSE,
                sep = "\t", row.names = FALSE, na = '')
    write.table(core, file = c_path, quote = FALSE,
                sep = "\t", row.names = FALSE, na = '')
}

args <- commandArgs(trailingOnly = TRUE)

accessory_genes <- divide_calls(fname = args[1])
divide_markers(accessory_genes, args[2])
