library(openxlsx)

args <- commandArgs(trailingOnly = T)
file <- args[1]
outfile <- args[2]

data <- read.delim(file, comment.char = "#", sep = "\t", na.strings = "NA")
write.xlsx(data, outfile)
