#!/home/hal9000/library/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
plotData <- function(filename) {
#	 print(filename)
	 tmp <- strsplit(filename, "/")[[1]]
#	 print(tmp)
	 fname <- tail(tmp, 1)
#	 print(fname)
	 data<-read.table(filename,header=TRUE)

	 # Dump pdf
	 outfile<-paste(filename, ".pdf", sep="")
	 pdf(outfile)
	 plot(data)
	 dev.off()

	 # Save data into RData
#	 load(".RData", .GlobalEnv)
#	 assign(fname, data, envir = .GlobalEnv)
#	 save.image(file = ".RData");
#	 save(fname, file = ".RData")

}


print(args[1])
plotData(args[1])
