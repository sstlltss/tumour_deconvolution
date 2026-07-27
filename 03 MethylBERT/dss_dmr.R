library(bsseq)
library(DSS)
library(HDF5Array)

ctypes <- list("DLBCL","GC","PC")
file_path <- "path/to/cov/files"
labels <- read.csv(
  paste0(file_path,"train_labels.txt"),
  header = FALSE,
  col.names = c("sample", "type"),
  stringsAsFactors = FALSE
)

normal_files <- labels$sample[labels$type == "normal"]
normal_names <- normal_files

for (n in 1:22){
	chr = paste0("chr",n)
	print(paste0("===chr: ",chr,"==="))
	# load all normal samples on single chromosome
	normal_list <- c()
	print(paste0("Loading normal file:",normal_files[1],"_",chr,'.csv.gz'))
	normal_list <- read.bismark(
			files = paste0(file_path,normal_files[1],"/",normal_files[1],"_",chr,'.csv.gz'),
			rmZeroCov = FALSE,
			verbose = TRUE,
			BACKEND = "HDF5Array"
			)
	for (f in normal_files[2:length(normal_files)]){
		print(paste0("Loading normal file:",f,"_",chr,'.csv.gz'))
		next_normal <- read.bismark(
			files = paste0(file_path,f,"/",f,"_",chr,'.csv.gz'),
			rmZeroCov = FALSE,
			verbose = TRUE,
			BACKEND = "HDF5Array"
			)
		normal_list <- combine(normal_list, next_normal)
	}
	print("Set sample name: normal cells.")
	sampleNames(normal_list) <- normal_names
	print("Normal samples:")
	print(normal_list)

	for (ctype in ctypes){
		print(paste0("===ctype: ",ctype,"==="))
		cancer_files <- labels$sample[labels$type == ctype]
		print(paste0("Loading cancer file:",cancer_files[1],"_",chr,'.csv.gz'))
		cancer_list <- read.bismark(
				files = paste0(file_path,cancer_files[1],"/",cancer_files[1],"_",chr,'.csv.gz'),
				rmZeroCov = FALSE,
				verbose = TRUE,
				BACKEND = "HDF5Array"
				)
		for (f in cancer_files[2:length(cancer_files)]){
			if (file.exists(paste0(file_path,f,"/",f,"_",chr,'.csv.gz'))){
				print(paste0("Loading cancer file:",f,"_",chr,'.csv.gz'))
				next_cancer <- read.bismark(
					files = paste0(file_path,f,"/",f,"_",chr,'.csv.gz'),
					rmZeroCov = FALSE,
					verbose = TRUE,
					BACKEND = "HDF5Array"
					)
				cancer_list <- combine(cancer_list, next_cancer)
			}
		}

		cancer_names <- cancer_files

		print("Set sample name: cancer cells.")
		sampleNames(cancer_list) <- cancer_names

		print("Creating BSobj...")
		print("Cancer samples:")
		print(cancer_list)
		BSobj <- combine(cancer_list, normal_list)
		print("Adding group to BSobj...")
		pData(BSobj)$Group <- c(rep(ctype, length(cancer_files)),
					rep("normal", length(normal_files)))
		print("BSobj done. Start DMLtest...")

		dmlTest <- DMLtest(BSobj,
				group1 = cancer_names,
				group2 = normal_names,
				smoothing = TRUE)

		dmls <- callDML(dmlTest, p.threshold = 0.001)

		dmrs <- callDMR(dmlTest, delta = 0.1, p.threshold = 0.05,
						minCG = 3, dis.merge = 500)
		write.csv(dmrs, file = paste0("path/to/save/dmr/results", ctype,"_",chr,".csv"), row.names = FALSE)
	}
}