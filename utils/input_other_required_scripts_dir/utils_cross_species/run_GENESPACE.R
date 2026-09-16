library(GENESPACE)


args <- commandArgs(TRUE)
workding_dir <- as.character(args[1])
#ipt_path2mcscanx <- as.character(args[2])

#gpar <- init_genespace(
#  wd = workding_dir,
#  path2mcscanx = ipt_path2mcscanx)

##it does not work for not specifing the path2mcscanx
#gpar <- init_genespace(
#  wd = workding_dir)



gpar <- init_genespace(
  wd = workding_dir,
  path2mcscanx = file.path(Sys.getenv("CONDA_PREFIX"), "bin"))

gpar$shellCalls$orthofinder <- 'orthofinder'

#gpar$shellCalls$mcscanx_h <- 'MCScanX_h'

gpar <- run_genespace(gsParam = gpar)



#gpar <- init_genespace(
#  wd = "/public2/home/yanhaidong/working_dir/jinyarong/ACRs/Pm10",
#  path2mcscanx = "/public2/home/yanhaidong/software/MCScanX-master/")
#gpar <- run_genespace(gsParam = gpar)

