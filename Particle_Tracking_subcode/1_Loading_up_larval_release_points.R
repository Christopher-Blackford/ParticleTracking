########################################################################
########################################################################
########################################################################
########################################################################
### [1] Loading up larval release points

#Acquiring files
filenames <- list.files(path = "./cuke_present/ReleaseLocations", pattern="rl_.", full.names=TRUE,recursive=T)

# load all files into a list, read_csv is much faster than read.csv
rllist <- lapply(filenames, read_csv,
                 col_names = c("long0","lat0","Z0","delay","site0"),
                 col_types = cols("d","d","i","i","i")
)

# set the names of the items in the list, so that you know which file it came from
rllist <- setNames(rllist,filenames)

# rbind the list
rl <- rbindlist(rllist, idcol="filename")

rl$bin <- as.numeric(gsub(".*rl_|.txt.*", "",rl$filename))
head(rl)
rm(rllist, filenames)
