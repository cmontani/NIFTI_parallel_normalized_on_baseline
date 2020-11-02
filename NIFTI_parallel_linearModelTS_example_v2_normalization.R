#Set-up librarier
library(R.utils)
library(tidyr)
library(dplyr)
library(Hmisc)
library(lme4)
library(nifti.io)
library(doParallel)

library(purrr)
library(stringr)

#  install.packages("pracma", repos="http://R-Forge.R-project.org")
#   library(pracma)

source("/media/DATA/cmontani/IIT_20May_Oxtcre-Gq_BOLD_PV6/Script_and_template/R/NIFTI_parallel_v2/gunzip.R")
source("/media/DATA/cmontani/IIT_20May_Oxtcre-Gq_BOLD_PV6/Script_and_template/R/NIFTI_parallel_v2/table.to.nii.R")

numcores <- 10 #Change this depending on your workstation
len_timeseries <- 4570

#Create (or load) data frame with information about each subject
base_path <- "/media/DATA/cmontani/IIT_20May_Oxtcre-Gq_BOLD_PV6/37_modeling_ROIs_or_voxelwise_and_sliding_window/Linear_model/on_normalized_BOLD_entire_ts/"
setwd(base_path)
subjects <- readLines(con="list_subjects.txt") 
df <- data.frame(subject=subjects, stringsAsFactors = FALSE) %>%
         mutate(BOLD_file =subject,
         genotype = stringr::str_split(subject, "_") %>% map_chr(., 3))

#Before starting all files need to be unzipped for fast reading
scratch_dir <- paste0(base_path,"temporary/")
if (!dir.exists(scratch_dir)) {
  dir.create(path=scratch_dir, recursive=TRUE)
}
data_nii <- df$BOLD_file

 for (i in 1:length(data_nii)) {
  new_fname <- paste(scratch_dir, sub("\\.gz$","", basename(data_nii[i])), sep="/")
  if (!file.exists(new_fname))  {
    gunzip(filename=data_nii[i], destname=new_fname, overwrite = FALSE, remove=FALSE) # remove = FALSE is important
   }
 data_nii[i] <- new_fname
 }


#Create a dataframe that has information about each subject needed for the model
df_model <- data.frame(tr=1:len_timeseries) %>% 
  merge(df, by=NULL) %>% select(subject, tr,genotype)

#Mask file, unzip if needed
data_mask <- paste0(base_path,"chd8_functional_template_mask_wo_cerebellum.nii.gz")
new_fname <- paste(scratch_dir, sub("\\.gz$", "", basename(data_mask)), sep="/")
gunzip(filename=data_mask, destname=new_fname, overwrite = FALSE, remove=FALSE) # remove = FALSE is important
data_mask <- new_fname

#Directory to put the new data


save_dir <- "/media/DATA/cmontani/IIT_20May_Oxtcre-Gq_BOLD_PV6/37_modeling_ROIs_or_voxelwise_and_sliding_window/Linear_model/on_normalized_BOLD_entire_ts/output_groups/"
if (!dir.exists(save_dir)) {
  dir.create(path=save_dir, recursive=TRUE)
}


### Load parameters (size and orientation) of desired output  
eg_file <- data_nii[1]
pixdim <- unlist(nii.hdr(eg_file, "pixdim"))
orient <- nii.orient(eg_file)
dims <- nii.dims(eg_file)

mask_voxels <- load.mask(data_mask)

#This is just for debugging
# concat_timeseries <- c()
# for(i in seq(length(data_nii))){
#   thistimeseries <- read.nii.voxel(nii.file = data_nii[i], coords = c(mask_voxels[[1]][200,],Inf))
#   concat_timeseries <- c(concat_timeseries, thistimeseries)
# }
# timeseries <- concat_timeseries


#Function for calculating lmer

myfunct <- function(timeseries, df_model){
    df_model$r <- timeseries
    mymice <- unique(df_model$subject)
    for (i in seq(length(mymice))){
    sub_df <- df_model %>% filter(subject == mymice[i])
    baseline <- sub_df$r[1:900]
    mean_baseline <- mean(baseline)
    df_model$r_norm[df_model$subject == mymice[i]] <- sub_df$r - mean_baseline   #normalize on the baseline 
    }
   
    mylm <- lm(r_norm~tr*genotype, df_model)     #stat between groups
    myOut <- summary(mylm)$coefficients
    return(myOut)
     }
  
 

  


#Start parallel clustering
cl <- makeCluster(numcores)
registerDoParallel(cl)

#Run loop
result <- foreach(j = seq(dim(mask_voxels[[1]])[1]),
                  .packages=c('nifti.io','dplyr')) %dopar% {
  concat_timeseries <- c()
  for(i in seq(length(data_nii))){
    thistimeseries <- read.nii.voxel(nii.file = data_nii[i], coords = c(mask_voxels[[1]][j,],Inf))
    concat_timeseries <- c(concat_timeseries, thistimeseries)
  }
#   thismodel <- lm(concat_timeseries ~ tr*genotype, df_model)
  #   summary(thismodel)$coefficients
  myfunct(timeseries = concat_timeseries, df_model = df_model)
}

stopCluster(cl) 

num_coeffs <- dim(result[[1]])[1] 


result2 <- do.call("rbind", result) %>% as.data.frame()


result3 <- tibble::rownames_to_column(result2, "VALUE") %>%
  mutate(num = gsub('.*\\.', '', VALUE))


for (i in seq(num_coeffs)) {
  result3$num[i] <- 0
}

result3 <- result3 %>% mutate(tmp = str_remove(VALUE,num))
result3$coeff <- gsub(".","",result3$tmp, fixed=TRUE)


result4 <- result3 %>%
  pivot_wider(id_cols = num, names_from = coeff, values_from = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")) %>%
  select(-num)





#Write results from liner model to multiple files
table.to.nii(in.table = result4, coords = mask_voxels[[1]], img.dims = dims, save.dir = save_dir, 
             prefix = "test", model = "concat_timeseries ~ tr*genotype", pixdim = pixdim, orient = orient)

#Everything is done running
#Time to clean up

#Remove extra data frames
rm(result, result2, result3, result4)

#Zip all the output files
thisfilelist <- list.files(save_dir)
for (i in 2:length(thisfilelist)) {
  fname <- paste(save_dir, thisfilelist[i], sep="/")
  gzip(fname)
}

#Delete scratch directory
unlink(scratch_dir, recursive = TRUE)

