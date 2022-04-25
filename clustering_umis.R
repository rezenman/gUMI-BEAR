start = Sys.time()
getwd()
library(ShortRead)
library(stringr)
library(dplyr)
print("finished loading packages")

##reading file name from what provided by commandline
args = commandArgs(trailingOnly=TRUE)
print(paste0("processing file ", args))

##getting experiment name
split_wd = unlist(strsplit(getwd(), "/"))
exp_name = split_wd[length(split_wd)]


##loading seeds umi file
print("loading seeds umi file")
umi_seeds_loc = list.files("../umis/", pattern = paste0(exp_name, "_umis_Seeds.fastq"), full.names = T)
umi_seeds = readFastq(umi_seeds_loc) # change to the location of your umi_seeds file
df_umis_seeds = data.frame(id = umi_seeds@id, UMI = umi_seeds@sread)
new_ID_col_umis_seeds = sapply(strsplit(df_umis_seeds$id, " "), `[`, 1)
df_umis_seeds$id = new_ID_col_umis_seeds

##looking which lib name matches sample name
sample_name = strsplit(args, "-|\\.1\\.")[[1]][2] 
umi_file = list.files("../umis/", pattern = sample_name, full.names = T)
print(paste0("using umi file ", umi_file))

##loading umi data frame
umis = readFastq(umi_file)
df_umis = data.frame(id = umis@id, UMI = umis@sread)
new_ID_col_umis = sapply(strsplit(df_umis$id, " "), `[`, 1)
df_umis$id = new_ID_col_umis
print("finished loading umis df")
head(df_umis)
dim(df_umis)

##binding the umis seeds table with the sample umis table
df_umis = rbind(df_umis, df_umis_seeds)
df_umis = distinct(df_umis)

##getting seeds data frame
seeds_file = list.files(pattern = "seeds.clst$", full.names = T)
print(paste0("loading seeds file: ", seeds_file)) 
seeds = readFastq(seeds_file)
df_seeds = data.frame(id = seeds@id, read_1 = seeds@sread)
new_id_col_seeds = sapply(strsplit(df_seeds$id, " "), `[`, 1)
df_seeds$id = new_id_col_seeds
print("finished loading seeds df")
print("processing sample")
        current_df = args
        ##loading clusterd data frame
        df = read.delim(current_df)               
        print(head(df))
        print(dim(df))
        print(paste0("1. finished loading clusters df:", current_df))
        
        ##merging umi table with clusters table
        merged = merge(df_umis, df, by = "id")
        merged = merged[order(merged$clstr),]
        merged = merged[merged$clstr_size != 1,]
        print(head(merged))
        print(dim(merged))
        print(paste0("2. finished merging umi's with:", current_df))
        
        ##loking for redundant umis
        dis_1 = merged[merged$clstr_rep == 1, ]
        print(head(dis_1))
        print(dim(dis_1))
        dis_2 = distinct(merged, clstr , UMI, .keep_all = T)
        print(head(dis_2))
        print(dim(dis_2))
        df_non_unique_umis = as.data.frame(table(dis_2$clstr))
        dis_1$real_barcodes = df_non_unique_umis$Freq
        print(head(dis_1))

        ##Adding seed sequences to each dis_1 data frame
        final_df = merge(dis_1, df_seeds, by = "id")        
        
        save(final_df, merged, dis_2, file = paste0("final_df", current_df, ".Rdata"))
        print(paste0("finished counting umis in:", current_df))
end = Sys.time()
print(paste0("time took: ", end-start))

