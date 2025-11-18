# Setup
rm(list = ls(all.names = TRUE))
setwd("~/Dropbox/PAPERS/projects/MFDA")

source("helper_functions.R")
     
# Inputs
band  <- "metallica"

year_data  <- c(
     1983,
     1984,
     1986,
     1987,
     1988,
     1991,
     1996,
     1997,
     1998,
     2003,
     2008,
     2016,
     2023
)

album_data <- c(
     "kill_em_all",
     "ride_the_lightning",
     "master_of_puppets",
     "garage_days_rerevisited",
     "and_justice_for_all",
     "metallica",
     "load",
     "reload",
     "garage_inc",
     "st_anger",
     "death_magnetic",
     "hardwired_to_self_destruct",
     "72_seasons"
)

for (i in 1:length(album_data)){
     # Year and album 
     year  <- year_data[i]
     album <- album_data[i]
     
     # Build base path and input file
     base_stub <- sprintf("%d_%s_%s", year, band, album)
     path_file <- file.path(band, base_stub)
     
     # Read tracklist
     df <- read.csv(file = paste0(path_file, ".txt"), stringsAsFactors = FALSE)
     
     # Preallocate output list
     n <- nrow(df)
     dat <- vector("list", n)
     
     # Process each song
     for (j in seq_len(n)) {
          song <- df$song[j]
          
     cat(sprintf("Working on song %d: %s | %s | %s | %s\n",
                      j, year, band, album, song))
          
          dat[[j]] <- process_song(band, album, song)
     }
     
     # Save
     save(dat, file = paste0(path_file, ".RData"))    
}