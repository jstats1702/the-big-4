# Setup
rm(list = ls(all.names = TRUE))
setwd("~/Dropbox/PAPERS/projects/MFDA")

source("helper_functions.R")
     
# Inputs
band  <- "slayer"

year_data <- c(
     1983,
     1985,
     1986,
     1988,
     1990,
     1994,
     1998,
     2001,
     2006,
     2009,
     2015
)

album_data <- c(
     "show_no_mercy",
     "hell_awaits",
     "reign_in_blood",
     "south_of_heaven",
     "seasons_in_the_abyss",
     "divine_intervention",
     "diabolus_in_musica",
     "god_hates_us_all",
     "christ_illusion",
     "world_painted_blood",
     "repentless"
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
   
     