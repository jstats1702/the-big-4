# Setup
rm(list = ls(all.names = TRUE))
setwd("~/Dropbox/PAPERS/projects/MFDA")

source("helper_functions.R")
     
# Inputs
band  <- "anthrax"

year_data  <- c(
     1984,
     1985,
     1987,
     1988,
     1990,
     1993,
     1995,
     1998,
     2003,
     2011,
     2016
)

album_data <- c(
     "fistful_of_metal",
     "spreading_the_disease",
     "among_the_living",
     "state_of_euphoria",
     "persistence_of_time",
     "sound_of_white_noise",
     "stomp_442",
     "volume_8_the_threat_is_real",
     "weve_come_for_you_all",
     "worship_music",
     "for_all_kings"
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
   
     