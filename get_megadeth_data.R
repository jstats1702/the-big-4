# Setup
rm(list = ls(all.names = TRUE))
setwd("~/Dropbox/PAPERS/projects/MFDA")

source("helper_functions.R")
     
# Inputs
band  <- "megadeth"

year_data  <- c(
     1985,
     1986,
     1988,
     1990,
     1992,
     1994,
     1997,
     1999,
     2001,
     2004,
     2007,
     2009,
     2011,
     2013,
     2016,
     2022
)

album_data <- c(
     "killing_is_my_business_and_business_is_good",
     "peace_sells_but_whos_buying",
     "so_far_so_good_so_what",
     "rust_in_peace",
     "countdown_to_extinction",
     "youthanasia",
     "cryptic_writings",
     "risk",
     "the_world_needs_a_hero",
     "the_system_has_failed",
     "united_abominations",
     "endgame",
     "th1rt3en",
     "super_collider",
     "dystopia",
     "the_sick_the_dying_and_the_dead"
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
   
     