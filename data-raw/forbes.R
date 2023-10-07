## code to prepare `forbes` dataset goes here


forbes <- read.csv("./data-raw/forbes.csv", row.names = 1)

usethis::use_data(forbes, overwrite = TRUE)
