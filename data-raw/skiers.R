
skiers <- data.frame(
  cost = c(32, 61, 59, 36, 62),
  lift = c(64, 37, 40, 62, 46),
  depth = c(65, 62, 45, 34, 43),
  powder = c(67, 65, 43, 35, 40)
)

usethis::use_data(skiers, overwrite = TRUE)
