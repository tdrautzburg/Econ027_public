
pacman::p_load(haven, readr)
oj_df <- as.data.frame(FrozenJuice)
oj_df$time_index <- as.numeric(time(FrozenJuice))
oj_df <- janitor::clean_names(oj_df)
write_dta(oj_df, here("..","Data", "FrozenJuice.dta"))
write.csv(oj_df, here("..", "Data", "FrozenJuice.csv"), row.names = FALSE)