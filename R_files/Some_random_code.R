disease_df <- disease_df %>%
  mutate(x1 = 0.5*cases, x2 = 0.2*cases)
write.csv(disease_df, "HepatitisC2.csv", row.names = FALSE)
