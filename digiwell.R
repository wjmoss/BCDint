library(dplyr)


# load the csv file
df <- read.csv("4HAIE.csv", header=T, sep=";")
data <- df %>% 
  select(-c(34:42)) %>%
  filter(!is.na(qol_1)) %>%
  mutate_at(vars(9,12), ~ as.POSIXct(. , format = "%m/%d/%Y %H:%M:%S"))


# check if data_expiration_start - date_received >= -1  (9,12,24)
diff = (as.Date(as.POSIXct(data[,9], format = "%m/%d/%Y %H:%M:%S")) - 
  as.Date(as.POSIXct(data[,12], format = "%m/%d/%Y %H:%M:%S")))

as.Date(as.POSIXct(df[,9], format = "%m/%d/%Y %H:%M:%S"))[901]
df
