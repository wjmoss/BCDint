library(dplyr)
library(tidyverse)
library(ggplot2)


# load the csv file
df <- read.csv("4HAIE.csv", header=T, sep=";")
data <- df %>% 
  select(-c(34:42)) %>%
  filter(!is.na(qol_1)) %>%
  mutate_at(vars(9,12), ~ as.POSIXct(. , format = "%m/%d/%Y %H:%M:%S")) %>%
  mutate(date = as.Date(date))


# check if data_expiration_start - date_sent >= -1  (9,12,24)
# diff = as.Date(data[,9]) - as.Date(data[,11])


df1 = data %>% 
  filter(id == "501") %>%
  arrange(date)

df_long <- df1 %>%
  pivot_longer(cols = c(13:16,19:22), names_to = "variable", values_to = "value")

ggplot(df_long, aes(x = date, y = value, color = variable)) +
  geom_line() +
  theme_minimal()

ggplot(df1, aes(x=date, y=steps_wake_daysum)) +
  geom_line() + 
  labs(x = "Time", y = "Value") +
  theme_minimal()


df_summary <- data %>%
  group_by(date) %>%
  summarise(
    mean_y = mean(qol_1),
    lower = quantile(qol_1, 0.25),
    upper = quantile(qol_1, 0.75)
  )

ggplot(df_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean_y), color = "blue") +
  theme_minimal()

tmp=data %>% filter(date==data$date[1])

df_summary <- data %>%
  filter(!is.na(activity)) %>%
  mutate(activity=factor(activity)) %>%
  group_by(date, activity) %>%
  summarise(
    mean_y = mean(qol_1, na.rm = TRUE),
    lower = quantile(qol_1, 0.25, na.rm = TRUE),
    upper = quantile(qol_1, 0.75, na.rm = TRUE)
  )

ggplot(df_summary, aes(x = date, y = mean_y, color = activity)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = activity), alpha = 0.2) +
  theme_minimal()

ggplot(df_summary, aes(x = date, y = mean_y, color = activity, group = activity)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = activity), alpha = 0.2, color = NA) +
  theme_minimal()



df_summary <- data %>%
  filter(!is.na(activity)) %>%
  mutate(activity=factor(activity)) %>%
  group_by(date, activity) %>%
  summarise(
    mean_y = mean(steps_wake_daysum, na.rm = TRUE),
    lower = quantile(steps_wake_daysum, 0.25, na.rm = TRUE),
    upper = quantile(steps_wake_daysum, 0.75, na.rm = TRUE)
  )

ggplot(df_summary, aes(x = date, y = mean_y, color = activity, group = activity)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = activity), alpha = 0.2, color = NA) +
  theme_minimal()


cor(data$weartime_day_hr, data$weartime_day_steps, method = "spearman", use = "complete.obs")

cor(data$qol_1, data$steps_wake_daysum, method = "spearman", use = "complete.obs")

library(changepoint)
cpt <- cpt.mean((df_summary%>% filter(!is.na(mean_y)))$mean_y, method = "PELT")
plot(cpt)


data1 <- data %>% filter(id=="501")
cor(data1$qol_1, data1$qol_2, method = "spearman", use = "complete.obs")



ggplot(df_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean_y), color = "blue") +
  theme_minimal()


df_long <- data %>%
  pivot_longer(cols = c(13:16,19:22), names_to = "variable", values_to = "value")

ggplot(df_long, aes(x = date, y = value, color = variable)) +
  geom_line() +
  theme_minimal()

ggplot(df_long, aes(x=date, y=steps_wake_daysum)) +
  geom_line() + 
  labs(x = "Time", y = "Value") +
  theme_minimal()



