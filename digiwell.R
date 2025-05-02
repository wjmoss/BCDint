library(dplyr)
library(tidyverse)
library(ggplot2)


## load the csv file
# filter non-NA rows and transform time format
df <- read.csv("4HAIE.csv", header=T, sep=";")
data <- df %>% 
  select(-c(34:42)) %>%
  filter(if_all(all_of(c("qol_1", "qol_2", "Sleep5")), ~ !is.na(.))) %>%
  mutate_at(vars(9,12), ~ as.POSIXct(. , format = "%m/%d/%Y %H:%M:%S")) %>%
  mutate(date = as.Date(date))

## histograms of some variables
ggplot(data, aes(x = qol_1)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

ggplot(data, aes(x = qol_2)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

ggplot(data, aes(x = steps_wake_daysum)) +
  geom_histogram(binwidth = 100, fill = "steelblue", color = "black") +
  theme_minimal()

ggplot(data, aes(x = Sleep5)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

data <- data %>% mutate(activity=factor(activity))
ggplot(data, aes(x = qol_1, fill = activity)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  theme_minimal()

ggplot(data, aes(x = qol_2, fill = activity)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  theme_minimal()


## visualization for one id
df1 = data %>% 
  filter(id == "501") %>%
  arrange(date)

df_long <- df1 %>%
  pivot_longer(cols = c(13:16,19:22), names_to = "variable", values_to = "value")

# variables in range [0,100]
ggplot(df_long, aes(x = date, y = value, color = variable)) +
  geom_point() +
  theme_minimal()

# steps
ggplot(df1, aes(x=date, y=steps_wake_daysum)) +
  geom_point() + 
  labs(x = "Time", y = "Value") +
  theme_minimal()



## visualization of mean of each date
df_summary <- data %>%
  group_by(date) %>%
  summarise(
    mean_y = mean(qol_1, na.rm=T),
    lower = quantile(qol_1, 0.25, na.rm=T),
    upper = quantile(qol_1, 0.75, na.rm=T)
  )


ggplot(df_summary, aes(x = date)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean_y), color = "blue") +
  geom_vline(xintercept = as.Date("2020-03-12"), color = "red", linetype = "solid", size = 1) +
  annotate("text",
           x = as.Date("2020-03-12"),
           y = min(df_summary$lower, na.rm = TRUE) + 0.05 * diff(range(df_summary$lower, df_summary$upper, na.rm = TRUE)),
           label = "2020-03-12",
           color = "red",
           angle = 0,
           vjust = -0.5,
           hjust = -0.2) +
  theme_minimal()

#tmp=data %>% filter(date==data$date[1])

# mean plot, with group 1 = active and 2 = non-active
df_summary <- data %>%
  filter(!is.na(activity)) %>%
  mutate(activity=factor(activity)) %>%
  group_by(date, activity) %>%
  summarise(
    mean_y = mean(qol_1, na.rm = TRUE),
    lower = quantile(qol_1, 0.25, na.rm = TRUE),
    upper = quantile(qol_1, 0.75, na.rm = TRUE)
  )

#ggplot(df_summary, aes(x = date, y = mean_y, color = activity)) +
#  geom_line() +
#  geom_ribbon(aes(ymin = lower, ymax = upper, fill = activity), alpha = 0.2) +
#  theme_minimal()

ggplot(df_summary, aes(x = date, y = mean_y, color = activity, group = activity)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = activity), alpha = 0.2, color = NA) +
  theme_minimal()


## steps, group 1 = active and 2 = non-active
df_summary <- data %>%
  filter(!is.na(activity)) %>%
  mutate(activity=factor(activity)) %>%
  group_by(date, activity) %>%
  summarise(
    mean_steps = mean(steps_wake_daysum, na.rm = TRUE),
    lower = quantile(steps_wake_daysum, 0.25, na.rm = TRUE),
    upper = quantile(steps_wake_daysum, 0.75, na.rm = TRUE)
  )

ggplot(df_summary, aes(x = date, y = mean_steps, color = activity, group = activity)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = activity), alpha = 0.2, color = NA) +
  geom_vline(xintercept = as.Date("2020-03-12"), color = "blue", linetype = "solid", size = 1) +
  annotate("text",
           x = as.Date("2020-03-12"),
           y = min(df_summary$lower, na.rm = TRUE) + 0.05 * diff(range(df_summary$lower, df_summary$upper, na.rm = TRUE)),
           label = "2020-03-12",
           color = "blue",
           angle = 0,
           vjust = -0.5,
           hjust = -0.2) +
  theme_minimal()


















cor(data$stresd_1, data$steps_wake_daysum, method = "spearman", use = "complete.obs")



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



