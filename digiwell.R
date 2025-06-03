library(dplyr)
library(tidyverse)
library(ggplot2)
library(MASS)
#library(bestNormalize)
source("ricf_int.R")
source("ricf_dg.R")


## helper functions
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

repair_50_bias <- function(x, jitter_sd = 5) {
  # Drop 0 and 100
  x <- x[!(x %in% c(0, 100))]
  
  # Jitter the values that are exactly 50
  idx_50 <- which(x == 50)
  n_50 <- length(idx_50)
  
  if (n_50 > 0) {
    # Replace 50 with values jittered around 50 (e.g., N(50, 3))
    x[idx_50] <- (rnorm(n_50, mean = 50, sd = jitter_sd))
    # Ensure jittered values stay in bounds
    x[idx_50] <- pmin(pmax(x[idx_50], 1), 99)
  }
  
  return(x)
}

repair_bias <- function(x, jitter_sd = 5) {
  # Drop 0 and 100
  x <- x[!(x %in% c(0, 100))]
  
  # Jitter the values that are between 48 to 52
  #is_heaped <- x %in% 48:52
  #x_jittered <- x  
  #x_jittered[is_heaped]  <- rnorm(sum(is_heaped), mean = 50, sd = jitter_sd)
  #x_jittered[is_heaped] <- runif(sum(is_heaped), min = 47.5, max = 52.5)
  res = sapply(x, \(t) t + rnorm(1, mean=0, sd=jitter_sd) * (t <=52) * (t>=48) )
  
  return(res)
}


## load the csv file
# filter non-NA rows and transform time format
df <- read.csv("data/4HAIE.csv", header=T, sep=";")
#na_check <- c("qol_1", "qol_2", "stresd_1", "res_1", "Sleep5")
na_check <- names(df)[c(13:23)]
data <- df %>% 
  dplyr::select(-(34:42)) %>%
  #filter(if_all(all_of(na_check), ~ !is.na(.))) %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  filter(Sleep5 < 100 & Sleep5 > 0) %>%
  filter(steps_wake_daysum > 0) %>%
  filter(qol_1 < 100 & qol_1 > 0) %>%
  filter(qol_2 < 100 & qol_2 > 0) %>%
  mutate(qol_1 = repair_bias(qol_1, jitter_sd = 2),
         qol_2 = repair_bias(qol_2, jitter_sd = 2),
         Sleep5 = repair_bias(Sleep5, jitter_sd = 2)) %>%
  #filter(if_all(all_of(na_check), ~ !is.na(.))) %>%
  mutate_at(vars(9,12), ~ as.POSIXct(. , format = "%m/%d/%Y %H:%M:%S")) %>%
  mutate(date = as.Date(date)) %>%
  mutate(qol = (qol_1 + qol_2) / 2)
  #mutate(qol_1_t=qnorm((qol_1+1)/102))

## histograms of main variables
ggplot(data, aes(x = qol_1)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

ggplot(data, aes(x = qol_2)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

ggplot(data, aes(x = qol)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

ggplot(data, aes(x = steps_wake_daysum)) +
  geom_histogram(binwidth = 100, fill = "steelblue", color = "black") +
  theme_minimal()

#ggplot(data%>%filter(abs(qol-100)<2), aes(x = steps_t)) +
#  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
#  theme_minimal()

# other variables
ggplot(data, aes(x = Sleep5)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

data <- data %>% mutate(activity=factor(activity))
ggplot(data, aes(x = qol_t, fill = activity)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  theme_minimal()

ggplot(data, aes(x = steps_t, fill = activity)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  theme_minimal()


d1_o = filter(data_obs, activity==1)
d2_o = filter(data_obs, activity==2)
cor(d1_o$qol, d1_o$steps_t)
cor(d2_o$qol, d2_o$steps_t)
cor(d1_o$Intention, d1_o$steps_t)
cor(d2_o$Intention, d2_o$steps_t)
sum(d1_o$Intention) / nrow(d1_o)
sum(d2_o$Intention) / nrow(d2_o)

d1_i = filter(data_int, activity==1)
d2_i = filter(data_int, activity==2)
cor(d1_i$qol, d1_i$steps_t)
cor(d2_i$qol, d2_i$steps_t)
cor(d1_i$Intention, d1_i$steps_t)
cor(d2_i$Intention, d2_i$steps_t)
sum(d1_i$Intention) / nrow(d1_i)
sum(d2_i$Intention) / nrow(d2_i)

#### cyclic model

## tranform qol
# conflict with bestNormalize
bc <- boxcox(data$qol ~ 1)  # recommends lambda
lambda <- bc$x[which.max(bc$y)]
data = data %>% mutate(qol_t = (qol^lambda - 1) / lambda)
#x_trans <- if (lambda == 0) log(x) else (x^lambda - 1) / lambda
ggplot(data, aes(x = qol_t)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

# for multiple cols
#df <- df %>%
#  mutate(across(c(x1, x2), ~ qnorm((rank(.) - 0.5) / length(.))))

ecdf_qol <- ecdf(data$qol)
eps <- 1e-6
probs <- ecdf_qol(data$qol)
probs <- pmin(pmax(probs, eps), 1 - eps)

data = data %>% 
  mutate(qol_t = (\(x) qnorm(pmin(pmax(ecdf_qol(x), eps), 1 - eps))) (qol))
ggplot(data, aes(x = qol_t)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black") +
  theme_minimal()


normalizer_qol = bestNormalize(data$qol)
data = data %>% mutate(qol_t = normalizer_qol$x.t)
normalizer_qol$chosen_transform
attributes(normalizer_qol)
class(normalizer_qol)
methods(class = class(normalizer_qol))

# transform steps
bc <- boxcox(data$steps_wake_daysum ~ 1)  # recommends lambda
lambda <- bc$x[which.max(bc$y)]
data = data %>% 
  mutate(steps_t = (steps_wake_daysum^lambda - 1) / lambda) #%>%
  #mutate(steps_t = scale(steps_t))
ggplot(data, aes(x = steps_t)) +
  geom_histogram(binwidth = (max(x)-min(x))/100, fill = "steelblue", color = "black") +
  theme_minimal()
ggplot(data, aes(x = steps_t)) +
  geom_histogram(binwidth = 10, fill = "steelblue", color = "black") +
  theme_minimal()


normalizer_steps = bestNormalize(data$steps_wake_daysum)
normalizer_steps$method
data = data %>% mutate(steps_t = normalizer_steps$x.t)


# transform sleep (optional)
bc <- boxcox(data$Sleep5 ~ 1)  # recommends lambda
lambda <- bc$x[which.max(bc$y)]
data = data %>% mutate(sleep_t = (Sleep5^lambda - 1) / lambda)
#x_trans <- if (lambda == 0) log(x) else (x^lambda - 1) / lambda
ggplot(data, aes(x = sleep_t)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  theme_minimal()

data = data %>% 
  mutate(sleep_t = (\(x) qnorm(pmin(pmax(ecdf_qol(x), eps), 1 - eps))) (Sleep5))
ggplot(data, aes(x = sleep_t)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black") +
  theme_minimal()

normalizer_sleep = bestNormalize(data$Sleep5)
normalizer_sleep$method
data = data %>% mutate(sleep_t = normalizer_sleep$x.t)


## data before and after the covid lockdown
#data_before <- data %>% filter(date < as.Date("2020-03-12"))
data_int1 <- data %>% filter(date >= as.Date("2020-03-12") & date <= as.Date("2020-04-30"))
data_int2 <- data %>% filter(date >= as.Date("2020-10-22") & date <= as.Date("2020-11-22"))
data_int <- rbind(data_int1, data_int2)
data_obs <- data %>% 
  filter(date < as.Date("2020-03-12") | date > as.Date("2020-04-30")) %>%
  filter(date < as.Date("2020-10-22") | date > as.Date("2020-11-20"))
#12.03 - 30.04
#22.10 - 20.11
dim(data_obs)
dim(data_int)
data <- data %>% 
  mutate(ifint = (date>=as.Date("2020-03-12") & date <= as.Date("2020-04-30")) | 
           (date>=as.Date("2020-10-22") & date <= as.Date("2020-11-20"))  ) %>%
  mutate(ifint=factor(ifint))
ggplot(data, aes(x = steps_t, fill = ifint)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  theme_minimal() 
t.test(data_obs$steps_t, data_int$steps_t )

## steps -> qol_t -> steps, steps also influenced by lockdown
Y_o <- dplyr::select(data_before, c(qol_t, steps_t))
Y_i <- dplyr::select(data_after, c(qol_t, steps_t))
targets = list(numeric(0), c(2))
target.length <- c(nrow(Y_o), nrow(Y_i))
cat(dim(Y_o), dim(Y_i))
res <- ricf_int_(L = matrix(c(0,1,1,0),2,2), 
                 data = as.matrix(rbind(Y_o,Y_i)), 
                 targets = targets,
                 target.length = target.length)
res


## model with sleep
# to do: use 1-lag for sleep?
# sleep -> qol <-> steps
Y_o <- dplyr::select(data_obs, c(qol_t, steps_t, sleep_t))
Y_i <- dplyr::select(data_int, c(qol_t, steps_t, sleep_t))
targets = list(numeric(0), c(2))
target.length <- c(nrow(Y_o), nrow(Y_i))
cat(dim(Y_o), dim(Y_i))
L <- matrix(0, 3, 3)
L[1,2] = L[2,1] = L[3,1] = 1
res <- ricf_int_(L = L,
                 data = as.matrix(rbind(Y_o,Y_i)), 
                 targets = targets,
                 target.length = target.length)
res
#In t2 * solve(S[pa, pa], a) :
#Recycling array of length 1 in array-vector arithmetic is deprecated.
#Use c() or as.vector() instead.

## CI




## model as observational data
Y <- dplyr::select(data, c(qol_t, steps_t, sleep_t))
L <- matrix(0, 3, 3)
L[1,2] = L[2,1] = L[3,1] = 1
res <- ricf_int_(L = L, data = as.matrix(Y))
res

res1 <- ricf_int_(L = L, data = as.matrix(Y_o))
res2 <- ricf_int_(L = L, data = as.matrix(Y_i))
res1$Lambdahat
res2$Lambdahat









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



