rm(list = ls(all = TRUE)) ## efface les données
setwd('~/thib/projects/MG6/')
#renv::init()
# ---------------------------------------------------
# Libraries
# ---------------------------------------------------


library("svglite")
library("dplyr")
library("brms")
library("sjPlot")
library("lme4")
library('ggplot2')

# ---------------------------------------------------
## * Load and clean data
# ---------------------------------------------------


data <- read.csv('./data/mg6_all.csv', sep = ';')
## identify outliers (accuracy < .55)
## There are 3 outliers
data_full <- data
d <- data_full %>% 
    group_by(subject_id) %>%
    summarise(acc = mean(accuracy_gabor)) %>%
    filter(acc <.55)
outlier <- unique(d$subject_id)

## Remove outliers (3) and  "errors"
data <- data %>%
        filter(!(subject_id %in% outlier), conf %in% c("O", "N")) %>%
        mutate(accuracy2 = ifelse((accuracy_gabor == 1 & meta_evaluation == "C") | (accuracy_gabor == 0 & meta_evaluation == "I"), 1, 0)) %>%
        mutate(accuracy2 = ifelse(view_again == 1, accuracy2, NA)) %>%
        mutate(meta_accuracy = ifelse((view_again == 0 & accuracy_gabor == 1) | (view_again == 1 & accuracy_gabor == 0), 1, 0))
## 2.17% of errors (after exclusion of outliers) 
(nrow(data_full %>% filter(!(subject_id %in% outlier)))- nrow(data))/nrow(data_full %>% filter(!(subject_id %in% outlier)))
    
# ---------------------------------------------------
# define variables and  contrasts
# ---------------------------------------------------

data$cond_order = as.factor(data$cond_order)
data$rt_gabor_centered <- data$rt_gabor - mean(data$rt_gabor, na.rm = TRUE)
data$dt <-  data$OOZ_time
data$dt_centered <- data$OOZ_time - mean(data$OOZ_time, na.rm = TRUE)

data$acc_num <- data$accuracy_gabor
data$accuracy_gabor <- as.factor(data$accuracy_gabor)
contrasts(data$accuracy_gabor) <- - contr.sum(2) ## erreur: -1; correct: 1

data$size <- as.factor(data$size)
contrasts(data$size) <-  -contr.sum(2) ## big: -1; small: 1

data$position <- as.factor(data$circle_pos)
contrasts(data$position) <-  contr.sum(3)

data <- data %>%  mutate(view_again_num  =  ifelse(conf == 'O', 1, 0))
data$view_again <- as.factor(data$view_again)
contrasts(data$view_again) <-  -contr.sum(2) ## O=1, N=-1

data$acc2_num <- data$accuracy2
data$accuracy2 <- as.factor(data$accuracy2)
contrasts(data$accuracy2) <- -contr.sum(2) ## C=1, I=-1


data$meta_acc_num <- data$meta_accuracy
data$meta_accuracy <- as.factor(data$meta_accuracy)
contrasts(data$meta_accuracy) <- -contr.sum(2) ## C=1, I=-1

# ---------------------------------------------------
## * DESCRIPTIVE DATA
# ---------------------------------------------------

## rt, accuracy, conf mean by subject

d <- data  %>%
    group_by(size) %>%
    summarise(
        mean_rt = mean(rt_gabor),
        sd_rt = sd(rt_gabor),
        mean_dt = mean(OOZ_time),
        sd_dt = sd(OOZ_time),
        mean_accuracy = mean(acc_num), 
        sd_accuracy = sd(acc_num), 
        mean_view_again = mean(view_again_num),
        sd_view_again = sd(view_again_num),
        mean_meta_accuracy = mean(meta_acc_num),
        sd_meta_accuracy = sd(meta_acc_num))    
tab_df(d, file = './tables/descriptive.html')

# ---------------------------------------------------
## * fit confirmation time
# ---------------------------------------------------

fit.confirmation_time  <- brm(
        rt_gabor * 1000 ~ accuracy_gabor * size * position +
                (1 + accuracy_gabor * size * position | subject_id),
        family = exgaussian(link = "identity"),
        data = data,
        cores = 4, chains = 4,
        control = list(adapt_delta = .95, max_treedepth = 12),
        iter = 4000, warmup = 1000, seed = 123,
        save_pars = save_pars(all = TRUE)
)
saveRDS(fit.confirmation_time, "./results/fit_confirmation_time.rds")
tab_model(fit.confirmation_time, file = "./tables/fit_confirmation_time.html")

# ---------------------------------------------------
# Figure 2B
# ---------------------------------------------------
fit.confirmation_time <- readRDS("./results/fit_confirmation_time.rds")
size <-  conditional_effects(fit.confirmation_time, "size")

plot <- plot(size, plot = FALSE)[[1]] +
        labs(y = "Confirmation Time (ms)") +
        theme(text = element_text(size = 18))
ggsave('./plots/confirmation_time_size_2B.svg', plot)

# ---------------------------------------------------
# fit commitmment time 
# ---------------------------------------------------

fit.commitment_time <- brm(OOZ_time*1000 ~  accuracy_gabor * size * position  +
                     (1 + accuracy_gabor * size * position |subject_id) ,
           family=exgaussian(link="identity"),
           data = data,
           cores = 4, chains = 4,
           control = list(adapt_delta = .95,  max_treedepth = 12),
           iter = 4000,  warmup = 1000, seed = 123,
           save_model = 'dt.stan',
           save_pars = save_pars(all = TRUE)
           )
saveRDS(fit.commitment_time, "./results/fit_commitment_time.rds")
tab_model(fit.commitment_time, file = "./tables/fit_commitment_time.html")


# ---------------------------------------------------
## * Fit accuracy
# ---------------------------------------------------


fit_acc <- brm(accuracy_gabor  ~ size * position * dt_centered + (1 + size * position * dt_centered | subject_id),
               family = bernoulli(link = "logit"),
           data = data,
           cores = 4, chains = 4,
           control = list(adapt_delta = .9,  max_treedepth = 12),
           iter = 3000,  warmup = 1000, seed = 123,
           )
saveRDS(fit_acc, "./results/fit_acc.rds")
tab_model(fit_acc, file = "./tables/fit_accuracy.html")
         -0.11699656 0.05544853 -0.22540972 -0.008194695

# ---------------------------------------------------
# Figure 2A
# ---------------------------------------------------
fit_acc <- readRDS("./results/fit_acc.rds")

size <-  conditional_effects(fit_acc, "size")
plot <- plot(size, plot = FALSE)[[1]] +
    labs(y = "Accuracy") +
    theme(text = element_text(size = 18))        
        
ggsave('./plots/acc_size_2A.svg', plot)

# ---------------------------------------------------
## fit View again
#
# ---------------------------------------------------

fit.conf <- brm(view_again ~ accuracy_gabor * size * position * rt_gabor_centered + (1 + accuracy_gabor * size * position * rt_gabor_centered | subject_id),
        family = bernoulli(link = "logit"),
        data = data,
        init = 0,
        cores = 4, chains = 4,
        control = list(adapt_delta = .9, max_treedepth = 12),
        iter = 3000, warmup = 1000, seed = 123,
        save_pars = save_pars(all = TRUE)
)


saveRDS(fit.conf, "./results/fit_conf.rds")
tab_model(fit.conf, file = "./tables/fit_conf.html")

# ---------------------------------------------------
# Figure 2C
# ---------------------------------------------------

fit.conf <- readRDS("./results/fit_conf.rds")

acc_size <-  conditional_effects(fit.conf, "accuracy_gabor:size")
plot <- plot(acc_size, plot = FALSE)[[1]] +
    labs(y = "View again", x = "Accuracy") +
    scale_x_discrete(labels = c("error", "correct")) +
    theme(text = element_text(size = 18))        
ggsave('./plots/conf_acc_size_2C.svg', plot)

# ---------------------------------------------------
# Figure 2D
# ---------------------------------------------------

size_rt <-  conditional_effects(fit.conf, "rt_gabor_centered:size")
mean_rt <- mean(data$rt_gabor, na.rm = TRUE)
size_rt[[1]]$rt_gabor <- size_rt[[1]]$rt_gabor_centered + mean_rt
df <- size_rt[[1]]
plot <- ggplot(df, aes(x = rt_gabor, y = estimate__, color = size)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = size),
                alpha = 0.2, color = NA) +
        labs(x = "Confirmation time (s)", y = "View again") +
    theme(text = element_text(size = 18))          
ggsave("./plots/conf_size_rt_2D.svg", plot)

# ---------------------------------------------------
#  Figure S1
# ---------------------------------------------------

acc_rt <- conditional_effects(fit.conf, "rt_gabor_centered:accuracy_gabor")
acc_rt[[1]]$rt_gabor <- acc_rt[[1]]$rt_gabor_centered + mean_rt
df <- acc_rt[[1]]
plot <- ggplot(df, aes(
        x = rt_gabor, y = estimate__,
        color = factor(accuracy_gabor, levels = c(0,1), labels = c("error","correct")),
        fill  = factor(accuracy_gabor, levels = c(0,1), labels = c("error","correct")),
        group = accuracy_gabor
)) +
        geom_line(linewidth = 1) +
        geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, color = NA) +
        labs(
            x = "Confirmation time (s)",
            y = "View again",
            color = "Accuracy",
            fill = "Accuracy"
        ) +
    theme(text = element_text(size = 18))          
plot
ggsave('./plots/conf_acc_rt_S1.svg', plot)


# ---------------------------------------------------
## Additional analysis Revision I
## Effect of see_again on accuracy
# ---------------------------------------------------

## Compute final_accuracy : accuracy of the revised decision
dd <- data %>% filter(view_again == 1) %>% ## select trials for which view_again
        mutate(final_accuracy = ifelse((accuracy_gabor == 1 & meta_evaluation == "C") | (accuracy_gabor == 0 & meta_evaluation == "I"), 1, 0))

dd <- dd %>%
        group_by(subject_id) %>%
        summarise(accuracy = mean(as.numeric(accuracy_gabor)-1), final_accuracy = mean(as.numeric(final_accuracy)))
dd_summary <- dd %>% summarise(accuracy = mean(accuracy), final_accuracy = mean(final_accuracy))

##   accuracy final_accuracy
##      <dbl>          <dbl>
## 1    0.609          0.739


t.test(dd$accuracy, dd$final_accuracy, paired = TRUE)
## data:  dd$accuracy and dd$final_accuracy
## t = -5.9454, df = 23, p-value = 4.625e-06
## alternative hypothesis: true mean difference is not equal to 0
## 95 percent confidence interval:
##  -0.1749834 -0.0846466
## sample estimates:
## mean difference 
##       -0.129815 

# ---------------------------------------------------
## Additional analysis Revision II
## Reported errors
# ---------------------------------------------------

dd <- data_full %>%
    filter(!(subject_id %in% outlier)) %>%
    mutate(error = ifelse(conf == 0, 1, 0))
table(dd$error, dd$size)
##   big small
## 0 2246  2262
## 1   57    43
ll <- lm(error ~ accuracy_gabor*size*rt_gabor, data = dd)
summary(ll)

## Coefficients:
##                 Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     0.044139   0.004504   9.800  < 2e-16 ***
## accuracy_gabor -0.027144   0.004671  -5.811 6.62e-09 ***
## sizesmall      -0.006854   0.004280  -1.601    0.109    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.1452 on 4605 degrees of freedom
## Multiple R-squared:  0.007714,	Adjusted R-squared:  0.007283 
## F-statistic:  17.9 on 2 and 4605 DF,  p-value: 1.804e-08

## remove outliers and E trials
## accuracy2: whereas 2nd accuracy judg. is correct
## meta_accuracy: whereas it was correct to ask a 2nd stim



