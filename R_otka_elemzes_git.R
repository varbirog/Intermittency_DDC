## Libraries--------
library(readxl)
library(tidyverse)
library(ggpubr)
library(ade4)
library(cowplot)
library(gt)

# import excel--------
# # Read data
setwd("C:/Users/VGabor/Documents/R/R_Otka")
# import excel "otka_final.xlsx" spec table
source("Suppl_function.R")
spec_table <- read_excel("otka_final.xlsx", sheet = "species")
# import excel "otka_final.xlsx" spec table
dry <- read_excel("otka_final.xlsx", sheet = "dry")
sites <- read_excel("otka_final.xlsx", sheet = "sites")
data <- read_excel("otka_final.xlsx")

# Function sub.diag.R
sub.diag <- function(mat, method="gaining.turnover", relativize="J")
  #
  # Compute directional.response(mat); extract sub-diagonal of output matrices
{
  x <- directional.response(mat, method=method, relativize=relativize)
  n.sd <- nrow(x$mat.out) # Number of values in sub-diagonal of mat.out
  tmp <- cbind( 2:n.sd, 1:(n.sd-1) )
  # print(tmp)
  s.diag <- x$mat.out[tmp]
  if(is.na(x$total.t[[1]])) tt <- NA else tt <- x$total.t[tmp]
  if(is.na(x$total.n[[1]])) tn <- NA else tn <- x$total.n[tmp]
  if(is.na(x$total.strict.n[[1]])) tn2 <- NA else tn2 <- x$total.strict.n[tmp]
  #
  list(sub.diag=s.diag, total.t=tt, total.n=tn, total.strict.n=tn2)
}





spec_table$year_season=factor(spec_table$year_season, levels=c(
  "20NY",
  "20Ő",
  "21TÉL",
  "21TAV",
  "21NY",
  "21Ő",
  "22TÉL",
  "22TAV",
  "22NY",
  "22Ő",
  "23TÉL",
  "23TAV"
  
),ordered=TRUE)
spec_table$campaign=factor(spec_table$campaign, levels=c(
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  14,
  15,
  16
  
  
),ordered=TRUE)



# nested minta  year_season----------------
minta_nest= spec_table %>% 
  select(site,year_season,taxon,db)%>% 
  complete(site, year_season, fill = list(db = 0)) %>%
  group_by(site) %>% 
  nest() %>%
  mutate(table=map(data, ~.x %>% pivot_wider(names_from = taxon, values_from = db,values_fill = 0))) %>%
  mutate(kodes=map(table, ~.x[,1])) %>%
  mutate(kodes=map(kodes, ~as.data.frame(.x))) %>%
  mutate(kodex=map(kodes, ~.x[-1,])) %>%
  mutate(table=map(table, ~.x[,-1])) %>%
  mutate(table=map(table, ~ifelse(.x>0, 1, 0))) %>% 
  mutate(table=map(table, ~as.data.frame(.x)))


# nested minta  campaign ---------------
minta_nest_c= spec_table %>% 
  select(site,campaign,taxon,db)%>% 
  complete(site, campaign, fill = list(db = 0)) %>%
  group_by(site) %>% 
  nest() %>%
  mutate(table=map(data, ~.x %>% pivot_wider(names_from = taxon, values_from = db,values_fill = 0))) %>%
  mutate(kodes=map(table, ~.x[,1])) %>%
  mutate(kodes=map(kodes, ~as.data.frame(.x))) %>%
  mutate(kodex=map(kodes, ~.x[-1,])) %>%
  mutate(table=map(table, ~.x[,-1])) %>%
  mutate(table=map(table, ~ifelse(.x>0, 1, 0))) %>% 
  mutate(table=map(table, ~as.data.frame(.x)))


#Calculate indices - Jaccard----------

minta_calc3= minta_nest_c %>% 
  mutate(overlap=map(table, ~sub.diag(.x, method="overlap", relativize="J")$sub.diag)) %>% 
  mutate(gain=map(table, ~sub.diag(.x, method="gain", relativize="J")$sub.diag)) %>% 
  mutate(loss=map(table, ~sub.diag(.x, method="loss", relativize="J")$sub.diag)) %>% 
  mutate(gaining.turnover=map(table, ~sub.diag(.x, method="gaining.turnover", relativize="J")$sub.diag)) %>% 
  mutate(neutral.turnover=map(table, ~sub.diag(.x, method="neutral.turnover", relativize="J")$sub.diag)) %>% 
  mutate(losing.turnover=map(table, ~sub.diag(.x, method="losing.turnover", relativize="J")$sub.diag)) %>% 
  mutate(gaining.nestedness=map(table, ~sub.diag(.x, method="gaining.nestedness", relativize="J")$sub.diag)) %>% 
  mutate(neutral.nestedness=map(table, ~sub.diag(.x, method="neutral.nestedness", relativize="J")$sub.diag)) %>% 
  mutate(losing.nestedness=map(table, ~sub.diag(.x, method="losing.nestedness", relativize="J")$sub.diag)) 
#Calculate indices - Sorensen----------
minta_calc4= minta_nest %>% 
  mutate(overlap=map(table, ~sub.diag(.x, method="overlap", relativize="S")$sub.diag)) %>% 
  mutate(gain=map(table, ~sub.diag(.x, method="gain", relativize="S")$sub.diag)) %>% 
  mutate(loss=map(table, ~sub.diag(.x, method="loss", relativize="S")$sub.diag)) %>% 
  mutate(gaining.turnover=map(table, ~sub.diag(.x, method="gaining.turnover", relativize="S")$sub.diag)) %>% 
  mutate(neutral.turnover=map(table, ~sub.diag(.x, method="neutral.turnover", relativize="S")$sub.diag)) %>% 
  mutate(losing.turnover=map(table, ~sub.diag(.x, method="losing.turnover", relativize="S")$sub.diag)) %>% 
  mutate(gaining.nestedness=map(table, ~sub.diag(.x, method="gaining.nestedness", relativize="S")$sub.diag)) %>% 
  mutate(neutral.nestedness=map(table, ~sub.diag(.x, method="neutral.nestedness", relativize="S")$sub.diag)) %>% 
  mutate(losing.nestedness=map(table, ~sub.diag(.x, method="losing.nestedness", relativize="S")$sub.diag)) 

#unnested Jaccard

unnested_minta_J= minta_calc3 %>% 
  select(c(site,overlap,gain,loss,gaining.turnover,neutral.turnover,losing.turnover,gaining.nestedness,neutral.nestedness,losing.nestedness)) %>%
  unnest_longer(c(overlap,gain,loss,gaining.turnover,neutral.turnover,losing.turnover,gaining.nestedness,neutral.nestedness,losing.nestedness))


#unnested Sorensen

unnested_minta_S= minta_calc4 %>% 
  select(c(site,overlap,gain,loss,gaining.turnover,neutral.turnover,losing.turnover,gaining.nestedness,neutral.nestedness,losing.nestedness)) %>%
  unnest_longer(c(overlap,gain,loss,gaining.turnover,neutral.turnover,losing.turnover,gaining.nestedness,neutral.nestedness,losing.nestedness))




unnested_kodes= minta_calc4 %>%  
  select(c(site,kodex))  %>% 
  unnest(kodex) 
  
# Bind cols

dat_s= bind_cols(unnested_minta_S,unnested_kodes[,-1])  

dat_j= bind_cols(unnested_minta_J,unnested_kodes[,-1]) 




# create a matrix for the heatmap not used in the article
heatmap_data <- data.frame(
  id = dat_s$site,
  kodex = dat_s$kodex,
  value = dat_s$losing.turnover
)

# convert to matrix
heatmap_matrix <- as.matrix(heatmap_data)

# Heatmap 
ggplot(heatmap_data, aes(y=id, x=kodex, fill= value)) + 
  geom_tile()



# Figures --------------------------

my_comparisons <- list( c("p", "i"))



#join with sites --------------

jacc=dat_j %>% inner_join(sites,by="site")

sorr=dat_s %>% inner_join(sites,by="site")




# Make a longer datastructure

# longer Jaccard
df_long_jacc <- jacc %>% 
  pivot_longer(cols = c(overlap,gain,loss,
                        gaining.turnover,
                        neutral.turnover,
                        losing.turnover,
                        gaining.nestedness,
                        neutral.nestedness,
                        losing.nestedness,
  ),
  names_to = "metric",
  values_to = "value")


# longer Sorrensen

df_long_sorr <- sorr %>%
  pivot_longer(cols = c(overlap,gain,loss,
                        gaining.turnover,
                        neutral.turnover,
                        losing.turnover,
                        gaining.nestedness,
                        neutral.nestedness,
                        losing.nestedness,
                        ),
               names_to = "metric",
               values_to = "value")

### By sites Permanent---------
df_long_sorr %>% filter(p_i=="p", metric %in% c("gain","loss","overlap")) %>% 
ggplot(., aes(x = kodex, y = value, color = metric)) +
  geom_line(aes(group = interaction(site, metric)), alpha = 0.6) +
  facet_wrap(~site) +
  labs(x = "kodex", y = "Value", title = "Overlap, Gain, Loss by Camp and Site") +
  scale_x_discrete(breaks = unique(df_long_sorr$kodex)) +
  theme_minimal()

df_long_sorr %>% filter(p_i=="p", metric %in% c("gaining.turnover", "neutral.turnover" ,"losing.turnover" )) %>% 
  ggplot(., aes(x = kodex, y = value, color = metric)) +
  geom_line(aes(group = interaction(site, metric)), alpha = 0.6) +
  facet_wrap(~site) +
  labs(x = "kodex", y = "Value", title = "Turnover Overlap, Gain, Loss by Camp and Site") +
  scale_x_discrete(breaks = unique(df_long_sorr$kodex)) +
  theme_minimal()


df_long_sorr %>% filter(p_i=="p", metric %in% c("gaining.nestedness", "neutral.nestedness" ,"losing.nestedness" )) %>% 
  ggplot(., aes(x = kodex, y = value, color = metric)) +
  geom_line(aes(group = interaction(site, metric)), alpha = 0.6) +
  facet_wrap(~site) +
  labs(x = "kodex", y = "Value", title = "Nestedness Overlap, Gain, Loss by Camp and Site") +
  scale_x_discrete(breaks = unique(df_long_sorr$kodex)) +
  theme_minimal()


### By sites Intermittent---------
df_long_sorr %>% filter(p_i=="p", metric %in% c("gain","loss","overlap")) %>% 
  ggplot(., aes(x = kodex, y = value, color = metric)) +
  geom_line(aes(group = interaction(site, metric)), alpha = 0.6) +
  facet_wrap(~site) +
  labs(x = "kodex", y = "Value", title = "Overlap, Gain, Loss by Camp and Site") +
  scale_x_discrete(breaks = unique(df_long_sorr$kodex)) +
  theme_cowplot()

df_long_sorr %>% filter(p_i=="i", metric %in% c("gaining.turnover", "neutral.turnover" ,"losing.turnover" )) %>% 
  ggplot(., aes(x = kodex, y = value, color = metric)) +
  geom_line(aes(group = interaction(site, metric)), alpha = 0.6) +
  facet_wrap(~site) +
  labs(x = "kodex", y = "Value", title = "Turnover Overlap, Gain, Loss by Camp and Site") +
  scale_x_discrete(breaks = unique(df_long_sorr$kodex)) +
  theme_minimal()


df_long_sorr %>% filter(p_i=="i", metric %in% c("gaining.nestedness", "neutral.nestedness" ,"losing.nestedness" )) %>% 
  ggplot(., aes(x = kodex, y = value, color = metric)) +
  geom_line(aes(group = interaction(site, metric)), alpha = 0.6) +
  facet_wrap(~site) +
  labs(x = "kodex", y = "Value", title = "Nestedness Overlap, Gain, Loss by Camp and Site") +
  scale_x_discrete(breaks = unique(df_long_sorr$kodex)) +
  theme_minimal()

# Group by summaries camp, metric, p_i 

df_summary <- df_long_jacc %>%
  group_by(kodex, metric, p_i) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Plot
ggplot() +
  # Halvány egyedi vonalak
  geom_line(data = df_long_jacc, aes(x = as.numeric(kodex), y = value, group = interaction(site, metric), color = metric),
            alpha = 0.1, size = 0.5) +
  # Vastag csoportátlag vonalak p_i szerint
  geom_line(data = df_summary, aes(x = as.numeric(kodex), y = mean_value, color = metric, linetype = p_i),
            size = 1) +
  facet_wrap(~p_i,  scales = "free_y") +
  labs(x = "kodex", y = "Value", title = "Overlap, Gain, Loss with Group Averages by p_i") +
  theme_minimal()

# Seasonal summaries-Jaccard---------------


df_summary <- df_long_jacc %>%
  group_by(kodex, metric, p_i) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")


df_summary$p_i <- factor(df_summary$p_i, levels = c("i", "p"), labels = c("intermittent", "permanent"))

df_summary %>% 
  ggplot(aes(x = kodex, y = mean_value, color = p_i, group = interaction(p_i, metric))) +
  geom_line(size = 0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("#CC6600","#33CC99")) +
  labs(x = "kodex", y = "Value", 
       title = "Jaccard Metric with Group Averages") +
  theme_minimal() +
  facet_wrap(~metric) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(text = element_text(size = 14))


ggplot() +
  
  geom_line(data = df_summary, aes(x = round(as.numeric(kodex),0), y = mean_value, color = metric, linetype = p_i),
            size = 1) +
  facet_wrap(~p_i ,labeller = label_value) +
  labs(x = "kodex", y = "Value", title = "Jaccard Metric with Group Averages") +
  theme_minimal()

# ---- Scientific table with gt -Jaccard-----------

df <- df_long_jacc %>%
  mutate(flow_type = ifelse(p_i == "p", "Permanent", "Intermittent"))

df_summary <- df %>%
  group_by(flow_type, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(n()),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    range = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )


df_summary %>%
  arrange(flow_type, metric) %>%
  gt(groupname_col = "flow_type") %>%
  tab_header(title = "Summary Statistics by Flow Type and Metric") %>%
  fmt_number(columns = where(is.numeric), decimals = 3)

# ---- Visualization ----
ggplot(df, aes(x = metric, y = value, fill = flow_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  #geom_jitter(width = 0.1, alpha = 0.8, shape = 21) +
  labs(title = "Metric Distributions by Flow Type",
       x = "Metric", y = "Value") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Scientific table with gt -Sørensen-----------

df <- df_long_sorr %>%
  mutate(flow_type = ifelse(p_i == "p", "Permanent", "Intermittent"))

df_summary <- df %>%
  group_by(flow_type, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE) / sqrt(n()),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    range = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )


df_summary %>%
  arrange(flow_type, metric) %>%
  gt(groupname_col = "flow_type") %>%
  tab_header(title = "Summary Statistics by Flow Type and Metric") %>%
  fmt_number(columns = where(is.numeric), decimals = 3)

# ---- Visualization boxplot ----
ggplot(df, aes(x = metric, y = value, fill = flow_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  #geom_jitter(width = 0.1, alpha = 0.8, shape = 21) +
  labs(title = "Metric Distributions by Flow Type",
       x = "Metric", y = "Value") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Csoportátlag -Sørensen---------------

df_summary <- df_long_sorr %>%
  group_by(kodex, metric, p_i) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")


df_summary$p_i <- factor(df_summary$p_i, levels = c("i", "p"), labels = c("intermittent", "permanent"))


df_summary_s_range <- df_summary %>%
  group_by(metric, p_i) %>%
  summarise(
    mean = mean(mean_value, na.rm = TRUE),
  
    
    q25 = quantile(mean_value, 0.25, na.rm = TRUE),
    
    q75 = quantile(mean_value, 0.75, na.rm = TRUE),
    min = min(mean_value, na.rm = TRUE),
    max = max(mean_value, na.rm = TRUE),
    
    n = n(),
    .groups = "drop"
  )



df_summary_s_range %>%
  arrange(p_i, metric) %>%
  gt(groupname_col = "p_i") %>%
  tab_header(title = "Summary Statistics by Flow Type and Metric") %>%
  fmt_number(columns = where(is.numeric), decimals = 3)


df_summary %>% 
  ggplot(aes(x = kodex, y = mean_value, color = p_i, group = interaction(p_i, metric))) +
  geom_line(size = 0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("#CC6600","#33CC99")) +
  labs(x = "kodex", y = "Value", 
       title = "Sørensen Metric with Group Averages") +
  theme_minimal() +
  facet_wrap(~metric) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(text = element_text(size = 14))



ggplot() +
  
  geom_line(data = df_summary, aes(x = round(as.numeric(kodex),0), y = mean_value, color = metric, linetype = p_i),
            size = 1) +
  facet_wrap(~p_i ,labeller = label_value) +
  labs(x = "kodex", y = "Value", title = "Sørensen Metric with Group Averages") +
  theme_minimal()


#soorensen violin ---------------
my_comparisons <- list( c("p", "i"))

df_long_sorr$p_i2 <- factor(df_long_sorr$p_i, levels = c("i", "p"), labels = c("intermittent", "permanent"))

df_long_sorr %>% dplyr::filter(metric %in% c("overlap" ,"gain" ,"loss"  )) %>% 
ggviolin(., x = "p_i", y = "value", fill = "p_i", 
         palette = c("#33CC99", "#CC6600"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 1.3)+
  facet_wrap(~metric)+
  ggtitle("Sørensen Components")

df_long_sorr %>% dplyr::filter(metric %in% c("gaining.turnover" ,"neutral.turnover" ,"losing.turnover"  )) %>% 
  
  ggviolin(., x = "p_i", y = "value", fill = "p_i", 
           palette = c("#33CC99", "#CC6600"),
           add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 1.3)+
  facet_wrap(~metric)+
  ggtitle("Sørensen turnover Components")

df_long_sorr %>% dplyr::filter(metric %in% c("gaining.nestedness" ,"neutral.nestedness" ,"losing.nestedness"  )) %>% 
  
  ggviolin(., x = "p_i", y = "value", fill = "p_i", 
           palette = c("#33CC99", "#CC6600"),
           add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 1.3)+
  facet_wrap(~metric)+
  ggtitle("Sørensen nestedness Components")
#jaccard violin ---------------

df_long_jacc %>% dplyr::filter(metric %in% c("overlap" ,"gain" ,"loss"  )) %>% 
  
  ggviolin(., x = "p_i", y = "value", fill = "p_i", 
           palette = c("#33CC99", "#CC6600"),
           add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 1.3)+
  facet_wrap(~metric)+
  ggtitle("Jaccard Components")

df_long_jacc %>% dplyr::filter(metric %in% c("gaining.turnover" ,"neutral.turnover" ,"losing.turnover"  )) %>% 
  
  ggviolin(., x = "p_i", y = "value", fill = "p_i", 
           palette = c("#33CC99", "#CC6600"),
           add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 1.3)+
  facet_wrap(~metric)+
  ggtitle("Jaccard turnover Components")

df_long_jacc %>% dplyr::filter(metric %in% c("gaining.nestedness" ,"neutral.nestedness" ,"losing.nestedness"  )) %>% 
  
  ggviolin(., x = "p_i", y = "value", fill = "p_i", 
           palette = c("#33CC99", "#CC6600"),
           add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 1.3)+
  facet_wrap(~metric)+
  ggtitle("Jaccard nestedness Components")




