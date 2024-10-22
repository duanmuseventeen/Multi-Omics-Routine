# heatmap for confusion matrix

require(tidymodels)
require(tidyverse)
require(themis)

volcano_raw <- data.table::fread("13.volcano.csv")

volcano_raw %>%
  count(primary_volcano_type, sort = TRUE)

volcano_df <-
  volcano_raw %>%
  transmute(volcano_type = case_when(
    str_detect(primary_volcano_type, "Stratovolcano") ~ "Stratovolcano",
    str_detect(primary_volcano_type, "Shield") ~ "Shield",
    TRUE ~ "Other"),
    volcano_number, latitude, longitude, elevation,
    tectonic_settings, major_rock_1) %>%
  mutate_if(is.character, factor)

world <- map_data("world")

ggplot() +
  geom_map(data = world, map = world,
           aes(long, lat, map_id = region),
           color = "white", fill = "gray50", alpha = 0.2) +
  geom_point(data = volcano_df,
             aes(longitude, latitude, color = volcano_type),
             alpha = 0.8)

volcano_boot <- bootstraps(volcano_df)
volcano_boot

volcano_rec <- recipe(volcano_type ~ ., data = volcano_df) %>%
  update_role(volcano_number, new_role = "Id") %>%
  step_other(tectonic_settings) %>%
  step_other(major_rock_1) %>%
  step_dummy(tectonic_settings, major_rock_1) %>%
  step_zv(all_predictors()) %>%
  step_smote(volcano_type)

volcano_prep <- prep(volcano_rec)

rf_spec <- rand_forest(trees = 1000) %>%
  set_mode("classification") %>%
  set_engine("ranger")

volcano_wf <- workflow() %>%
  add_recipe(volcano_rec) %>%
  add_model(rf_spec)

volcano_wf

volcano_res <- fit_resamples(
  volcano_wf,
  resamples = volcano_boot,
  control = control_resamples(save_pred =  TRUE,
                              verbose = TRUE)
)

volcano_res %>% 
  collect_metrics()

tmp <- volcano_res %>% 
  collect_predictions() %>% 
  conf_mat(volcano_type, .pred_class)

heat <- tmp$table %>% 
  as_tibble() %>% 
  group_by(Truth) %>% 
  mutate(sum = sum(n),
         accuracy = n/sum) %>% 
  mutate(x = Truth %>% factor) %>% 
  mutate(x = x %>% as.numeric) %>% 
  mutate(y = Prediction %>% factor %>% as.numeric) %>%
  mutate(color = ifelse(x == y, 1, 0) %>% factor)

heat.highlight <- heat %>% 
  filter(color == 1)
ggplot() +
  geom_tile(data = heat, aes(x, y, fill = accuracy)) +
  geom_tile(data = heat.highlight, 
            aes(x, y), fill = NA,linewidth = 2, color = "black") +
  geom_text(data = heat, aes(x, y, label = n)) +
  scale_x_continuous(breaks = c(1:3), expand = c(0,0),
                     labels = c("Other (3438)",
                                "Shield (1082)",
                                "Stratovolcano (4245)"),
                     sec.axis = dup_axis(
                       labels = c("56.98% (1959/3438)",
                                  "48.71% (527/1082)",
                                  "76.14% (3232/4245)")
                     )) +
  scale_y_continuous(breaks = c(1:3), expand = c(0,0), 
                     labels = c("Other (3113)",
                                "Shield (944)",
                                "Stratovolcano (4708)"),
                     sec.axis = dup_axis(
                       labels = c("62.93% (1959/3113)",
                                  "55.83% (527/944)",
                                  "68.65% (3232/4708)")
                     )) +
  scale_fill_gradient2(low = "gray90", high = scales::muted("red")) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1),
        axis.text.x.top = element_text(angle = 45, hjust = 0),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        plot.margin = ggplot2::margin(40,40,40,40))

# vip plot----------------------------------------------------------------------
library(vip)

rf_spec %>% 
  set_engine("ranger", importance = "permutation") %>% 
  fit(
    volcano_type ~ .,
    data = juice(volcano_prep) %>% 
      select(-volcano_number) %>% 
      janitor::clean_names()
  ) %>% 
  vip(geom = "point")

# world map + stat--------------------------------------------------------------
volcano_pred <- volcano_res %>% 
  collect_predictions() %>% 
  mutate(correct = volcano_type == .pred_class) %>% 
  left_join(volcano_df %>% 
              mutate(.row = row_number()))

ggplot() +
  geom_map(data = world, map = world, 
           aes(long, lat, map_id = region),
           color = "white", fill = "gray50", alpha = 0.2) +
  stat_summary_hex(data = volcano_pred,
                   aes(longitude, latitude, z = as.integer(correct)),
                   fun = "mean",
                   alpha = 0.7,
                   bins = 60)
