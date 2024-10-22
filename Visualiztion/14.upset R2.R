# upset manually

require(ggplot2)

dat <- readxl::read_excel("8-38347143-Source Data Extended Fig8.xlsx", sheet = "ED_8a") %>% 
  transmute(PatientCode = PatientCode,
            `C/I (N=773)` = `Clinical stage/Subtype`,
            `T (N=752)` = `Transcriptomic`,
            `P (N=626)` = `Pathologic`,
            `M (N=453)` = `Metabolomic`,
            `R (N=419)` = `Radiologic`) %>% 
  as.data.frame()
# ------------------------------------------------------------------------------
dat.p1 <- dat %>% 
  tidyr::pivot_longer(cols = `C/I (N=773)`:`R (N=419)`, 
                      names_to = "variable", 
                      values_to = "Freq") %>% 
  filter(Freq == 1) %>% 
  group_by(variable) %>% 
  mutate(sum = sum(Freq)) %>% 
  arrange(desc(sum)) 
ord <- dat.p1 %>% distinct(sum) %>% dplyr::select(variable) %>% unlist
p1 <- dat.p1 %>% 
  mutate(variable = factor(variable, levels = rev(ord))) %>% 
  ggplot(aes(x = variable)) +
    geom_hline(yintercept = seq(0,800,100), color = "gray80") +
    geom_bar() +
    labs(y = "Set size", x = "") +
    scale_x_discrete(breaks = NULL) +
    scale_y_reverse(breaks = NULL, sec.axis = dup_axis(name = "", breaks = c(0,400,800))) +
    coord_flip() +
    theme_minimal() +
    theme(text = element_text(size = 16),
          axis.title.x = element_text(vjust = -1))
p1
# ------------------------------------------------------------------------------
# dat <- dat %>% 
#   mutate(PatientCode = factor(PatientCode, levels = ord))
num.row <- sum(choose(5,c(1:5)))
dat.p2 <- data.frame(
  sets = rep(NA, num.row),
  n = rep(NA, num.row),
  stringsAsFactors = F
)
n <- i <- 1
while(n <= length(ord)){
  tmp <- combn(ord, n) %>% as.matrix %>% as.data.frame
  for (j in 1:ncol(tmp)) {
    dat.p2$sets[i] <- paste0(tmp[[j]], collapse = "-")
    
    dat.p2$n[i] <- dat %>% 
      dplyr::select(tmp[[j]]) %>% 
      mutate(sum = rowSums(.)) %>% 
      filter(sum == n) %>% 
      nrow

    i <- i + 1 
  }
  n <- n + 1 
}

dat.p2[nrow(dat.p2) + 1, ] <- c("", 0)

dat.p2 <- dat.p2 %>% 
  arrange(desc(n)) %>% 
  mutate(n = as.numeric(n),
         row = row_number())

p2 <-
  ggplot(dat.p2, aes(x = row, y = n)) +
    geom_hline(yintercept = seq(0,800,100), color = "gray80") +
    geom_vline(aes(xintercept = seq(1,nrow(dat.p2),1)), color = "gray80") +
    geom_col() +
    scale_x_continuous(expand = c(0,0), breaks = NULL) +
    scale_y_continuous(breaks = seq(0,800,200)) +
    labs(x = "",  y = "Modality combination size") +
    # coord_cartesian(clip = "off") +
    theme_minimal() +
    theme(text = element_text(size = 16),
          axis.line.x = element_blank(),
          plot.margin = ggplot2::margin(0,30,0,10))
p2
# ------------------------------------------------------------------------------
dat.p3 <- data.frame(
  row = rep(dat.p2$row,length(ord)),
  x = rep(dat.p2$row,length(ord)),
  y = rep(ord,nrow(dat.p2)),
  stringsAsFactors = FALSE
) %>% 
  left_join(dat.p2, by = "row") %>% 
  mutate(y = factor(y, levels = rev(ord)),
         string = str_remove_all(sets, " \\("),
         string = str_remove_all(string, "\\)"),
         pattern = str_remove_all(y, " \\("),
         pattern = str_remove_all(pattern, "\\)"),
         col = case_when(str_detect(string, pattern) ~ 1,
                         TRUE ~ 0)) %>% 
  mutate(col = factor(col))
p3 <- 
  ggplot(dat.p3, aes(x = x, y = y, fill = col)) +
  geom_vline(aes(xintercept = x), color = "gray80") +
  geom_point(size = 4, pch = 21, color = "gray80") +
  scale_fill_manual(values = c("gray50","black")) +
  scale_x_continuous(limits = c(1,nrow(dat.p2)), expand = c(0.015,0), breaks = NULL) +
  labs(x = "Modality combination", y = "") +
  guides(fill = "none") +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        plot.margin = ggplot2::margin(0,30,30,10),
        axis.text.y = element_text(hjust = 0),
        axis.title.x = element_text(vjust = -1))
p3
# ------------------------------------------------------------------------------
require(patchwork)

design <- "
  #2222
  #2222
  #2222
  13333
  13333
"

p_upset <- p1 + p2 + p3 + plot_layout(design = design)



