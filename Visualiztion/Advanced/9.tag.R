# ggpp--------------------------------------------------------------------------
# https://cloud.tencent.com/developer/article/2413512
corner_letters.tb <- tibble(label = LETTERS[1:4],
                            x = "left", 
                            y = "top",
                            cyl = c(4,5,6,8))

ggplot(mpg, aes(displ,hwy)) +
  geom_point() +
  facet_wrap(~cyl, scales = "free") +
  ggpp::geom_text_npc(data = corner_letters.tb,
                      aes(npcx = x, npcy = y, label = label)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
# patchwork---------------------------------------------------------------------