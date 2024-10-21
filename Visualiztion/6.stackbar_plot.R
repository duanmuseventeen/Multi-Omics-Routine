# TGF-β attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells
# http://research-pub.gene.com/IMvigor210CoreBiologies/
# PMID: 29443960
# PRJNA482620, Braun_2020, IMvigor210, GSE93157
# 
# RNA sequencing
# All raw sequencing data required for RNAseq analyses have been deposited to 
# the European Genome-Phenome Archive under accession number EGAS#00001002556. 
# In addition, the source code and processed data used for all analyses presented 
# here have been made available in IMvigor210CoreBiologies, a fully documented 
# software and data package for the R statistical computing environment (R Core 
# Team, 2016. R: A language and environment for statistical computing. R 
# Foundation for Statistical Computing, Vienna, Austria. 
# https://www.r-project.org). This package is freely available under the Creative 
# Commons 3.0 license and can be downloaded from 

require(dplyr)
require(GEOquery)

load("IMvigor210CoreBiologies.Rdata")

pdata <- pdata %>%
  filter(complete.cases(binaryResponse)) %>%
  mutate(Response = ifelse(binaryResponse == "CR/PR","R","NR"),
         EVENT = censOS,
         OS = os)

dim(pdata)

IMvigor <- expMatrix_tpm %>% as.data.frame %>%
  mutate(mean = rowMeans(.)) %>%
  mutate(entrez_id = rownames(.)) %>%
  left_join(feature , by = "entrez_id") %>%
  arrange(desc(mean)) %>%
  filter(!duplicated(symbol)) %>%
  filter(complete.cases(symbol))
rownames(IMvigor) <- IMvigor$Symbol
IMvigor_response <- IMvigor %>% 
  dplyr::select(-c(mean, entrez_id, symbol, n_exons, length, source, Symbol)) %>%
  t %>% as.data.frame %>% 
  mutate(ID = rownames(.)) %>%
  left_join(pdata %>% mutate(ID = rownames(pdata)), by = "ID") %>%
  mutate(Race = ifelse(Race == "UNKNOWN", NA, Race)) %>%
  filter(complete.cases(`IC Level`)) %>%
  mutate(`Best Confirmed Overall Response` = factor(`Best Confirmed Overall Response`,
                                                    levels = c("PD","SD","PR","CR")))

dat02 <- data.frame(
  label = c("(83)","(112)","(102)"),
  `IC Level` = factor(c(1:3), labels = c("IC0","IC1","IC2+")),
  stringsAsFactors = F,
  check.names = F
)

# 1
p <- ggplot() +
  geom_bar(data = IMvigor_response, aes(x = `IC Level`, fill = `Best Confirmed Overall Response`),
           position = 'fill', color = "black") + 
  geom_text(data=dat02,
            aes(x=`IC Level`,y=1, label = label, vjust=-0.5, size = 2
            )) +
  scale_fill_manual(values = c("#006eb3","#86c6da","#f7a27c","#ca1f25")) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), expand = c(0,0)) + 
  xlab("") + 
  ylab("Fraction of patients") +
  ggtitle("PD-L1 IC") +
  guides(size = "none") + 
  coord_cartesian(clip = "off") + 
  theme_classic() +
  theme(text = element_text(size = 20), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 5),
        axis.line.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.background = element_blank(),
        plot.margin = ggplot2::margin(50,30,10,30), 
        panel.grid = element_blank(), legend.background = element_blank(),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4)
  )

ggsave(plot = p, filename = "29443960_fig1a.pdf", path = "D:/", device = "pdf",width = 4, height = 6, dpi = 300)

# 2
dat <- 
  table(IMvigor_response$`IC Level`,
        IMvigor_response$`Best Confirmed Overall Response`) %>%
  as.data.frame() %>%
  group_by(Var1) %>%
  mutate(sum = sum(Freq)) %>%
  mutate(percent = Freq/sum) %>%
  mutate(label = sprintf("%.2f", percent*100) %>% paste0("%"))

p2 <- ggplot(dat, aes(x = Var1, y = percent, fill = Var2, label = label)) +
  geom_bar(stat="identity",position = position_stack(), color = "black") +
  geom_text(stat="identity",position = position_stack(), aes(label=label),vjust = -0.1,size=4,color="black")+
  scale_fill_manual(values = c("#006eb3","#86c6da","#f7a27c","#ca1f25")) +
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), expand = c(0,0)) +
  xlab("") + 
  ylab("Fraction of patients") +
  ggtitle("PD-L1 IC") +
  guides(size = "none") + 
  coord_cartesian(clip = "off") + 
  theme_classic() +
  theme(text = element_text(size = 20), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 5),
        axis.line.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.background = element_blank(),
        plot.margin = ggplot2::margin(50,30,10,30), 
        panel.grid = element_blank(), legend.background = element_blank(),
        axis.title.x = element_text(vjust = -4),
        axis.title.y = element_text(vjust = 4)
  )
p2
ggsave(plot = p2, filename = "29443960_fig1a R2.pdf", path = "D:/", device = "pdf",width = 4.5, height = 6, dpi = 300)
