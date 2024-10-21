require(UpSetR)

# UpSetR is a method like venn plot with more groups
# Note: the zero doesn't mean there are no always intersection within these groups
#       but there are no exclusive intersection within these groups

readxl::read_excel("8-38347143-Source Data Extended Fig8.xlsx", sheet = "ED_8a") %>% 
  transmute(PatientCode = PatientCode,
            `C/I (N=773)` = `Clinical stage/Subtype`,
            `T (N=752)` = `Transcriptomic`,
            `P (N=626)` = `Pathologic`,
            `M (N=453)` = `Metabolomic`,
            `R (N=419)` = `Radiologic`) %>% 
  as.data.frame() %>% 
  upset(
    nsets = 5,
    sets = c("C/I (N=773)",	"T (N=752)",	"P (N=626)", "M (N=453)",	"R (N=419)"), 
    keep.order = TRUE,
    point.size = 3,
    line.size = 0,
    order.by = "freq",
    empty.intersections = TRUE,
    # show.numbers = FALSE,
    mainbar.y.label = "Modality Combination Size",
    cutoff = 0
  )