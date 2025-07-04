#usage: Rscript BayesPrism_gen.r --ID=<job_id>
#prerequesites: you have run good_main.py 

load("myinput_gen.gbm.rdata")

suppressWarnings(library(BayesPrism))
args = commandArgs(trailingOnly = TRUE)
ID <- "default_id"
if(length(args) > 0){
  for (arg in args){
    if (grepl("--ID=", arg)) {
        print(args)
            ID <- sub("--ID=", "", arg)
        }
  }
}


sc.dat.filtered <- cleanup.genes (input=sc.dat, input.type="count.matrix", species="hs", gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") , exp.cells=5)
myPrism <-new.prism(
  reference=sc.dat,#not filtered atm because the discreptancy is too large
  mixture=bk.dat,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key=NULL,
  outlier.cut=0.001,#changed these from default values
  outlier.fraction=0.05,#changed these two from default values
)
bp.res <- run.prism(prism = myPrism, n.cores=1)
theta <- get.fraction (bp=bp.res,
            which.theta="final",
            state.or.type="type")
#write.table(theta, file = "theta.csv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
#theta <- read.csv("theta.csv", sep = "\t", check.names = FALSE, row.names = 1)
#theta$Mixture <- rownames(theta)

# Now set the correct column names
colnames(theta) <- c(
  #"Mixture",
  "LFQ.intensity.imputed_B.memory_04_steady-state",
  "LFQ.intensity.imputed_B.naive_04_steady-state",
  "LFQ.intensity.imputed_B.plasma_04_steady-state",
  "LFQ.intensity.imputed_mTregs_04_steady-state",
  "LFQ.intensity.imputed_T4.naive_04_steady-state",
  "LFQ.intensity.imputed_nTregs_04_steady-state",
  "LFQ.intensity.imputed_T4.CM_04_steady-state",
  "LFQ.intensity.imputed_T4.EM_04_steady-state",
  "LFQ.intensity.imputed_T4.EMRA_04_steady-state",
  "LFQ.intensity.imputed_Th1_04_steady-state",
  "LFQ.intensity.imputed_Th17_04_steady-state",
  "LFQ.intensity.imputed_Th2_04_steady-state",
  "LFQ.intensity.imputed_T8.naive_04_steady-state",
  "LFQ.intensity.imputed_T8.CM_04_steady-state",
  "LFQ.intensity.imputed_T8.EM_04_steady-state",
  "LFQ.intensity.imputed_T8.EMRA_04_steady-state",
  "LFQ.intensity.imputed_mDC_04_steady-state",
  "LFQ.intensity.imputed_pDC_04_steady-state",
  "LFQ.intensity.imputed_Basophil_04_steady-state",
  "LFQ.intensity.imputed_Eosinophil_04_steady-state",
  "LFQ.intensity.imputed_Neutrophil_04_steady-state",
  "LFQ.intensity.imputed_MO.classical_04_steady-state",
  "LFQ.intensity.imputed_MO.intermediate_04_steady-state",
  "LFQ.intensity.imputed_MO.nonclassical_04_steady-state",
  "LFQ.intensity.imputed_NK.bright_04_steady-state",
  "LFQ.intensity.imputed_NK.dim_04_steady-state"
)
output_filename <- paste0("BAYESPRISM_gen-Results",ID,".txt")
updated_sig <- bp.res@reference.update@phi
output_update_filename <- paste0("BAYESPRISM-updated_sig_matrix",ID,".txt")
write.table(updated_sig, file = output_update_filename, sep ="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(theta, file = output_filename, sep ="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
save.image(file = "myinput3.gbm.rdata")