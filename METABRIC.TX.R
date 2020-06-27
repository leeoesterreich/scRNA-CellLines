# module load gcc/8.2.0
# module load r/3.6.0
library(GSVA)
setwd('/bgfs/alee/chelsea/projects/10X/CellLine/data')
df = read.csv('METABRIC/merged.METABRIC.csv', row.names = 1, check.names = F) # merged data cleaned in ~/PublicData/METABRIC/code, merged nolan.survival.Data and other clinical data from cbioportal https://www.cbioportal.org/study/clinicalData?id=brca_metabric 
expr = t(df[,13:19496])

parse_gmt = function(gmt_path){
  db = readLines(gmt_path)
  geneset = list()
  for (line in db){
    line = 
      words = as.vector(strsplit(line, "\\s{1,}")[[1]])
    set_name = words[1]
    genes = words[-c(1,2)]
    geneset[[set_name]] = genes
  }
  return(geneset)}

dbpath = 'GSVA/ApopSig.gmt'
geneset = parse_gmt(dbpath) 

resPrefix = 'GSVA/Metabric_ApopSig'
res = gsva(expr, geneset, method='ssgsea')
write.csv(res, file=paste(resPrefix, '.ssgsea.csv', sep=''), quote=F)

# plot.dat
library("survminer")
library(dplyr)
require("survival")
library(gtools)
coxConvert = function(x){ x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                        }



# gsva 
setwd('/bgfs/alee/chelsea/projects/10X/CellLine/data/GSVA')
ssgsea = as.matrix(read.csv('Metabric_ApopSig.ssgsea.csv', check.names = F, row.names = 1))

base = df[,1:12]
ssgsea = t(ssgsea)[rownames(base),]

base$ApopSig = t(ssgsea)[1,]

table(base[, c( "NOT_IN_OSLOVAL_Pam50Subtype","Hormone Therapy","Chemotherapy" )])

# , , Chemotherapy = NO

#                            Hormone Therapy
# NOT_IN_OSLOVAL_Pam50Subtype      NO YES
#                      Basal    0  99  59
#                      Her2     0  65  81
#                      LumA     0 221 438 \\ selected individually
#                      LumB     0  91 348
#                      NC       0   2   4
#                      Normal   0  59  94

# , , Chemotherapy = YES

#                            Hormone Therapy
# NOT_IN_OSLOVAL_Pam50Subtype      NO YES
#                      Basal    0 131  37
#                      Her2     0  64  28
#                      LumA     0   3  54 
#                      LumB     0   6  41
#                      NC       0   0   0
#                      Normal   0  13  31


# ER pos and LumA, no chemo, w/ hormone therapy
## select
plot.dat = base[base[,"ER_IHC_status"]=='pos' & base[,"NOT_IN_OSLOVAL_Pam50Subtype"]=='LumA' & base[,'Chemotherapy']=='NO' &  base[,'Hormone Therapy']=='YES',] 

########## ssgsea
x.cut= surv_cutpoint(plot.dat, time="DSS.time", event="DSS.status", variables="ApopSig")
x.cat <- surv_categorize(x.cut)
plot.dat$cat = x.cat$ApopSig

cox <- coxph(Surv(DSS.time, DSS.status) ~ ApopSig, data = plot.dat)
res.cox = coxConvert(cox)
print(summary(res.cox))
print(res.cox)
            #   beta HR (95% CI for HR)          wald.test            p.value 
            # "-2.3" "0.1 (0.0093-1.1)"              "3.6"            "0.059" 
fit <- survfit(Surv(DSS.time, DSS.status) ~ ApopSig, data = x.cat)
g = ggsurvplot(fit, risk.table = TRUE, pval = TRUE, conf.int = FALSE, risk.table.height = 0.4)
g$plot = g$plot + labs(y = "Disease Free Survival") 
pdf('../../Figure/Fig3/MetabricLumA.HormTx.DSS.ApopSig_ssgsea.pdf', width=4.8, height=4)
print(g$plot)
dev.off()

# ER pos and LumA, no chemo, w/o hormone therapy
## select
plot.dat = base[base[,"ER_IHC_status"]=='pos' & base[,"NOT_IN_OSLOVAL_Pam50Subtype"]=='LumA' & base[,'Chemotherapy']=='NO' &  base[,'Hormone Therapy']=='NO',] 

########## ssgsea
x.cut= surv_cutpoint(plot.dat, time="DSS.time", event="DSS.status", variables="ApopSig")
x.cat <- surv_categorize(x.cut)
plot.dat$cat = x.cat$ApopSig

cox <- coxph(Surv(DSS.time, DSS.status) ~ ApopSig, data = plot.dat)
res.cox = coxConvert(cox)
print(res.cox)
            #   beta HR (95% CI for HR)          wald.test            p.value 
            # "0.23"    "1.3 (0.02-78)"             "0.01"             "0.91" 
fit <- survfit(Surv(DSS.time, DSS.status) ~ ApopSig, data = x.cat)
g = ggsurvplot(fit, risk.table = TRUE, pval = TRUE, conf.int = FALSE, risk.table.height = 0.4)
g$plot = g$plot + labs(y = "Disease Free Survival") 
pdf('../../Figure/Fig3/MetabricLumA.NoHormTx.DSS.ApopSig_ssgsea.pdf', width=4.8, height=4)
print(g$plot)
dev.off()