# module load gcc/8.2.0
# module load r/3.6.0
library(GSVA)
setwd('/bgfs/alee/chelsea/projects/10X/CellLine/data')
df = read.csv('METABRIC/nolan.surv.expr.csv', row.names = 1, check.names = F)
expr = t(df[,8:19496])

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

# ER pos and LumA
## select
plot.dat = base[base[,"ER_IHC_status"]=='pos' & base[,"NOT_IN_OSLOVAL_Pam50Subtype"]=='LumA',] 

########## ssgsea
x.cut= surv_cutpoint(plot.dat, time="DSS.time", event="DSS.status", variables="ApopSig")
x.cat <- surv_categorize(x.cut)
plot.dat$cat = x.cat$ApopSig

cox <- coxph(Surv(DSS.time, DSS.status) ~ ApopSig, data = plot.dat)
print(summary(cox))
# Call:
# coxph(formula = Surv(DSS.time, DSS.status) ~ ApopSig, data = plot.dat)

#   n= 695, number of events= 141 

#            coef exp(coef) se(coef)      z Pr(>|z|)  
# ApopSig -1.7514    0.1735   1.0219 -1.714   0.0865 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#         exp(coef) exp(-coef) lower .95 upper .95
# ApopSig    0.1735      5.763   0.02342     1.286

# Concordance= 0.555  (se = 0.026 )
# Likelihood ratio test= 2.9  on 1 df,   p=0.09
# Wald test            = 2.94  on 1 df,   p=0.09
# Score (logrank) test = 2.94  on 1 df,   p=0.09
res.cox = coxConvert(cox)
print(res.cox)
            #   beta HR (95% CI for HR)          wald.test            p.value 
            # "-1.8" "0.17 (0.023-1.3)"              "2.9"            "0.087" 
fit <- survfit(Surv(DSS.time, DSS.status) ~ ApopSig, data = x.cat)
g = ggsurvplot(fit, risk.table = TRUE, pval = TRUE, conf.int = FALSE, risk.table.height = 0.4)
g$plot = g$plot + labs(y = "Disease Free Survival") 
pdf('../../Figure/Fig3/MetabricLumA.DSS.ApopSig_ssgsea.pdf', width=4.8, height=4)
print(g)
dev.off()

