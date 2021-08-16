# SRM - Peptideos ----

# Dados ----

# Tecidos:
# np: nasal polyps 
# p: parotid gland
# t: palatine tonsils
# pool: sample pool

# Metodos de preparacao :
#in-gel (IGD)
#in-solution (ISD)
#on-filter (OFD)
#on-pellet (OPD)


# Pacotes
library(ggplot2)
library(nlme)
library(visreg)
library(emmeans)
library(predictmeans)
library(writexl)


library(pcalg)
library(gRbase)
library(Rgraphviz)
library(RBGL)
library(robustbase)
library(Rcpp)
library(igraph)

# Carregar os dados em formato .RData

ggplot(d,aes(sort(Fn),log10(y)))+ geom_boxplot()+
  theme(legend.position="none")+
  labs(title="SRM - Quantificação de Peptídeo - log10(y)", x = "Peptídeo")

ggplot(d,aes(gr,log10(y),color=Tecido))+ geom_boxplot()+
  labs(title="SRM - Quantificação de Peptídeo - log10(y)", x = "Método(Tecido)")


ggplot(d,aes(gr,y,color=Tecido))+ geom_boxplot()+
  labs(title="SRM - Quantificação de Peptídeo - y", x = "Método(Tecido)")

ggplot(d,aes(x = Metodo, color = Tecido, group =Tecido, y = y)) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line") +
  labs(title="SRM - Perfis y", x = "Método")

ggplot(d,aes(x = Metodo, color = Tecido, group = Tecido, y = log10(y))) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line") +
  labs(title="SRM - Perfis log10(y)", x = "Método")

# Modelos de normalização ----

# Modelo Linear Misto Homocedastico ----
#y = (mu + u1 + u2) + M + e 



fit1 <- lme(fixed = log10(y) ~ M, na.action = na.exclude, random = ~ 1|u1/u2, data = d)
fit1
anova(fit1)
summary(fit1)

intervals(fit1)

#intervalos
round(intervals(fit1)[1]$fixed[,2]-intervals(fit1)[1]$fixed[,1],4)
round(intervals(fit1)[2]$fixed[,2]-intervals(fit1)[1]$fixed[,1],4)

round(intervals(fit1)[2]$reStruct$u1[,2] - intervals(fit1)[2]$reStruct$u1[,1],4)
round(intervals(fit1)[2]$reStruct$u2[,2] - intervals(fit1)[2]$reStruct$u2[,1],4)

round(intervals(fit1)[3]$sigma[2] - intervals(fit1)[3]$sigma[1],4)


# Estimativas das medias marginais na escala transformada
emm1.1 <- emmeans(fit1, "M", sigmaAdjust = TRUE, mode = "containment" ,data = d)
emm1.1

# Estimativas das medias marginais na escala original
emm1.1 <- emmeans(fit1, "M", sigmaAdjust = TRUE, transform = "response", mode = "containment" ,data = d)

# Contrastes
pairs(emm1.1)

# Avaliação do ajuste do modelo
residplot(fit1, level = 2)


# Modelo Linear Misto Heteroscedastico devido ao Tecido ----
#y = (mu + u1 + u2) + M + e Var(e) = 1|T


fT <- d$T
fit2 <- lme(log10(y) ~ M, random = ~1|u1/u2, na.action = na.exclude, weights = varIdent(form = ~1|fT), control = lmeControl(maxIter = 100, opt = "optim"), data = d)
fit2
anova(fit2)

#intervalos
round(intervals(fit2, which = "fixed")[1]$fixed[,2]-intervals(fit2, which = "fixed")[1]$fixed[,1],4)
round(intervals(fit2, which = "fixed")[1]$fixed[,3]-intervals(fit2, which = "fixed")[1]$fixed[,2],4)

# Extraindo os pesos utilizados na função de variância e calculando as componentes de variância estimadas
VarCorr(fit2)
summary(fit2)$sigma #sigma(NP) é a referencia
coef(fit2$modelStruct$varStruct, unconstrained=FALSE)

# Estimativas das componentes de variancia
(summary(fit2)$sigma*coef(fit2$modelStruct$varStruct, uncons=FALSE))^2


# Estimativas das medias marginais na escala transformada
emm2.1 <- emmeans(fit2, "M", sigmaAdjust=TRUE, mode="containment",data = d)
emm2.1

# Estimativas das medias marginais na escala original
emm2.1<-emmeans(fit2, "M", sigmaAdjust=TRUE,transform = "response", mode="containment",data = d)

# Contrastes
pairs(emm2.1)

# Avaliacao do modelo ajustado
residplot(fit2, level = 2)


# Modelo Linear Misto Heteroscedastico devido ao M(T) ----
#y = (mu + u1 + u2) + M + e Var(e) = 1|M(T)


fTM <- as.factor(d$u2)
fit3 <- lme(log10(y) ~ M, random = ~1|u1/u2, na.action = na.exclude, weights = varIdent(form = ~1|fTM), control = lmeControl(maxIter = 100, opt = "optim"), data = d)
fit3
anova(fit3)

#intervalos
round(intervals(fit3, which = "fixed")[1]$fixed[,2]-intervals(fit3, which = "fixed")[1]$fixed[,1],4)
round(intervals(fit3, which = "fixed")[1]$fixed[,3]-intervals(fit3, which = "fixed")[1]$fixed[,2],4)


# Extraindo os pesos utilizados na função de variância e calculando as componentes de variância estimadas
VarCorr(fit3)
summary(fit3)$sigma ## sigma(NP:IGD) é a referencia
coef(fit3$modelStruct$varStruct, unconstrained=FALSE)
summary(fit3)$sigma*coef(fit3$modelStruct$varStruct, uncons=FALSE)

# Estimativas das componentes de variancia
(summary(fit3)$sigma*coef(fit3$modelStruct$varStruct, uncons=FALSE))^2

# Estimativas das medias marginais na escala transformada
emm3.1 <- emmeans(fit3, "M", sigmaAdjust = TRUE, mode = "containment", data = d)
emm3.1

# Estimativas das medias marginais na escala original
emm3.1 <- emmeans(fit3, "M", sigmaAdjust = TRUE, transform = "response", mode = "containment", data = d)

# Contrastes
pairs(emm3.1)

# Avaliação do modelo ajustado
residplot(fit3, level = 2)

# Testes de Razao de Verossimilhança
anova(fit1,fit2)
anova(fit2,fit3)

# Estimacao dos modelos por ML
fit1.ml <- lme(log10(y) ~ M, random = ~1|u1/u2, na.action = na.exclude, method = "ML", data = d)
fit1.ml

fT <- as.factor(d$T)
fit2.ml <- lme(log10(y) ~ M, random = ~1|u1/u2, na.action = na.exclude, weights = varIdent(form = ~1|fT), control = lmeControl(maxIter = 100, opt = "optim"), method = "ML", data = d)
fit2.ml

fTM <- as.factor(d$u2)
fit3.ml <- lme(log10(y) ~ M, random = ~1|u1/u2, na.action = na.exclude, weights = varIdent(form = ~1|fTM), control = lmeControl(maxIter = 100, opt = "optim"), method = "ML", data = d)
fit3.ml

anova(fit1.ml,fit2.ml)
anova(fit2.ml,fit3.ml)


# Modelo de Significancia ----

drms <- cbind(d,data.frame(residuals(fit3)),data.frame(fitted(fit3)))
names(drms)
str(drms)
names(drms) <- c("gr","u1","u2","u3","M","T","Metodo","Tecido","S","Fn","F","y","y_norm","pred")

ggplot(drms,aes(sort(Fn),y_norm))+ geom_boxplot()+
  theme(legend.position="none")+
  labs(title="SRM - Quantificação Peptídeos - log10(y)", x = "Peptídeo")

ggplot(drms,aes(gr,y_norm,color=Tecido))+ geom_boxplot()+
  labs(title="SRM - Quantificação Peptídeos - normalizado", x = "Método(Tecido)")


# Renomenando y_norm para res
names(drms) <- c("gr","u1","u2","u3","Mf","Tf","M","T","S","Fn","F","y","res","pred")
str(drms)

pvalues2by2 <- list()
anova_pvalues <- c()
for (J in levels(drms$F)) {
  drmsJ <- drms[which(drms$F==J),]
  drmsJ$gr <- as.factor(drmsJ$gr)
  drmsJ <- drmsJ[complete.cases(drmsJ),] #retira os Na's
  fitresJ <- tryCatch({
    gls(res ~ T/M, weights = varIdent(form = ~1|gr), #considera erros heterocedasticos dados por FTM
        control = glsControl(maxIter = 100, opt = "optim"), na.action = na.exclude, data = drmsJ)
  }, error = function(e) {
    NULL
  })
  if (!is.null(fitresJ)) {
    aovJ <- anova(fitresJ)
    anova_pvalues <- rbind(anova_pvalues, c(J, aovJ$"p-value"))
    refJ <- ref_grid(fitresJ, mode = "df.error", data=drmsJ, nesting = "M %in% T")
    methJ <- emmeans(refJ, "M")
    meth2J <- pairs(methJ, by = "T", reverse = TRUE)
    summJ <- summary(meth2J)
    contrasts <- paste0(summJ$T, ": ", summJ$contrast)
    for (contrastId in 1:length(contrasts)) {
      pvalues2by2[[contrasts[contrastId]]] <-
        rbind(pvalues2by2[[contrasts[contrastId]]], c(J,summJ$estimate[contrastId], summJ$p.value[contrastId]))
      colnames(pvalues2by2[[contrasts[contrastId]]]) <- c("Peptideo","Estimativa", "p-valor")
    }
  }
}
colnames(anova_pvalues) <- c("Peptídeo", "Intercepto", "T", "T:M")
dim(anova_pvalues)

#Grafico de vulcão ----

plot(rep(1:58),-log10(as.numeric(anova_pvalues[,4])))
min(-log10(as.numeric(anova_pvalues[,4])))

anova_pvaluesMT<-as.numeric(anova_pvalues[,4])
p.adjust.MT <- p.adjust(anova_pvaluesMT, method="fdr")

plot(anova_pvaluesMT,p.adjust.MT)
plot(rep(1:58),-log10(p.adjust.MT))
    
names(pvalues2by2)
p.bf <- -log10(0.05/((length(anova_pvalues)/4)*6)) #ponto de corte

pep_signif <- NULL
for(c in 1:24){
  prots <- pvalues2by2[[c]][-log10(as.numeric(pvalues2by2[[c]][,"p-valor"])) > p.bf & abs(as.numeric(pvalues2by2[[c]][,"Estimativa"])) > 2,]
  
  if(length(prots) == 3){
    pep_signif <- rbind(pep_signif,c(names(pvalues2by2)[[c]],
                                     pvalues2by2[[c]][-log10(as.numeric(pvalues2by2[[c]][,"p-valor"])) > p.bf & abs(as.numeric(pvalues2by2[[c]][,"Estimativa"])) > 2,]))
  }
  else{pep_signif <- rbind(pep_signif,cbind(names(pvalues2by2)[[c]],
                                            pvalues2by2[[c]][-log10(as.numeric(pvalues2by2[[c]][,"p-valor"])) > p.bf & abs(as.numeric(pvalues2by2[[c]][,"Estimativa"])) > 2,]))
  }
}


pdf("SRM Peptideos - Significance.pdf", width = 12, height = 8)

par(mfrow=c(4,6))
for (L in 1:6){ 
  log10.pvalue <- log10(as.numeric(pvalues2by2[[L]][,"p-valor"]))
  estimate <- as.numeric(pvalues2by2[[L]][,"Estimativa"])
  plot(estimate,-log10.pvalue,main=names(pvalues2by2)[L],col=ifelse(-log10.pvalue>p.bf & abs(estimate)>2,"red","black"))
  abline(h = p.bf)
  abline(v = c(-2,2))
}

for (L in 7:12){
  log10.pvalue <- log10(as.numeric(pvalues2by2[[L]][,"p-valor"]))
  estimate <- as.numeric(pvalues2by2[[L]][,"Estimativa"])
  plot(estimate,-log10.pvalue,main=names(pvalues2by2)[L],col=ifelse(-log10.pvalue>p.bf & abs(estimate)>2,"red","black"))
  abline(h = p.bf)
  abline(v = c(-2,2))
}

for (L in 13:18){
  log10.pvalue <- log10(as.numeric(pvalues2by2[[L]][,"p-valor"]))
  estimate<-as.numeric(pvalues2by2[[L]][,"Estimativa"])
  plot(estimate,-log10.pvalue,main=names(pvalues2by2)[L],col=ifelse(-log10.pvalue>p.bf & abs(estimate)>2,"red","black"))
  abline(h = p.bf)
  abline(v = c(-2,2))
}

for (L in 19:24){
  log10.pvalue <- log10(as.numeric(pvalues2by2[[L]][,"p-valor"]))
  estimate <- as.numeric(pvalues2by2[[L]][,"Estimativa"])
  plot(estimate,-log10.pvalue,main=names(pvalues2by2)[L],col=ifelse(-log10.pvalue>p.bf & abs(estimate)>2,"red","black"))
  abline(h = p.bf)
  abline(v = c(-2,2))
}
par(mfrow=c(1,1))
dev.off()

# Base para analise multivariada 
# Dados completos

df <- as.data.frame(table(d$F))
variaveis <- df$Var1[df$Freq == 80]  
dcomp <- d[d$F %in% variaveis,]
  
names(dcomp)<-c("grc","u1c","u2c","u3c","Mc","Tc", "Mfc","Tfc","Sc","Fnc","Fc", "yc")
str(dcomp)


fitdc <- lme(log10(yc) ~ Mc, random = ~1|u1c/u2c, weights = varIdent(form = ~1|u2c)
             ,control = lmeControl(maxIter = 100, opt = "optim"), na.action = na.exclude, data = dcomp)

names(dcomp)
drcomp<-cbind(dcomp,data.frame(residuals(fitdc)),data.frame(fitted(fitdc)))
names(drcomp)
dim(drcomp)
str(drcomp)
names(drcomp) <- c("grc","u1c","u2c","u3c","Mc","Tc","Mfc","Tfc", "Sc","Fnc","Fc", "yc","resc","predc")

sdresc <- c()
for (J in variaveis) {
  drcJ <- drcomp[which(drcomp$Fc==J),]
  dru2 <- as.factor(drcJ$u2c)
  fitrescJ <- gls(resc ~ Tc/Mc, weights = varIdent(form = ~1|dru2),
                control = glsControl(maxIter = 100, opt = "optim"), na.action = na.exclude, data = drcJ)
  sdrcJ <- summary(fitrescJ)$sigma*c(1,coef(fitrescJ$modelStruct$varStruct, uncons=FALSE)) #estimativas do desvio padrão de TxM para cada F
  sdresc <- rbind(sdresc, c(J,sdrcJ))
}

dim(sdresc)
sdresc<-data.frame(sdresc)
names(sdresc)
str(sdresc)
names(sdresc) <- c("F","NP:IGD","NP:ISD","NP:OFD","NP:OPD",
                 "P:IGD","P:ISD","P:OFD","P:OPD",
                 "T:IGD","T:ISD","T:OFD","T:OPD",
                 "Pool:IGD","Pool:ISD","Pool:OFD","Pool:OPD")

sdresc$F <- as.factor(sdresc$F)

# Selecionando as colunas de interesse F, Sc, Tc, Mc, resc
ynorm<-drcomp[,c(11,9,8,7,13)] 
str(ynorm)

# Pesos dados pelos sigmas obtidos para TxM para cada F
weights <- sdresc

ynorm2 <- reshape(ynorm, direction = "wide", timevar = "Fc", idvar = c("Tfc","Mfc","Sc"))
ynorm.homo <- ynorm2
tec <- c("NP", "P", "T", "Pool")
met <- c("ISD", "IGD", "OPD", "OFD")

#padronizando os residuos para cada F
for (teci in tec) {
  for (meti in met) {
    ids <- which(ynorm2$Tfc == teci & ynorm2$Mfc ==  meti)
    for (var in 1:46) {
      wei <- as.numeric(as.character(weights[var, paste0(teci, ":", meti)]))
      ynorm.homo[ids,paste0("resc.",weights$F[var])] <- ynorm.homo[ids, paste0("resc.",
                                                                     weights$F[var])]/wei
    }
  }
}

dim(ynorm.homo)
names(ynorm.homo)

# Decomposicao da SOMA de Quadrados Total
#SSTotal = SS_T + SS_MT + SS_Res

X <- data.frame(ynorm.homo)
dim(X)
n.F <- ncol(ynorm.homo)-3 #numero de variaveis
n.S <- length(unique(X$Sc)) #numero de amostras
T.id <- c("NP", "P", "T", "Pool")
M.id <- c("IGD", "ISD", "OFD", "OPD")
n.T <- length(T.id) #numero de tecidos
n.M <- length(M.id) #number de metodos
mean.TM <- matrix(0,n.F,1)
mean.T <- matrix(0,n.F,1)


SSMT <- SST <- matrix(0,n.F,n.F)
mean.Total <- as.vector(colMeans(X[,4:49])) #media de cada F
for (tis in T.id) {
  for (meth in M.id) {
    idTM <- which(ynorm.homo$Tfc == tis & ynorm.homo$Mfc ==  meth) #posições em que temos T = t e M = m
    idT <- which(ynorm.homo$Tfc == tis) #posições em que temos T = t
    mean.TM <- as.vector(colMeans(X[idTM,4:49])) #media das observações T = t e M = m para cada F
    mean.T <- as.vector(colMeans(X[idT,4:49])) #media das observações T = t para cada F
    SSMT <- SSMT+n.S*(mean.TM-mean.T)%*%t(mean.TM-mean.T)
    SST <- SST+(n.M*n.S)*(mean.T-mean.Total)%*%t(mean.T-mean.Total)
  }
}

dim(SSMT)
dim(SST)

SSTotal <- matrix(0,n.F,n.F)
ids <- c(rep(1:nrow(X)))
for (grus in ids) {
  SSTotal <- SSTotal+(t(X[grus,4:49])-mean.Total)%*%t(t(X[grus,4:49])-mean.Total)
}
dim(SSTotal)

SSRes <- SSTotal-(SSMT+SST)

# MANOVA - Dados aninhados M(T)

netMT <- cbind(ynorm.homo[,1:3],matrix(0,80,46))
dim(netMT)
T.id <- c("NP", "P", "T", "Pool")
M.id <- c("ISD", "IGD", "OPD", "OFD")
for (tis in T.id) {
  for (meth in M.id) {
    idTM <- which(ynorm.homo$Tfc == tis & ynorm.homo$Mfc ==  meth)
    idT <-  which(ynorm.homo$Tfc == tis)
    mean.TM <- as.vector(colMeans(X[idTM,4:49]))
    mean.T <- as.vector(colMeans(X[idT,4:49]))
    difMT <- t(as.matrix(mean.TM-mean.T))
    netMT[which(netMT$Tfc == tis & netMT$Mfc ==  meth),c(4:49)] <- rbind(difMT,difMT,difMT,difMT,difMT)
  }
}


# Analises Multivariadas ----
# PCA ----

pcnet <- prcomp(netMT[,4:49])
pcnet
summary(pcnet)
names(pcnet)

plot(pcnet$x[,1],pcnet$x[,2], xlab = "PC1 (53,73%)", ylab = "PC2 (28,34%)", type = "n", xlim = c(-300,300))
text(pcnet$x[,1],pcnet$x[,2],labels=as.factor(paste0(netMT$T,":",netMT$M)))
plot(pcnet$rotation[,1],pcnet$rotation[,2], xlab = "PC1 (53,73%)", ylab = "PC2 (28,34%)", type = "n")
text(pcnet$rotation[,1],pcnet$rotation[,2])

# Grafos ----

S <- cov.wt(netMT[,4:49], method = "ML")$cov
C <- cov2cor(S)
suffStat <- list(C = C,n = nrow(netMT))
indepTest <- gaussCItest

cpdag <- pc(suffStat, gaussCItest, p = ncol(netMT[,4:49]),
            alpha = 0.01)
nodes(cpdag@graph) <- names(netMT[,4:49])

# Passos necessarios para obter grafos conexos
A5 <- cpdag@graph		
A5 <- as(A5, "matrix")	
g5 <- graph.adjacency(A5, mode='directed')
V(g5)$name <- names(netMT[,4:49]) 

subgraph <- decompose.graph(g5)
clusters(g5)
grafos <- which(clusters(g5)$csize >= 1)

pdf("Grafos_SRM.pdf")
for(i in grafos){
  
  l = layout_with_lgl(subgraph[[i]])
  
  plot(subgraph[[i]], vertex.size = 12, edge.width = 1, layout = l, vertex.color = "white",
       edge.arrow.size = .4, asp = 0, vertex.label.cex = 0.8, vertex.label.color = "black")
}
dev.off()