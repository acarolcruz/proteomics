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

# DDA TIC Peptideos ----

#carregar os dados em formato .RData

d$gr <- as.factor(d$gr)
d$T <- as.factor(d$T)

# Analise descritiva

ggplot(d,aes(gr,log10(y),color=T))+ geom_boxplot()+
  labs(title="DDA TIC - Quantificação Peptídeo - log10(y)", x = "Método(Tecido)")

ggplot(d,aes(gr,y,color=T))+ geom_boxplot()+
  labs(title="DDA TIC - Quantificação Peptídeo - y", x = "Método(Tecido)")

ggplot(d,aes(x = M, color = T, group =T, y = y)) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line") +
  labs(title="DDA TIC - Perfis - y", x = "Método")

ggplot(d,aes(x = M, color = T, group =T, y = log10(y))) +
  stat_summary(fun.y = mean, geom = "point") +
  stat_summary(fun.y = mean, geom = "line") +
  labs(title="DDA TIC - Perfis - log10(y)", x = "Método")


# Modelo de Significância  -------

ynorm <- log10(d$y)

#Base de analise
drms <- cbind(d,ynorm)

names(drms)
str(drms)
names(drms) <- c("gr","u1","u2","u3","F", "S", "T", "M","y","y_norm")

# Grafico boxplot das quantificações normalizadas
ggplot(drms,aes(gr,y_norm,color=as.factor(T)))+ geom_boxplot()+
  theme(legend.position="none")+
  labs(title="DDA noNorm - Peptide Quantification - y_Normalized", x = "Method(Tissue)")

# Renomenando y_norm para res
head(drms)
names(drms) <- c("gr","u1","u2","u3","F", "S", "T", "M","y","res")
str(drms)

# Estimação dos valores-p para cada constraste e para cada variável
pvalues2by2 <- list()
anova_pvalues <- c()
for (J in levels(drms$F)) {
  drmsJ <- drms[which(drms$F == J),]
  drmsJ$gr <- as.factor(drmsJ$gr)
  drmsJ <- drmsJ[complete.cases(drmsJ),] #retira os Na's
  fitresJ <- tryCatch({
    gls(res ~ T/M, weights = varIdent(form = ~ 1|gr), #considera erros heterocedasticos dados por FTM
        control = glsControl(maxIter = 100, opt = "optim"), na.action = na.exclude, data = drmsJ)
  }, error = function(e) {
    NULL
  })
  if (!is.null(fitresJ)) {
    aovJ <- anova(fitresJ)
    anova_pvalues <- rbind(anova_pvalues, c(J, aovJ$"p-value"))
    refJ <- ref_grid(fitresJ, mode = "df.error", data = drmsJ, nesting = "M %in% T")
    methJ <- emmeans(refJ, "M")
    meth2J <- pairs(methJ, by = "T", reverse = TRUE)
    summJ <- summary(meth2J)
    contrasts <- paste0(summJ$T, ": ", summJ$contrast)
    for (contrastId in 1:length(contrasts)) {
      pvalues2by2[[contrasts[contrastId]]] <- rbind(pvalues2by2[[contrasts[contrastId]]], c(J,summJ$estimate[contrastId], summJ$p.value[contrastId]))
      colnames(pvalues2by2[[contrasts[contrastId]]]) <- c("Peptideo","Estimativa", "p-valor")
    }
  }
}
colnames(anova_pvalues) <- c("Peptideo", "Intercepto", "T", "T:M")
dim(anova_pvalues)


# Grafico de vulcão ----

plot(rep(1:12607), -log10(as.numeric(anova_pvalues[,4])))
min(-log10(as.numeric(anova_pvalues[,4])))

anova_pvaluesMT <- as.numeric(anova_pvalues[,4])

# valores-p corrigidos por FDR
p.adjust.MT <- p.adjust(anova_pvaluesMT, method = "fdr")

plot(anova_pvaluesMT, p.adjust.MT)
plot(rep(1:12607), -log10(p.adjust.MT))

names(pvalues2by2)
p.bf <- -log10(0.05/((length(anova_pvalues)/4)*6)) #ponto de corte

pdf("DDA_peptideo_TIC_siginificance.pdf", width = 12, height = 8)

par(mfrow=c(4,6))
for (L in 1:6){ 
  log10.pvalue <- log10(as.numeric(pvalues2by2[[L]][,"p-valor"]))
  estimate <- as.numeric(pvalues2by2[[L]][,"Estimativa"])
  plot(estimate, -log10.pvalue, main = names(pvalues2by2)[L], col = ifelse(-log10.pvalue > p.bf & abs(estimate) > 2, "red","black"))
  abline(h = p.bf)
  abline(v = c(-2,2))
}

for (L in 7:12){
  log10.pvalue <- log10(as.numeric(pvalues2by2[[L]][,"p-valor"]))
  estimate <- as.numeric(pvalues2by2[[L]][,"Estimativa"])
  plot(estimate, -log10.pvalue, main = names(pvalues2by2)[L], col = ifelse(-log10.pvalue > p.bf & abs(estimate) > 2, "red","black"))
  abline(h = p.bf)
  abline(v = c(-2,2))
}

for (L in 13:18){
  log10.pvalue <- log10(as.numeric(pvalues2by2[[L]][,"p-valor"]))
  estimate <- as.numeric(pvalues2by2[[L]][,"Estimativa"])
  plot(estimate, -log10.pvalue, main = names(pvalues2by2)[L], col = ifelse(-log10.pvalue > p.bf & abs(estimate) > 2, "red","black"))
  abline(h = p.bf)
  abline(v = c(-2,2))
}

for (L in 19:24){
  log10.pvalue <- log10(as.numeric(pvalues2by2[[L]][,"p-valor"]))
  estimate <- as.numeric(pvalues2by2[[L]][,"Estimativa"])
  plot(estimate, -log10.pvalue, main = names(pvalues2by2)[L], col = ifelse(-log10.pvalue > p.bf & abs(estimate) > 2, "red","black"))
  abline(h = p.bf)
  abline(v = c(-2,2))
}
par(mfrow=c(1,1))

dev.off()

# Gera lista de peptideos significantes
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

write_xlsx(as.data.frame(pep_signif), "Resultados/DDA-TIC-Peptideos-signif.xlsx")

# Tabelas baseadas no arquivo pep_signif

pep_signif <- as.data.frame(pep_signif)

# Peptideos identificados como significantes em algum contraste
length(unique(pep_signif$Peptideo))

# Quantidade de peptídeos significantes em cada contraste

table(pep_signif$V1)

# Quantas vezes cada peptideo foi considerado significativo
sort(table(pep_signif$Peptideo), decreasing =  TRUE)

sum(table(pep_signif$Peptideo) == 1)
sum(table(pep_signif$Peptideo) == 2)
sum(table(pep_signif$Peptideo) == 3)

# Padronizando os resíduos para realizaçao das analises multivariadas
# Considerando somente os peptideos presentes em todas as observações

df <- as.data.frame(table(d$F[!is.na(d$y)]))
variaveis <- df$Var1[df$Freq == 80]  
dcomp <- d[d$F %in% variaveis,]
dim(dcomp)

names(dcomp) <- c("grc","u1c","u2c","u3c","Fc","Sc","Tc", "Mc","yc")
str(dcomp)

dcomp$resc <- log10(dcomp$yc)
names(dcomp) <- c("grc","u1c","u2c","u3c","Fc","Sc","Tc", "Mc","yc","resc")

#Base de analise
drcomp <- dcomp

# Desvios-padrao
sdresc <- c()
for (J in variaveis) {
  drcJ <- drcomp[which(drcomp$Fc==J),]
  dru2 <- as.factor(drcJ$u2c)
  fitrescJ <- gls(resc ~ Tc/Mc, weights = varIdent(form = ~ 1|dru2),
                control = glsControl(maxIter = 100, opt = "optim"), na.action = na.exclude, data = drcJ)
  sdrcJ <- summary(fitrescJ)$sigma*c(1,coef(fitrescJ$modelStruct$varStruct, uncons = FALSE)) #estimativas do desvio padrão de TxM para cada F
  sdresc <- rbind(sdresc, c(J,sdrcJ))
}

dim(sdresc)
head(sdresc)
sdresc <- data.frame(sdresc)
names(sdresc)
str(sdresc)
names(sdresc) <- c("F","NP:IGD","P:IGD","T:IGD","pool:IGD",
                 "NP:ISD","P:ISD","T:ISD","pool:ISD",
                 "NP:OFD","P:OFD","T:OFD","pool:OFD",
                 "NP:OPD","P:OPD","T:OPD","pool:OPD")

sdresc$F <- as.factor(sdresc$F)

# Selecionando as colunas F, Sc, Tc, Mc, resc
ynorm <- drcomp[,c(5,6,7,8,10)] 

# Pesos dados pelos sigmas estimados para TxM para cada F
weights <- sdresc

ynorm2 <- reshape(ynorm, direction = "wide", timevar = "Fc", idvar = c("Tc","Mc", "Sc"))

ynorm.homo <- ynorm2
tec <- c("NP", "P", "T", "pool")
met <- c("ISD", "IGD", "OPD", "OFD")

# Padronizando os residuos para cada F
for (teci in tec) {
  for (meti in met) {
    ids <- which(ynorm2$Tc == teci & ynorm2$Mc ==  meti)
    for (var in 1:4824) {
      wei <- as.numeric(as.character(weights[var, paste0(teci, ":", meti)]))
      ynorm.homo[ids,paste0("resc.", weights$F[var])] <- ynorm.homo[ids, paste0("resc.",
                                                                                weights$F[var])]/wei
    }
  }
}

dim(ynorm.homo)
names(ynorm.homo)

# ASCA ----

# Descomposição da SOMA de Quadrados 
# SSTotal = SS_Tissue + SS_Method(Tissue) + SS_Residual

X <- data.frame(ynorm.homo)
dim(X)
n.F <- ncol(ynorm.homo)-3 #numero de variaveis
n.S <- length(unique(X$S)) #number de amostras (5)
T.id <- c("NP", "P", "T", "pool")
M.id <- c("ISD", "IGD", "OPD", "OFD")
n.T <- length(T.id) #numero de tecidos
n.M <- length(M.id) #numero de metodos
mean.TM <- matrix(0,n.F,1)
mean.T <- matrix(0,n.F,1)


SSMT <- SST <- matrix(0,n.F,n.F)
mean.Total <- as.vector(colMeans(X[,4:4827])) #media de cada F
for (tis in T.id) {
  for (meth in M.id) {
    idTM <- which(ynorm.homo$Tc == tis & ynorm.homo$Mc ==  meth) #posições em que temos T = t e M = m
    idT <- which(ynorm.homo$Tc == tis) #posições em que temos T = t
    mean.TM <- as.vector(colMeans(X[idTM,4:4827])) #media das observações T = t e M = m para cada F
    mean.T <- as.vector(colMeans(X[idT,4:4827])) #media das observações T = t para cada F
    SSMT <- SSMT+n.S*(mean.TM-mean.T)%*%t(mean.TM-mean.T)
    SST <- SST+(n.M*n.S)*(mean.T-mean.Total)%*%t(mean.T-mean.Total)
  }
}

head(SSMT)
dim(SSMT)
head(SST)
dim(SST)

SSTotal <- matrix(0,n.F,n.F)
ids <- c(rep(1:nrow(X)))
for (grus in ids) {
  SSTotal <- SSTotal + (t(X[grus,4:4827]) - mean.Total)%*%t(t(X[grus,4:4827]) - mean.Total)
}

head(SSTotal)
dim(SSTotal)

SSRes <- SSTotal - (SSMT+SST)
head(SSRes)


# MANOVA Dados aninhados M(T)

netMT <- cbind(ynorm.homo[,1:3],matrix(0,80,4824))
dim(netMT)

T.id <- c("NP", "P", "T", "pool")
M.id <- c("ISD", "IGD", "OPD", "OFD")
for (tis in T.id) {
  for (meth in M.id) {
    idTM <- which(ynorm.homo$Tc == tis & ynorm.homo$Mc ==  meth)
    idT <- which(ynorm.homo$Tc == tis)
    mean.TM <- as.vector(colMeans(X[idTM,4:4827]))
    mean.T <- as.vector(colMeans(X[idT,4:4827]))
    difMT <- t(as.matrix(mean.TM-mean.T))
    netMT[which(netMT$Tc == tis & netMT$Mc ==  meth),c(4:4827)] <- rbind(difMT,difMT,difMT,difMT,difMT)
  }
}


# Análises Multivariadas ----

# PCA ----
pcnet <- prcomp(netMT[,4:4827])
summary(pcnet)
plot(pcnet$x[,1],pcnet$x[,2], xlab = "PC1 (23,90%)", ylab = "PC2 (14,31%)", type = "n", xlim = c(-1000,1000), ylim = c(-1000,1000))
text(pcnet$x[,1],pcnet$x[,2],labels = as.factor(paste0(netMT$T,":",netMT$M)))
plot(pcnet$rotation[,1],pcnet$rotation[,2], xlab = "PC1 (23,90%)", ylab = "PC2 (14,31%)", type = "n", xlim = c(-0.1,0.15), ylim = c(-0.3,0.3))
text(pcnet$rotation[,1],pcnet$rotation[,2])


# Grafo ----

S <- cov.wt(netMT[,4:4827], method = "ML")$cov
C <- cov2cor(S)
suffStat <- list(C = C, n = nrow(netMT))
indepTest <-  gaussCItest

cpdag <- pc(suffStat, gaussCItest, p = ncol(netMT[,4:4827]),
            alpha = 0.01)
nodes(cpdag@graph) <- names(netMT[,4:4827])

# Passos necessarios para obtencaos dos grafos conexos
A5 <- cpdag@graph		
A5 <- as(A5, "matrix")	
g5 <- graph.adjacency(A5, mode='directed')	
V(g5)$name <- names(netMT[,4:4827]) 		

subgraph <- decompose(g5)
clusters(g5)
table(clusters(g5)$csize)

pdf("Grafos_DDA_peptideo_TIC.pdf")

for(i in grafos){
  l <- layout_with_lgl(subgraph[[i]])
  
  plot(subgraph[[i]], vertex.size = 12, edge.width = 1, layout = l, vertex.color = "white",
       edge.arrow.size =.4, asp = 0, vertex.label.cex = 0.8, vertex.label.color = "black")
}
dev.off()

# Pode-se realizar previamente a analise de signficancia a normalização pela metodologia proposta neste trabalho, bastando seguir os passos dos experimentos DDA noNorm