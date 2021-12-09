##############################################################
# Initialisation
##############################################################
# Loading all data files into variables

PAll_k10 = "stat_tab_PARK_S1_rois_k10_p001.csv"
PAll_k10 = read.table(PAll_k10,header=TRUE,sep=";",dec=".")

PAll_k10_evp = "signal_per_event_tab_PARK_S1_rois_k10_p001.txt"
PAll_k10_evp = read.table(PAll_k10_evp,header=TRUE,sep=";",dec=".")

#PS1_s = "stat_tab_PARK_S1.txt"
#PS1_s = read.table(PS1_s,header=TRUE,sep=";",dec=".")
#
#PS1_evp = "signal_per_event_tab_PARK_S1.txt"
#PS1_evp = read.table(PS1_evp,header=TRUE,sep=";",dec=".")
#
#PS2_s = "stat_tab_PARK_S2.txt"
#PS2_s = read.table(PS2_s,header=TRUE,sep=";",dec=".")
#if (length(PS2_s$T)>length(PS1_s$T)) PS2_s <- PS2_s[1:(length(PS2_s$T)/2),]

#PS2_evp = "signal_per_event_tab_PARK_S2.txt"
#PS2_evp = read.table(PS2_evp,header=TRUE,sep=";",dec=".")
#if (length(PS2_evp$roi)>length(PS1_evp$roi)) PS2_evp <- PS2_evp[1:(length(PS2_evp$roi)/2),]


#PS1_s$Condition <- paste(PS1_s$contrast, PS1_s$roi, sep="+")
#PS2_s$Condition <- paste(PS2_s$contrast, PS2_s$roi, sep="+")

PAll_k10$Condition <- paste(PAll_k10$contrast, PAll_k10$roi, sep="+")

##############################################################
## Merging all variables for facilitation purposes
##id_gen = as.factor(substr(PS1_s$id, 1, 6))
##PS1_s = cbind(id_gen, PS1_s)
##id_gen = as.factor(substr(PS2_s$id, 1, 6))
##PS2_s = cbind(id_gen, PS2_s)
##PS1_s$session = as.factor(PS1_s$session)
##PS2_s$session = as.factor(PS2_s$session)

#PARKGAME = rbind(PS1_s,PS2_s)
#subj = as.factor(substr(PARKGAME$id, 1, 6))
#net = as.factor(substr(PARKGAME$roi,10,18))
#PARKGAME = cbind(subj,PARKGAME,net)
#
#PARKGAME$session = as.factor(PARKGAME$session)
#PARKGAME$Condition = as.factor(PARKGAME$Condition)

PARKGAME = PAll_k10
subj = as.factor(PARKGAME$id)
net = as.factor(substr(PARKGAME$roi,10,21))
PARKGAME = cbind(subj,PARKGAME,net)

PARKGAME$session = as.factor(PARKGAME$session)
PARKGAME$Condition = as.factor(PARKGAME$Condition)

# Libraries
library(reshape)
library(inferr)

library(ggplot2) # pour graphes
library(gridExtra) # pour graphes
library(cowplot) # pour regrouper les graphes
library(xlsx) # pour exporter
library(MASS) # Pour Wilcoxon

#nPARK = melt(PARKGAME, id.vars = c("subj","session","Condition"), measure.vars = c("value","pval","pvalC"))
id_S2 = PARKGAME$subj[PARKGAME$session == 2]
PARKGAME = subset(PARKGAME,is.element(PARKGAME$subj, id_S2) == TRUE)
PARKGAME$subj <- droplevels(PARKGAME$subj)

PARKGAME = subset(PARKGAME,PARKGAME$pval <= 0.001)
PARKGAME$roi <- droplevels(PARKGAME$roi)
PARKGAME$Condition <- droplevels(PARKGAME$Condition)

#IMAGINARY_con = subset(PARKGAME,PARKGAME$contrast == "IMAGINARY_Left" | PARKGAME$contrast == "IMAGINARY_Right")
#IMAGINARY_con$contrast <- droplevels(IMAGINARY_con$contrast)
#REAL_con = subset(PARKGAME,PARKGAME$contrast == "REAL_Left" | PARKGAME$contrast == "REAL_Right")
#REAL_con$contrast <- droplevels(REAL_con$contrast)

#IMAGINARY_net = subset(PARKGAME,PARKGAME$net == "Main_spe_RIGHT_IMA" | PARKGAME$net == "Main_spe_LEFT_IMAG")
#IMAGINARY_net$net <- droplevels(IMAGINARY_net$net)
#REAL_net = subset(PARKGAME,PARKGAME$net == "Main_spe_RIGHT_REA" | PARKGAME$net == "Main_spe_LEFT_REAL")
#REAL_net$net <- droplevels(REAL_net$net)
#IMAREAL_net = subset(PARKGAME,PARKGAME$net == "Main_conj_RIGHT_IM" | PARKGAME$net == "Main_conj_LEFT_IMA")
#IMAREAL_net$net <- droplevels(IMAREAL_net$net)


IMAGINARY_con = subset(PARKGAME,PARKGAME$contrast == "IMAGINARY_Left - Rest" | PARKGAME$contrast == "IMAGINARY_Right - Rest")
IMAGINARY_con$contrast <- droplevels(IMAGINARY_con$contrast)
REAL_con = subset(PARKGAME,PARKGAME$contrast == "REAL_Left - Rest" | PARKGAME$contrast == "REAL_Right - Rest")
REAL_con$contrast <- droplevels(REAL_con$contrast)

IMAGINARY_net = subset(PARKGAME,PARKGAME$net == "Mask_RIGHT_I" | PARKGAME$net == "Mask_LEFT_IM")
IMAGINARY_net$net <- droplevels(IMAGINARY_net$net)
REAL_net = subset(PARKGAME,PARKGAME$net == "Mask_RIGHT_R" | PARKGAME$net == "Mask_LEFT_RE")
REAL_net$net <- droplevels(REAL_net$net)
IMAREAL_net = subset(PARKGAME,PARKGAME$net == "Conj_RIGHT_R" | PARKGAME$net == "Conj_LEFT_RE")
IMAREAL_net$net <- droplevels(IMAREAL_net$net)

#############################################################################################################
nPARK = melt(PARKGAME, id.vars = c("subj","session","Condition"), measure.vars = "value", na.rm = TRUE)
nPARK = nPARK[order(nPARK$session), ]
net_table = melt(PARKGAME, id.vars = c("subj","session","net"), measure.vars = "value", na.rm = TRUE)

infer_levene_test(nPARK, value, group_var=session)


nIMAGINARY_net = melt(IMAGINARY_net, id.vars = c('subj','session','Condition'), measure.vars = "value")
nREAL_net = melt(REAL_net, id.vars = c('subj','session','Condition'), measure.vars = "value")
nIMAREAL_net = melt(IMAREAL_net, id.vars = c('subj','session','Condition'), measure.vars = "value")

infer_levene_test(nIMAGINARY_net, value, group_var=session)
infer_levene_test(nREAL_net, value, group_var=session)
infer_levene_test(nIMAREAL_net, value, group_var=session)


nIMAGINARY_con = melt(IMAGINARY_con, id.vars = c('subj','session','roi'), measure.vars = "value")
nREAL_con = melt(REAL_con, id.vars = c('subj','session','roi'), measure.vars = "value")

infer_levene_test(nIMAGINARY_con, value, group_var=session)
infer_levene_test(nREAL_con, value, group_var=session)

## ANOVA on these new good-looking tables value~session
an.all <- with(nPARK, aov(value~subj*session))

an.IMA_con <- with(IMAGINARY_con, aov(value~subj*session))
an.REAL_con <- with(REAL_con, aov(value~subj*session))
an.REAL_net <- with(REAL_net, aov(value~subj*session))
an.IMA_net <- with(IMAGINARY_net, aov(value~subj*session))
an.net_all <- with(net_table, aov(value~subj*session))
 
plot(nPARK$value~nPARK$subj)
plot(nPARK$value~nPARK$session)
summary(aov(value~session, data = nPARK))
