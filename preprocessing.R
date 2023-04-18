library(tidyverse)
library(mice)

#devtools::install_github("anujkhare/iregnet")
library(glmnet)
library(survival)
library(iregnet)

# installing penAFT

#install.packages("devtools")
# library(devtools)
# devtools::install_github("ajmolstad/penAFT")
library(penAFT)


data = read.csv('METABRIC_RNA_Mutation.csv')
head(data)
dim(data)

sum(is.na(data))

unique(data$cancer_type) # drop the patient with sarcoma

sum(data$cancer_type == "Breast Sarcoma")


unique(data$cancer_type_detailed)

sum(data$cancer_type_detailed == 'Breast Invasive Ductal Carcinoma')

data = data[data$cancer_type_detailed == 'Breast Invasive Ductal Carcinoma',]
data = data[data$pam50_._claudin.low_subtype != 'NC',]
# drop: patient_id
data = data %>% select(-c('patient_id', 'cancer_type', 'cancer_type_detailed', 'cohort', 'overall_survival'))


sum(is.na(data))
# select only the ones for which we know the surgery type

unique(data$type_of_breast_surgery)

data_sel = data[data$type_of_breast_surgery != "",]


# drop the patients with stage 0 or stage 4 tumor
#data_sel = data_sel[data_sel$tumor_stage != 0 & data_sel$tumor_stage != 4,]

sum(is.na(data_sel))

unique(data$overall_survival_months)

sum(is.na(data_sel$overall_survival_months))


data_sel = Filter(function(x)(length(unique(x))>1), data_sel)
data_BIDC = data_sel

names(data_sel)
 

print(names(data_BIDC))

sum(is.na(data_BIDC[,26:514]))

sum(is.na(data_BIDC[,26:514]))

# 0 - breast conserving, 1 - mastectomy
data_BIDC$treatment = 1*(data_BIDC$type_of_breast_surgery == "MASTECTOMY")
unique(data_BIDC$type_of_breast_surgery)

data_BIDC$primary_tumor_laterality[data_BIDC$primary_tumor_laterality == ""] = NA
data_BIDC$inferred_menopausal_state[data_BIDC$inferred_menopausal_state == ""] = NA
data_BIDC$er_status_measured_by_ihc[data_BIDC$er_status_measured_by_ihc == ""] = NA


sum(is.na(data_BIDC$primary_tumor_laterality))
sum(is.na(data_BIDC$inferred_menopausal_state))
sum(is.na(data_BIDC$er_status_measured_by_ihc))

sum(is.na(data_BIDC[,26:514]))

names(data_BIDC)[26:514]

categorical_vars <- c("cellularity", "pam50_._claudin.low_subtype",
                      "neoplasm_histologic_grade", "tumor_other_histologic_subtype",
                      "integrative_cluster", "X3.gene_classifier_subtype")

other_categorical = c("inferred_menopausal_state", 
                      "primary_tumor_laterality",
                      "pr_status",
                      "tumor_stage",
                      "her2_status",
                      "er_status_measured_by_ihc",
                      "radio_therapy", "hormone_therapy", "chemotherapy")

data_mice <- data_BIDC %>% mutate(across(all_of(unlist(c(categorical_vars, other_categorical))), factor))

data_mice = data_mice[,1:24]
#sum(is.na(data_mice))
names(data_mice)

mice_obj <- mice(data = data_mice, m = 1)
data_full <- complete(mice_obj)

data_BIDC[,1:24] = data_full

data_BIDC$tumor_stage = 0*(data_BIDC$tumor_stage == 0) + 1*(data_BIDC$tumor_stage == 1) + 2*(data_BIDC$tumor_stage == 2) + 3*(data_BIDC$tumor_stage == 3) + 4*(data_BIDC$tumor_stage == 4)

sum(data_BIDC$tumor_stage == 0)

sum(data$tumor_stage == 4,na.rm = TRUE)

# Drop patients with stage 0 or stage 4 tumor

data_BIDC = data_BIDC[data_BIDC$tumor_stage != 0 & data_BIDC$tumor_stage != 4,]

sum(is.na(data_BIDC))
# 0 - pre, 1 - post menopause
data_BIDC$menopause = 1*(data_BIDC$inferred_menopausal_state == "Post")

unique(data_BIDC$inferred_menopausal_state)

# 0 - left, 1 - right
data_BIDC$tumor_laterality = 1*(data_BIDC$primary_tumor_laterality == "Right")

unique(data_BIDC$tumor_laterality)

data_BIDC$progesterone_status = 1*(data_BIDC$pr_status == "Positive")

unique(data_BIDC$pr_status)

data_BIDC$HER2_status = 1*(data_BIDC$her2_status == "Positive")

unique(data_BIDC$her2_status)

data_BIDC$er_status_ihc = 1*(data_BIDC$er_status_measured_by_ihc == "Positive")

unique(data_BIDC$er_status_measured_by_ihc)

unique(data_BIDC$pam50_._claudin.low_subtype)

names(data_BIDC)

sum(is.na(data_BIDC))
#sum(is.na(data_sel[,515:687]))
#categorical variables: cellularity, chemotherapy, pam50_._claudin.low_subtype, er_status_measured_by_ihc, er_status
# neoplasm_histologic_grade, her2_status_measured_by_snp6, tumor_other_histologic_subtype, hormone_therapy,
# menopause, integrative_cluster, tumor_laterality, progesterone_status, radio_therapy, X3.gene_classifier_subtype,

#numerical variables: age_at_diagnosis, lymph_nodes_examined_positive, tumor_size, tumor_stage

#deleted variables: mutation_count, nottingham_prognostic_index, her2_status_measured_by_snp6



# unique(data_BIDC$tumor_stage)
# unique(data_BIDC$radio_therapy)
# 
# unique(data_BIDC$primary_tumor_laterality)
# sum(data_BIDC$primary_tumor_laterality == "")


# Load the required library
#library(dplyr)

# # Define the gene expression variable names
# gene_expr_mes <- c("brca1", "brca2", "palb2", "pten", "tp53", "st7", "star", "tnk2")
# 
# # Define the mutation variable names
# mutation_names <- c("brca1_mut", "brca2_mut", "palb_mut", "pten_mut")


gene_expr_names = names(data_BIDC)[26:514]

mutation_names = names(data_BIDC)[515:687]

sum(is.na(data_BIDC))

p = ncol(data_BIDC)
for (i in 1:p)
{
  if (sum(is.na(data_BIDC[,i])) > 0)
  {
    print(names(data_BIDC)[i])
    print(sum(is.na(data_BIDC[,i])))
  }
}

names(data_BIDC)[693]
sum(is.na(data_BIDC[,693]))

# Mutation indicators
for (var in mutation_names) {
  data_BIDC[[var]] <- as.integer(data_BIDC[[var]] != 0)
}


# Exclude patients with rare gene mutations  - below 20 in the population

n = nrow(data_BIDC)
include = rep(FALSE, n)
for (i in 1:n)
{
  if (sum(data_BIDC[i,c(608, 611:614, 616:687)]) == 0)
  {
    include[i] = TRUE
  }
}

sum(include)

data_BIDC_common = data_BIDC[include,]

# DROP STAGE 0 AND 4 TUMOR PATIENTS

#data_BIDC_common = data_BIDC_common[data_BIDC_common$tumor_stage != 0 & data_BIDC_common$tumor_stage != 4,]

nrow(data_BIDC_common)

# Exclude patients with rare subtype - NC
#data_BIDC_common = data_BIDC_common[data_BIDC_common$pam50_._claudin.low_subtype != 'NC',]
#nrow(data_BIDC_common)

# delete mutations which did not occur
mutation_names = setdiff(mutation_names, c("hras_mut", "siah1_mut", "smarcb1_mut", "stmn2_mut", "foxo1_mut", "ccnd3_mut", "ctnna1_mut", "klrg1_mut"))

# delete rare mutations
mutation_names = setdiff(mutation_names, names(data_BIDC)[c(608, 611:614, 616:685)])


for (i in 515:687)
{
  #print(sum(data_BIDC[,i]))
  cat(i, ', ', sum(data_BIDC_common[-data_BIDC_common,i]*(data_BIDC$treatment == 0)), '\n')
}


for (i in 515:687)
{
  #print(sum(data_BIDC[,i]))
  cat(i, ', ', sum(data_BIDC_common[,i]*(data_BIDC$treatment == 1)), '\n')
}

sum(data_BIDC_common[,14])

for (i in 515:687)
{
  #print(sum(data_BIDC[,i]))
  cat(i, ', ', sum(data_BIDC[,i]), '\n')
}



for (i in 1:24)
{
  #print(sum(data_BIDC[,i]))
  cat(i, ', ', sum(data_BIDC[,i]*(data_BIDC$treatment == 1)), '\n')
}


colnames(model_matrix)[1:24]

for (i in 1:30)
{
  #print(sum(data_BIDC[,i]))
  cat(i, ', ', sum(model_matrix[,i]*(data_BIDC$treatment == 0)), '\n')
}

# 14, 18, 19

colnames(model_matrix)[14]


# Rename type_of_breast_surgery to treatment
# data_BIDC <- data_BIDC %>%
#   rename(treatment = type_of_breast_surgery)

# Convert the mutation variables, categorical variables, and treatment to factors

categorical_vars <- c("cellularity", "pam50_._claudin.low_subtype",
                      "neoplasm_histologic_grade", "tumor_other_histologic_subtype",
                      "integrative_cluster", "X3.gene_classifier_subtype", 'tumor_stage')

categorical_recoded <- c("menopause", "HER2_status", "progesterone_status", "tumor_laterality", "radio_therapy", "hormone_therapy", "er_status_ihc", "chemotherapy")
#data_BIDC_2 <- mutate_at(data_BIDC, vars(mutation_names, categorical_vars, "treatment"), as.factor)

data_BIDC_common <- data_BIDC_common %>% mutate(across(all_of(unlist(c(categorical_vars))), as.factor))

sum(data_BIDC_common$pam50_._claudin.low_subtype == 'NC')
unique(data_BIDC_common$pam50_._claudin.low_subtype)

unique(data_BIDC$cellularity)
unique(data_BIDC$pam50_._claudin.low_subtype)
unique(data_BIDC$neoplasm_histologic_grade)
sum(is.na(data_BIDC$neoplasm_histologic_grade))

unique(data_BIDC$tumor_other_histologic_subtype)

unique(data_BIDC$integrative_cluster)

unique(data_BIDC$X3.gene_classifier_subtype)


# Create the interaction terms
interaction_terms <- lapply(1:length(gene_expr_names), function(i) {
  if (paste0(gene_expr_names[i], "_mut") %in% mutation_names) {
    return(paste0(gene_expr_names[i], ":", gene_expr_names[i], "_mut"))
  }
})
interaction_terms <- unlist(interaction_terms)

# Define the main effects variables
main_effects <- c("age_at_diagnosis", "lymph_nodes_examined_positive", "tumor_size")

# Create the formula for the model matrix
# formula <- paste(" ~ ", paste(main_effects, collapse = " + "), " + ",
#                  paste(categorical_vars, collapse = " + "), " + ",
#                  paste(gene_expr_names, collapse = " + "), " + ",
#                  paste(interaction_terms, collapse = " + "), " + ",
#                  paste("treatment", "(", paste(main_effects, categorical_vars, gene_expr_names, collapse = " + "), ")", sep = ""),
#                  " -1")

# Create the model matrix
#model_matrix <- model.matrix(as.formula(formula), data = data_BIDC)

#cat(formula, "\n")


vars1 = paste(main_effects, collapse = " + ")
vars2 = paste(categorical_vars, collapse = " + ")
vars3 = paste(categorical_recoded, collapse = " + ")
vars3 = paste(gene_expr_names, collapse = " + ")
vars4 = paste(interaction_terms, collapse = " + ")


vars_all = paste(vars1, vars2, vars3, vars4, sep = " + ")
form = paste(" ~  ", "treatment + ", vars_all, " + treatment:(", vars_all, ")")

model_matrix <- model.matrix(as.formula(form), data = data_BIDC_common)
model_matrix[1,1:10]

dim(model_matrix)

model_matrix[1,1:10]

dim(model_matrix)

###########################################

dataX = model_matrix[,2:ncol(model_matrix)]
logY = log(data_BIDC_common$overall_survival_months)
delta = 1*(data_BIDC_common$death_from_cancer == "Died of Disease")

x = model_matrix
p = ncol(dataX)

#weight.set <- list("w" = c(0, rep(1, p-1)))
set.seed(3)
fit.en.cv <- penAFT.cv(dataX, logY, delta, alpha = 0.5, nlambda = 30, nfolds = 5)

saveRDS(fit.en.cv, file = 'penAFT_fit_alpha_0_6_nlambda_20.RDS')

fit.en.cv$full.fit

fit.en.cv$full.fit[4]

fit.en.cv$lambda.min

fit.en.cv$cv.err.linPred

fit.en.cv$full.fit$lambda

fit.en.cv$cv.err.obj

preds <- penAFT.predict(fit.en.cv, Xnew = dataX, lambda = fit.en$lambda[20])

penAFT.plot(fit.en.cv)

fit.en = penAFT(dataX, logY, delta, alpha = 0.5, nlambda = 10)

p = penAFT.trace(fit.en)

fit.en.cv$full.fit$lambda

log10(0.2247)

penAFT.plot(fit.en)
fit.en$full.fit$lambda

log10(0.51)

############################################################
###### Change after fixing package #########################

lambda = fit.en.cv$full.fit$lambda
cv.err.linPred = fit.en.cv$cv.err.linPred
lambda[which(cv.err.linPred == min(cv.err.linPred))]

best.ind = which.min(cv.err.linPred)
lambda.min = lambda[which.min(cv.err.linPred)]
lambda.gehan = lambda.min

beta.gehan = fit.en.cv$full.fit$beta[,best.ind]
############################################################
str(fit.en.cv)


str(penAFT.plot(fit.en.cv))


p = ncol(dataX)
for (i in 1:p)
{
  if (sum((dataX[,i])) == 0)
  {
    print(i)
  }
}

colnames(dataX)[651]

cor(data_BIDC$tumor_size, data_BIDC$tumor_stage, use = "complete.obs")



####################################
# PENALIZED COX

y = Surv(data_BIDC_common$overall_survival_months, delta)
x = dataX
#out.cox = cox.path(x, y,nlambda = 10, alpha = 0.5)

# Elastic net Cox
out.cox = glmnet(x, y, family = "cox", nlambda = 20, alpha = 0.5)

set.seed(4)
out.cox.cv = cv.glmnet(x, y, family = "cox", nlambda = 20, alpha = 0.5, nfolds = 5)

lambda.cox = out.cox.cv$lambda.min

out.cox.cv$lambda.1se

sum(abs(out.cox$beta) < 1e-10) / length(out.cox$beta)

beta.cox = out.cox.cv$glmnet.fit$beta

dim(x%*%beta.cox)
dim(x%*%beta.gehan)


########################################
####### OUT OF SAMPLE PERFORMANCE ######
########################################

get.concordance = function(pred_test, truth_test, death)
{
  nvalid = length(pred_test)
  agree.count = 0
  pair.count = 0
  for (i in 2:nvalid)
  {
    for (j in 1:(i-1))
    {
      pair.count = pair.count + death[j]
      agree.count = agree.count + death[j]*((pred_test[i] >= pred_test[j]) == (truth_test[i] >= truth_test[j]))
    }
  }
  
  concord = agree.count / pair.count
  return(concord)
}

mut = dataX[,633]
dataX_2 = dataX[mut == 0,]
dataX_2 = dataX[mut == 0, setdiff(1:ncol(dataX), c(633))]

ndata = nrow(data_BIDC_common)
nvalid = floor(0.3*nrow(data_BIDC_common))

#ndata = nrow(dataX)
#nvalid = floor(0.3*nrow(dataX))

#set.seed(4)

set.seed(4)
index = sample(1:ndata, size = nvalid, replace = FALSE)
test_data = data_BIDC_common[index,]
train_data = data_BIDC_common[-index,]

#en.gehan = penAFT(dataX[-index,], logY[-index], delta[-index], alpha = 0.6, lambda = c(lambda.gehan))

set.seed(8)
en.gehan.cv = penAFT.cv(dataX[-index,], logY[-index], delta[-index], alpha = 0.5, nlambda = 30, nfold = 5, standardize = FALSE)

lambda.gehan = en.gehan.cv$lambda.min

saveRDS(en.gehan.cv, file = 'gehan_train.RDS')
#en.gehan = penAFT(dataX[-index,], logY[-index], delta[-index], alpha = 0.6, nlambda = 10)

x = model_matrix

set.seed(6)
en.cox.cv =  cv.glmnet(x[-index,], y[-index,], family = "cox", alpha = 0.5, nlambda = 30, nfold = 5)

sum(abs(en.gehan$beta))

en.gehan$beta[,1]['treatment']

en.gehan$beta[1,1]


which(colnames(dataX) == 'treatment')

sum(abs(en.cox$beta))

preds.gehan <- penAFT.predict(en.gehan.cv, Xnew = dataX[index,], lambda = lambda.gehan)


dim(en.cox.cv$glmnet.fit$beta)
en.cox.cv$index

beta.cox = en.cox.cv$glmnet.fit$beta[,4]
lambda.cox = en.cox.cv$glmnet.fit$lambda[4]

beta.cox = en.cox.cv$glmnet.fit$beta[,4]

preds.cox <- x[index,]%*%beta.cox

truth_test = logY[index]
death = delta[index]

get.concordance(-preds.cox, truth_test, death)
get.concordance(preds.gehan, truth_test, death)

# Cox - 0.7155
# Semiparametric AFT - 0.6437


###################################################################





####################################################################
####################################################################
####################################################################



# Parametric AFT
set.seed(4)
aft_param = cv.iregnet(x = dataX, y = y, family = "gaussian", num_lambda = 20, alpha = 0.5, nfolds = 5L)

set.seed(4)
valid.ind = sample(1:nrow(dataX), floor(0.7*nrow(dataX)), replace = FALSE)
trainX = dataX[valid.ind,]
trainY = y[valid.ind,]

validX = (dataX[-index,])[-valid.ind,]
validY = (logY[-index])[-valid.ind]
validDeath = (delta[-index])[-valid.ind]

aft_param_path = iregnet(x = trainX, y = trainY, family = "gaussian", num_lambda = 20, alpha = 0.5, intercept = FALSE)

for (i in 1:20)
{
  beta = aft_param_path$beta[,i]
  pred = validX%*%beta
  conc = get.concordance(pred_test = pred, truth_test = validY, death = validDeath)
  print(conc)
}

lambda.aft = aft_param_path$lambda[18]
aft_param_fin = iregnet(x = dataX, y = y, family = "gaussian", num_lambda = 20, alpha = 0.5, intercept = FALSE)

aft_param_cv = cv.iregnet(x = dataX[-index,], y = y[-index], family = "loglogistic", alpha = 0.5, num_lambda = 20, nfolds = 5L, intercept = FALSE)

lam = aft_param_cv$lambda

num.check  = 300
lam = exp(seq(from = log(5.0), to = log(1e-5), length = num.check))

conc.vec = vector(length = num.check)
for (i in 1:num.check)
{
# play with intercept
aft_param_fin = iregnet(x = dataX[-index,], y = y[-index], family = "weibull", alpha = 0.5, lambda = lam[i], intercept = TRUE)
pred = predict(aft_param_fin, newx = dataX[index,], lambda = lam[i])
#beta = aft_param_fin$beta[,1]
#pred = dataX[index,]%*%beta
conc = get.concordance(pred_test = pred, truth_test = logY[index], death = delta[index])
conc.vec[i] = conc
}

plot(lam[1:150], conc.vec[1:150], type = 'l')
max(conc.vec)

which.max(conc.vec)

lam.opt = lam[which.max(conc.vec)]
lam.opt

# loglogistic - decent performance with lambda = 1.848593e-05 

# Gaussian - 0.5380261
# Weibull - 0.5283745
# Logistic - 0.5413
# Lognormal (loggaussian) - 0.5317

#####################################################
#####################################################
######## CAUSAL INFERENCE ###########################

x = model_matrix

delta = 1*(data_BIDC_common$death_from_cancer == "Died of Disease")
y = Surv(data_BIDC_common$overall_survival_months, delta)

# PROPENSITY SCORES 

vars1 = paste(main_effects, collapse = " + ")
vars2 = paste(categorical_vars, collapse = " + ")
vars3 = paste(categorical_recoded, collapse = " + ")
vars3 = paste(gene_expr_names, collapse = " + ")
vars4 = paste(interaction_terms, collapse = " + ")


vars_all = paste(vars1, vars2, vars3, vars4, sep = " + ")
form_2 = paste("treatment ~  ", vars_all)

x_prop <- model.matrix(as.formula(form), data = data_BIDC_common)
treatment = data_BIDC_common$treatment

# ridge regression
set.seed(3)
prop.cv  = cv.glmnet(x_prop, treatment, family = "binomial", alpha = 0, nlambda = 30, nfold = 5)
e.vec = predict(prop.cv, newx = x, s = prop.cv$lambda.min, type = 'response')
omega = 1
W.vec = omega / (treatment*e.vec + (1-treatment)*e.vec)

set.seed(4)
cox.cv  = cv.glmnet(x, y, weights = W.vec, family = "cox", alpha = 0.5, nlambda = 30, nfold = 5)
lambda.cox = cox.cv$lambda.min

cox.fit =  glmnet(x, y, family = "cox", alpha = 0.5, lambda = c(lambda.cox))

cox.fit$beta[,1]['treatment']


sum(abs(cox.fit$beta[,1]) > 1e-10)

cox.fit$beta[,1][abs(cox.fit$beta[,1]) > 1e-10]

# Exclude treatment from the penalized predictors



#########################################################
######################################################3
#########################################################


sum(data_BIDC_common$integrative_cluster == 5)
sum(data_BIDC_common$tumor_stage == 4 & data_BIDC_common$treatment == 1)

sum(data_BIDC_common$tumor_stage == 0)

unique(data_BIDC_common$tumor_stage)

# FIRST WE PERFORM COMPLETE DATA ANALYSIS - LATER WE LOOK AT MICE
data_complete = na.omit(data_BIDC)
model_matrix <- model.matrix(as.formula(form), data = data_BIDC_common)



cox.cv  = cv.glmnet(dataX, y, family = "cox", alpha = 0.1, nlambda = 100, nfold = 5)
lambda.cox = cox.cv$lambda.min

cox.fit =  glmnet(dataX, y, family = "cox", alpha = 0.1, lambda = c(lambda.cox))
cox.fit$beta[,1]['treatment']

beta.cox = cox.fit$beta
beta.cox['treatment']

sum(beta.cox!= 0)

survival::survfit(cox.fit, s = 0.05, x = dataX, y = y)

plot(survival::survfit(cox.fit, s = 0.05, x = x, y = y))

plot(survival::survfit(cox.fit, s = 0.05, x = dataX[-index,], y = y[-index], newx = dataX[index,]))

plot(survival::survfit(cox.fit, s = 0.05, x = dataX[-index,], y = y[-index], newx = dataX))

plot(survival::survfit(cox.fit, s = 0.05, x = dataX[-index,], y = y[-index], newx = dataX[1:5,]))

plot(survival::survfit(cox.fit, s = 0.05, newx = dataX[1:3,]))


scurve = survival::survfit(cox.fit, s = lambda.cox, x = dataX, y = y, newx = dataX[2,])

# USE THESE TO COMPUTE THE AREA UNDER THE CURVE FOR EACH INDIVIDUAL AND THEN AVERAGE IT OUT !!!!
# OR AVERAGE IT OUT BEFOREHAND - I THINK THE EVENT TIMES SHOULD BE THE SAME

beta.cox['treatment']

sum(beta.cox != 0)

scurve$surv
tim2 = scurve$time

for (i in 1:)

survfit(cox.fit, newdata = dataX)



####################################################
aft_param$plot.data

aft_param$selected


dim(aft_param$beta)

aft_param$lambda
aft_param$selected['min']

  
plot(aft_param)

x%*%(array(aft_param$beta[,1]))

# aft_param = iregnet(x, y, family = "gaussian", num_lambda = 20, alpha = 0.6)

aft_param$lambda

dim(aft_param$beta)

sum(abs(aft_param$beta) < 1e-10) / length(aft_param$beta)

summary(aft_param)

#pred_test_cox <- predict(out.cox, newdata=x, type='response', se=FALSE)

data("neuroblastomaProcessed", package = "penaltyLearning")
X = neuroblastomaProcessed$feature.mat
Y = neuroblastomaProcessed$target.mat
cv_fit <- cv.iregnet(x, y, family = "gaussian", nfolds = 5L)

class(X)
class(x)
