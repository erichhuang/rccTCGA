## firstPass.R

## Erich S. Huang
## Sage Bionetworks
## Seattle, Washington
## erich.huang@sagebase.org

## REQUIRE
require(synapseClient)
require(annotate)
require(org.Hs.eg.db)
require(Biobase)
require(randomForest)
require(ggplot2)
require(ROCR)
require(reshape)

## LOAD DATA
rccEnt <- loadEntity('syn594913')
rccEset <- rccEnt$objects$eset
rccRnaSeq <- exprs(rccEset)

  # Read in VHL indicator
vhlInd <- read.table('Mutations/vhlIndicator.txt')
rownames(vhlInd) <- vhlInd[ , 1]

colnames(rccRnaSeq) <- sapply(strsplit(colnames(rccRnaSeq), '-'), function(x){
  paste(x[1:3], collapse="-")
  })

## INTERSECT THE DATA
idIntersect <- intersect(rownames(vhlInd), colnames(rccRnaSeq))
vhlSubset <- vhlInd[idIntersect, ]
rccRnaSubset <- rccRnaSeq[ , idIntersect]

## RANDOMLY DIVIDE THE DATA INTO TRAINING AND VALIDATION COHORTS
set.seed(121207)
cohortInd <- sample(1:dim(rccRnaSubset)[2], dim(rccRnaSubset)[2]/2)

trainClass <- vhlSubset[cohortInd, 2]
trainExpress <- rccRnaSubset[ , cohortInd]

validClass <- vhlSubset[-cohortInd, 2]
validExpress <- rccRnaSubset[ , -cohortInd]

## GENERATE THE NEW RANDOM FOREST MODEL PHASE 1
vhlClassModel <- randomForest(t(trainExpress),
                              as.factor(trainClass),
                              ntree = 500,
                              do.trace = 2,
                              importance = TRUE,
                              proximity = TRUE)

## Obtain train hat predictions
trainScoreHat <- predict(vhlClassModel, t(trainExpress), type = 'prob')

## VALIDATE THE MODEL
validScoreHat <- predict(vhlClassModel, t(validExpress), type = 'prob')

# Nonparametric ranksum
validNumeric <- as.numeric(validClass)
rankSum <- wilcox.test(validScoreHat[validNumeric == 0, 2], 
                       validScoreHat[validNumeric == 1, 2])

##### THIS, AND THE OOB ERROR SUGGEST THE MODEL MIGHT BE REFINED, THEREFORE
impMat <- importance(vhlClassModel)
impDF <- as.data.frame(impMat)
impDF <- data.frame(rownames(impDF), impDF)
impDF <- arrange(impDF, -MeanDecreaseAccuracy)

topMdaProbes <- impDF[impDF$MeanDecreaseAccuracy >= 1, 1]
topMdgProbes <- impDF[impDF$MeanDecreaseGini >= 0.05, 1]

selectProbes <- intersect(as.character(topMdaProbes), as.character(topMdgProbes))

    # New training and validation matrices using the 'select' probesets

trainExpress <- rccRnaSubset[selectProbes , cohortInd]
validExpress <- rccRnaSubset[selectProbes , -cohortInd]

## RERUN THE RANDOM FOREST MODEL NOW ON THE SELECT PROBESETS
vhlSelectModel <- randomForest(t(trainExpress),
                              as.factor(trainClass),
                              ntree = 500,
                              do.trace = 2,
                              importance = TRUE,
                              proximity = TRUE)

## Obtain train hat predictions
trainScoreHat <- predict(vhlSelectModel, t(trainExpress), type = 'prob')

## VALIDATE THE MODEL
validScoreHat <- predict(vhlSelectModel, t(validExpress), type = 'prob')

# Nonparametric ranksum
validNumeric <- as.numeric(validClass)
newRankSum <- wilcox.test(validScoreHat[validNumeric == 0, 2], 
                       validScoreHat[validNumeric == 1, 2])


# Visualize results
trainScoreDF <- data.frame(trainClass, trainScoreHat[ , 2])
colnames(trainScoreDF) <- c('vhlMutation', 'classPrediction')
trainBoxPlot <- ggplot(trainScoreDF, aes(factor(vhlMutation), classPrediction)) +
  geom_boxplot() +
  geom_jitter(aes(colour = as.factor(vhlMutation), size = 4)) +
  scale_size(guide = 'none') +
  ggtitle('Training Cohort Boxplot (Model Fit)\n')

validScoreDF <- data.frame(validClass, validScoreHat[ , 2])
colnames(validScoreDF) <- (c('vhlMutation', 'classPrediction'))
validBoxPlot <- ggplot(validScoreDF, aes(factor(vhlMutation), classPrediction)) +
  geom_boxplot() +
  geom_jitter(aes(colour = as.factor(vhlMutation), size = 4)) +
  ggtitle('Validation Cohort Boxplot\n') +
  scale_size(guide = 'none')

validDensPlot <- ggplot(validScoreDF, aes(classPrediction, fill = factor(vhlMutation))) +
  geom_density(alpha = 0.3) +
  ggtitle('Validation Cohort Density Plot\n') +
  scale_size(guide = 'none')


## EVALUATE VALIDATION MODEL PERFORMANCE
vhlPred <- prediction(as.numeric(validScoreHat[, 2]), as.numeric(validClass))
vhlPerf <- performance(vhlPred, "tpr", "fpr")
vhlAUC <- performance(vhlPred, "auc")

# FIND YOUDEN'S J POINT AND OPTIMAL SENSITIVITY AND SPECIFICITY
vhlSsPerf <- performance(vhlPred, "sens", "spec")
youdensJ <- vhlSsPerf@x.values[[1]] + vhlSsPerf@y.values[[1]] - 1
jMax <- which.max(youdensJ)
optCut <- vhlPerf@alpha.values[[1]][jMax]

optSens <- unlist(vhlSsPerf@x.values)[jMax]
optSpec <- unlist(vhlSsPerf@y.values)[jMax]
auc <- unlist(vhlAUC@y.values)

## CREATE A ROC CURVE USING GGPLOT
dfPerf <- as.data.frame(cbind(unlist(vhlPerf@x.values),
                              unlist(vhlPerf@y.values)))
colnames(dfPerf) <- c("FalsePositiveRate", "TruePositiveRate")

rocCurve <- ggplot(dfPerf, aes(FalsePositiveRate, TruePositiveRate)) +
  geom_line() + 
  geom_abline(slope = 1, colour = "red") +
  ggtitle('Figure 9: Validation Cohort ROC Curve\n') +
  ylab("False Positive Rate") +
  xlab("True Positive Rate") +
  annotate('text', label = paste('AUC =', signif(auc, 3), sep = ' '), 
           x = 0.75, y = 0.25, size = 5, colour = "red") +
             annotate('text', label = paste('Sens =', signif(optSens, 2), sep = ' '), 
                      x = 0.75, y = 0.20, size = 5, colour = "red") +
                        annotate('text', label = paste('Spec =', signif(optSpec, 2), sep = ' '), 
                                 x = 0.75, y = 0.15, size = 5, colour = "red")

## 2D DENSITY PLOT VISUALIZATION
svdValid <- svd(validExpress)

densDF <- data.frame(validScoreHat[ , 2], validClass, svdValid$v)
colnames(densDF) <- c('prediction', 'vhlMutation', paste('pc', 1:19, sep = ''))
contourPlot <- ggplot(densDF, aes(x = pc1, y = prediction)) +
  geom_point(aes(colour = as.factor(vhlMutation), 
                 size = prediction,
                 shape = as.factor(vhlMutation))) +
                   geom_density2d() +
                   ggtitle('VHL Landscape: Validaton Cohort\n') +
                   scale_size(guide = 'none')

