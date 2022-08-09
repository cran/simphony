## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE, collapse = TRUE, comment = '#>', fig.retina = 2)

## -----------------------------------------------------------------------------
library('data.table')
library('ggplot2')
library('kableExtra')
library('knitr')
library('limma')
library('precrec')
library('simphony')

## -----------------------------------------------------------------------------
set.seed(44)
featureGroups = data.table(fracFeatures = c(0.75, 0.25), amp = c(0, 0.3))
simData = simphony(featureGroups, nFeatures = 200, interval = 2, nReps = 1, family = 'negbinom')

## -----------------------------------------------------------------------------
kable(simData$abundData[1:3, 1:3])

## -----------------------------------------------------------------------------
kable(simData$sampleMetadata[1:3])

## -----------------------------------------------------------------------------
kable(simData$featureMetadata[149:151, !'dispFunc']) %>%
  kable_styling(font_size = 12)

## -----------------------------------------------------------------------------
fmExample = simData$featureMetadata[feature %in% c('feature_150', 'feature_151')]
dExample = mergeSimData(simData, fmExample$feature)

## -----------------------------------------------------------------------------
dExpect = getExpectedAbund(fmExample, 24, times = seq(0, 48, 0.25))

## ---- fig.width = 6, fig.height = 2.75----------------------------------------
dExample[, featureLabel := paste(feature, ifelse(amp0 == 0, '(non-rhythmic)', '(rhythmic)'))]
dExpect[, featureLabel := paste(feature, ifelse(amp0 == 0, '(non-rhythmic)', '(rhythmic)'))]

ggplot(dExample) +
  facet_wrap(~ featureLabel, nrow = 1) +
  geom_line(aes(x = time, y = log2(2^mu + 1)), size = 0.25, data = dExpect) +
  geom_point(aes(x = time, y = log2(abund + 1)), shape = 21, size = 2.5) +
  labs(x = 'Time (h)', y = expression(log[2](counts + 1))) +
  scale_x_continuous(limits = c(0, 48), breaks = seq(0, 48, 8))

## -----------------------------------------------------------------------------
sampleMetadata = copy(simData$sampleMetadata)
sampleMetadata[, timeCos := cos(time * 2 * pi / 24)]
sampleMetadata[, timeSin := sin(time * 2 * pi / 24)]
design = model.matrix(~ timeCos + timeSin, data = sampleMetadata)

## -----------------------------------------------------------------------------
fit = lmFit(log2(simData$abundData + 1), design)
fit = eBayes(fit, trend = TRUE)
rhyLimma = topTable(fit, coef = 2:3, number = Inf)

## -----------------------------------------------------------------------------
rhyLimma$feature = rownames(rhyLimma)
rhyLimma = merge(data.table(rhyLimma), simData$featureMetadata[, .(feature, amp0)], by = 'feature')

## ---- fig.width = 3.5, fig.height = 3-----------------------------------------
ggplot(rhyLimma) +
  geom_jitter(aes(x = factor(amp0), y = P.Value), shape = 21, width = 0.2) +
  labs(x = expression('Rhythm amplitude ' * (log[2] ~ counts)), y = 'P-value of rhythmicity')

## ---- fig.width = 3, fig.height = 3-------------------------------------------
rocprc = evalmod(scores = -log(rhyLimma$P.Value), labels = rhyLimma$amp0 > 0)
autoplot(rocprc, 'ROC')

