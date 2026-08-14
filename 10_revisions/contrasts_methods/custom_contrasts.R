### CONTRASTS INVESTIGATION ###

rm(list = ls())
.rs.restartR()

# For a few GeoMx analyses, I used the following method for creating linear contrasts to test: 
# Start with a model matrix, denoted mm
# mm is constructed in the following manner: 
# > model.matrix(~condition+patient, data=metadata)
# Subset mm to the observations in the reference level of condition, denoted L1, and collapse to column means: 
# > L1 <- mm[metadata$condition == "L1",] |> colMeans()
# Subset mm to the observations in the experimental level of condition, denoted L2, and collapse to column means: 
# > L2 <- mm[metadata$condition == "L2",] |> colMeans()
# Now, manually compute the comparison as (L2-L1). 

# Note that the column means for the dummy variables controlling L1 and L2 should collapse to 0 and 1.
# The patients, however, are not balanced across these levels, L1 and L2, and therefore the column means
# for those dummy variables will not collapse to 0 and 1. 

# As a result, if a patient is over-represented in the L1 group, then L2-L1 will be a contrast containing
# a negative value for that patient's dummy variable. 

# My main question: is there a way to understand this algebraically? Let us try...

# Let's use an easier example to start... the Ellisen TNBC data. We will look at 
# tumor AOIs, which can be hot or cold, as determined by their immune infiltration.
meta <- openxlsx::read.xlsx("~/Documents/projects/geomx/ellisen-breast/meta_cleaned_v5.xlsx")
rownames(meta) <- meta$aoi_id
counts <- read.csv(file = "~/Documents/projects/geomx/ellisen-breast/counts.csv", row.names = 1)

tumor <- meta[meta$AOI_code == "tumor",]
tumor$type <- ifelse(test = tumor$infiltration_category %in% c("excluded", "desert"), yes = "cold", no = "hot")

# Take a look at the study design:
table(tumor$Patient_number, tumor$type)
#   cold hot
# 1   15   0
# 2    5   4
# 3    6   6
# 4    0   9
# 5    2  10
# 6   10   2
# 7   10   0
# 8    5   3

# We have 8 patients, and they are sampled repeatedly and in an unbalanced way. 
# This is poor study design, but sometimes this is hard to avoid for these projects.

# We are left in a difficult spot for a few reasons.
# First, we hope to account for pseudo-replication. Repeated sampling produces correlated
# technical replicates and this should be recognized. However, since the sampling spans the
# two levels of the condition (i.e. for some patients we have data for both hot and cold)
# this is not purely a blocking variable that can be accounted for with a random intercept
# for patient. Second, even if we wished to do this with a mixed model, tools for RNA-seq, 
# such as DESeq2 and limma do not accommodate doing so. Thrird, another option is to set a 
# fixed effect for the patient variable. This works very well if we have balanced (e.g. paired) 
# sampling. Of course, here, we do not have this. 

# Something problematic is that if we have no balance, except for in a small number of
# patients (maybe even just 1), then the effect size we test will be entirely determined by 
# those patients (as opposed factoring in all patients).

# However, if we do what I described above, then the imbalances are baked into the effect
# size, hence permitting all patients to contribute. This is more of an intuitive 
# understanding though, and I want to see this algebraically. 

# So here is the design matrix:
mm <- model.matrix(~type+Patient_number, data = tumor)

# Contrast for "hot"
hot <- mm[tumor$type == "hot",] |> colMeans()

# Contrast for "cold"
cold <- mm[tumor$type == "cold",] |> colMeans()

# Contrast for (hot-cold)
hot_minus_cold <- hot-cold
hot_minus_cold |> unname() |> round(4) # Examine the values
# [1]  0.0000  1.0000  0.0233  0.0633  0.2647  0.2564 -0.1299 -0.1887 -0.0061

# Contrasts for each patient
pt1 <- mm[tumor$Patient_number == "1",] |> colMeans()
pt2 <- mm[tumor$Patient_number == "2",] |> colMeans()
pt3 <- mm[tumor$Patient_number == "3",] |> colMeans()
pt4 <- mm[tumor$Patient_number == "4",] |> colMeans()
pt5 <- mm[tumor$Patient_number == "5",] |> colMeans()
pt6 <- mm[tumor$Patient_number == "6",] |> colMeans()
pt7 <- mm[tumor$Patient_number == "7",] |> colMeans()
pt8 <- mm[tumor$Patient_number == "8",] |> colMeans()

# Reconstructing hot vs. cold by using the patients: attempt 1
# Here, I use the study design
HOT <- (0*pt1 + 4*pt2 + 6*pt3 + 9*pt4 + 10*pt5 + 2*pt6 + 0*pt7 + 3*pt8)/34
COLD <- (15*pt1 + 5*pt2 + 6*pt3 + 0*pt4 + 2*pt5 + 10*pt6 + 10*pt7 + 5*pt8)/53

HOT_MINUS_COLD <- HOT-COLD
HOT_MINUS_COLD |> unname() |> round(4) # Examine the values
# [1]  0.0000  0.4964  0.0233  0.0633  0.2647  0.2564 -0.1299 -0.1887 -0.0061

# Not quite right, though very close.
(hot_minus_cold-HOT_MINUS_COLD) |> unname() |> round(4)
# [1] 0.0000 0.5036 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000

# The dummy variable for hot is off...

# Attempt 2: admittedly, I used Claude for help.
pt1 <- mm[tumor$Patient_number == "1",] |> colMeans()
pt2_hot <- mm[tumor$Patient_number == "2" & tumor$type == "hot",] |> colMeans()
pt2_cold <- mm[tumor$Patient_number == "2" & tumor$type == "cold",] |> colMeans()
pt3_hot <- mm[tumor$Patient_number == "3" & tumor$type == "hot",] |> colMeans()
pt3_cold <- mm[tumor$Patient_number == "3" & tumor$type == "cold",] |> colMeans()
pt4 <- mm[tumor$Patient_number == "4",] |> colMeans()
pt5_hot <- mm[tumor$Patient_number == "5" & tumor$type == "hot",] |> colMeans()
pt5_cold <- mm[tumor$Patient_number == "5" & tumor$type == "cold",] |> colMeans()
pt6_hot <- mm[tumor$Patient_number == "6" & tumor$type == "hot",] |> colMeans()
pt6_cold <- mm[tumor$Patient_number == "6" & tumor$type == "cold",] |> colMeans()
pt7 <- mm[tumor$Patient_number == "7",] |> colMeans()
pt8_hot <- mm[tumor$Patient_number == "8" & tumor$type == "hot",] |> colMeans()
pt8_cold <- mm[tumor$Patient_number == "8" & tumor$type == "cold",] |> colMeans()

HOT <- (0*pt1 + 4*pt2_hot + 6*pt3_hot + 9*pt4 + 10*pt5_hot + 2*pt6_hot + 0*pt7 + 3*pt8_hot)/34
COLD <- (15*pt1 + 5*pt2_cold + 6*pt3_cold + 0*pt4 + 2*pt5_cold + 10*pt6_cold + 10*pt7 + 5*pt8_cold)/53

HOT_MINUS_COLD <- HOT-COLD
HOT_MINUS_COLD |> unname() |> round(4) # Examine the values
# [1]  0.0000  1.0000  0.0233  0.0633  0.2647  0.2564 -0.1299 -0.1887 -0.0061

# Exactly
(hot_minus_cold-HOT_MINUS_COLD) |> unname() |> round(4)
# [1] 0 0 0 0 0 0 0 0 0

# Reviewing what we did, we basically had marginalize across hot and cold.

# This is a method called OM weighting (I think). This stands for observed marginal 
# weighting. This may exist in emmeans. 
norm <- edgeR::DGEList(counts = counts) |> edgeR::calcNormFactors(method = "upperquartile") |> limma::voom()
Y <- norm$E["TACSTD2", rownames(mm)]
d <- cbind("TROP2" = Y, tumor)
mod <- lm(data = d, formula = TROP2 ~ type + Patient_number)

EMM <- emmeans::emmeans(object = mod, specs = c("type"), weights = "cells")
(EMM@linfct[2,] - EMM@linfct[1,]) |> unname() |> round(4)
# [1]  0.0000  1.0000  0.0233  0.0633  0.2647  0.2564 -0.1299 -0.1887 -0.0061

pairs(EMM, reverse = T)
# contrast   estimate    SE df t.ratio p.value
# hot - cold    0.963 0.182 78   5.293  <.0001
# 
# Results are averaged over the levels of: Patient_number 

sum(hot_minus_cold * coef(mod))
# [1] 0.9631304

multcomp::glht(model = mod, linfct = matrix(data = hot_minus_cold, ncol = length(hot_minus_cold))) |> summary(test = multcomp::univariate())
# Simultaneous Tests for General Linear Hypotheses
# 
# Fit: lm(formula = TROP2 ~ type + Patient_number, data = d)
# 
# Linear Hypotheses:
#          Estimate Std. Error t value Pr(>|t|)    
#   1 == 0   0.9631     0.1820   5.293 1.08e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)

# These results match exactly. Therefore, the method employed by weights = "cells"
# is what are doing. This is not exactly OM weighting. Instead it is cell frequency
# weighting. Setting weights="proportional" would do OM weighting. Nonetheless, it 
# is important that I put my finger on what exactly my contrasts are, and I will 
# make sure to note this in my methods from now on.

# This is what the emmeans docs says about this method: 
# Using "cells" weights gives each prediction the same weight as occurs in the model; 
# applied to a reference grid for a model with all interactions, "cells"-weighted 
# EMMs are the same as the ordinary marginal means of the data.

# Looks reasonable too!
plot(EMM) + 
  ggthemes::theme_gdocs() + 
  ggplot2::labs(title = "TROP2")

# In the paper for which I did this, I used limma-voom, but let's do some comparing 
# with limma and variancePartition, which is the LMM extension of the same framework.
norm <- edgeR::DGEList(counts = counts[rownames(counts) != "NegProbe-WTX", rownames(tumor)]) |> 
  edgeR::calcNormFactors(method = "upperquartile") |> 
  limma::voom(design = mm, plot = T) # This looks quite good.
lvfit <- limma::lmFit(object = norm, design = mm)
fixed_out <- lvfit |> limma::contrasts.fit(contrasts = hot_minus_cold) |> limma::eBayes() |> limma::topTable(number = Inf)
fixed_out$target <- rownames(fixed_out)
plot(x = fixed_out$logFC, y = -1*log10(fixed_out$adj.P.Val), pch = 16, col = "grey")

#renv::install("bioc::variancePartition")
library(variancePartition)
norm <- edgeR::DGEList(counts = counts[rownames(counts) != "NegProbe-WTX", rownames(tumor)]) |> 
  edgeR::calcNormFactors(method = "upperquartile") |> 
  limma::voom(design = mm, plot = T) # This looks quite good.
lvfit <- limma::lmFit(object = norm, design = mm)
fixed_out <- lvfit |> limma::contrasts.fit(contrasts = hot_minus_cold) |> limma::eBayes() |> limma::topTable(number = Inf)
plot(x = fixed_out$logFC, y = -1*log10(fixed_out$adj.P.Val), pch = 16, col = "grey")

y <- edgeR::DGEList(counts = counts[rownames(counts) != "NegProbe-WTX", rownames(tumor)]) |> 
  edgeR::calcNormFactors(method = "upperquartile")

param <- BiocParallel::SnowParam(8, "SOCK", progressbar = T)
v <- variancePartition::voomWithDreamWeights(counts = y,
                                             formula = ~type+(1|Patient_number),
                                             data = tumor,
                                             plot = T, save.plot = T,
                                             BPPARAM = param)
attr(v, 'errors')
fitmm <- variancePartition::dream(exprObj = v,
                                  formula = ~type+(1|Patient_number),
                                  data = tumor,
                                  BPPARAM = param
                                  )
attr(fitmm, 'errors')

mixed_out <- variancePartition::eBayes(fit = fitmm) |> 
  variancePartition::topTable(coef = "typehot", number = Inf)
mixed_out$target <- rownames(mixed_out)

plot(x = mixed_out$logFC, y = -1*log10(mixed_out$adj.P.Val), pch = 16, col = "grey")

#pdf(file = "~/Desktop/cell-weighted_vs_mixed.pdf", width = 6, height = 6)
plot(y = mixed_out[rownames(fixed_out),]$logFC, x = fixed_out$logFC, 
     ylab = "mixed effects logFC", xlab = "fixed effects logFC", main = "Cell-Weighted Fixed Effects Estimates\nvs.\nMixed Effects Estimates", cex = 0.2)
abline(0, 1, col = "blue", lwd = 3)
text(-1.5, 1, "y=x", col = "blue", cex = 2, adj = c(0, 0))
text(-1.5, 1.5, 
     paste0("r=", round(cor(y = mixed_out[rownames(fixed_out),]$logFC, x = fixed_out$logFC), 3)), 
     col = "red", cex = 2, adj = c(0, 0))
#dev.off()

fixed_out <- lvfit |> limma::eBayes() |> limma::topTable(number = Inf, coef = "typehot")
fixed_out$target <- rownames(fixed_out)
plot(x = fixed_out$logFC, y = -1*log10(fixed_out$adj.P.Val), pch = 16, col = "grey")

#pdf(file = "~/Desktop/equal-weighted_vs_mixed.pdf", width = 6, height = 6)
plot(y = mixed_out[rownames(fixed_out),]$logFC, x = fixed_out$logFC, 
     ylab = "mixed effects logFC", xlab = "fixed effects logFC", main = "Equal-Weighted Fixed Effects Estimates\nvs.\nMixed Effects Estimates", cex = 0.2)
abline(0, 1, col = "blue", lwd = 3)
text(-1, 1, "y=x", col = "blue", cex = 2, adj = c(0, 0))
text(-1, 1.5, 
     paste0("r=", round(cor(y = mixed_out[rownames(fixed_out),]$logFC, x = fixed_out$logFC), 3)), 
     col = "red", cex = 2, adj = c(0, 0))
#dev.off()

# At least in this case, the equal-weighting produces effect sizes that tend to 
# resemble those produced by mixed modeling to a greater extent. 

# However, this toy example still had 5 patients with data that spanned the two
# levels of type. What happens in a more unbalanced design in which only one 
# patient with data that spans? 

# Let's do the same exact investigation but with the prostate dataset.
meta <- read.csv(file = "~/Documents/projects/geomx/wu-gu-prostate-bladder/meta_amended.csv", row.names = 1)
cts <- read.csv(file = "~/Documents/projects/geomx/wu-gu-prostate-bladder/counts.csv", row.names = 1)

prostate <- meta[(meta$cancer_type == "prostate"),] |> rownames()
cts <- cts[,prostate]
meta <- meta[prostate,]

idx <- meta$sub_types %in% c("Ductal", "Acinar_crib") # Just one toy comparison
cts <- cts[,idx]
meta <- meta[idx,]

mm <- model.matrix(~sub_types_v2+patient_deid, data = meta)
contr <- (mm[meta$sub_types_v2 == "Ductal", ] |> colMeans()) - 
  (mm[meta$sub_types_v2 == "Acinar_crib", ] |> colMeans())
contr |> unname() |> round(4)
# [1]  0.0000  1.0000 -0.1429  0.1111  0.1111  0.1111 -0.1429 -0.1429  0.1111 -0.0317 -0.1429  0.1111 -0.1429  0.1111 -0.1429

norm <- edgeR::DGEList(counts = cts) |> edgeR::calcNormFactors(method = "upperquartile") |> limma::voom()
Y <- norm$E["TACSTD2", rownames(mm)]
d <- cbind("TROP2" = Y, meta)
mod <- lm(data = d, formula = TROP2 ~ sub_types_v2 + patient_deid)

EMM <- emmeans::emmeans(object = mod, specs = c("sub_types_v2"), weights = "cells")
(EMM@linfct[2,] - EMM@linfct[1,]) |> unname() |> round(4)
# [1]  0.0000  1.0000 -0.1429  0.1111  0.1111  0.1111 -0.1429 -0.1429  0.1111 -0.0317 -0.1429  0.1111 -0.1429  0.1111 -0.1429

pairs(EMM, reverse = T)
sum(contr * coef(mod))
multcomp::glht(model = mod, linfct = matrix(data = contr, ncol = length(contr))) |> summary(test = multcomp::univariate())

norm <- edgeR::DGEList(counts = cts[rownames(cts) != "NegProbe-WTX", rownames(meta)]) |> 
  edgeR::calcNormFactors(method = "upperquartile") |> 
  limma::voom(design = mm, plot = T) # Not great
lvfit <- limma::lmFit(object = norm, design = mm)
fixed_out <- lvfit |> limma::contrasts.fit(contrasts = contr) |> limma::eBayes() |> limma::topTable(number = Inf)
fixed_out$target <- rownames(fixed_out)
plot(x = fixed_out$logFC, y = -1*log10(fixed_out$adj.P.Val), pch = 16, col = "grey")

y <- edgeR::DGEList(counts = cts[rownames(cts) != "NegProbe-WTX", rownames(meta)]) |> 
  edgeR::calcNormFactors(method = "upperquartile")

param <- BiocParallel::SnowParam(8, "SOCK", progressbar = T)
v <- variancePartition::voomWithDreamWeights(counts = y,
                                             formula = ~sub_types_v2+(1|patient_deid),
                                             data = meta,
                                             plot = T, save.plot = T,
                                             BPPARAM = param)
attr(v, 'errors')
fitmm <- variancePartition::dream(exprObj = v,
                                  formula = ~sub_types_v2+(1|patient_deid),
                                  data = meta,
                                  BPPARAM = param
)
attr(fitmm, 'errors')

mixed_out <- variancePartition::eBayes(fit = fitmm) |> 
  variancePartition::topTable(coef = "sub_types_v2Ductal", number = Inf)
mixed_out$target <- rownames(mixed_out)

plot(x = mixed_out$logFC, y = -1*log10(mixed_out$adj.P.Val), pch = 16, col = "grey") # This literally shows nothing

#pdf(file = "~/Desktop/cell-weighted_vs_mixed_ex2.pdf", width = 6, height = 6)
plot(y = mixed_out[rownames(fixed_out),]$logFC, x = fixed_out$logFC, 
     ylab = "mixed effects logFC", xlab = "fixed effects logFC", main = "Cell-Weighted Fixed Effects Estimates\nvs.\nMixed Effects Estimates", cex = 0.2)
abline(0, 1, col = "blue", lwd = 3)
text(-1.5, 1, "y=x", col = "blue", cex = 2, adj = c(0, 0))
text(-1.5, 1.5, 
     paste0("r=", round(cor(y = mixed_out[rownames(fixed_out),]$logFC, x = fixed_out$logFC), 3)), 
     col = "red", cex = 2, adj = c(0, 0))
#dev.off()

fixed_out <- lvfit |> limma::eBayes() |> limma::topTable(number = Inf, coef = "sub_types_v2Ductal")
fixed_out$target <- rownames(fixed_out)
plot(x = fixed_out$logFC, y = -1*log10(fixed_out$adj.P.Val), pch = 16, col = "grey")

#pdf(file = "~/Desktop/equal-weighted_vs_mixed_ex2.pdf", width = 6, height = 6)
plot(y = mixed_out[rownames(fixed_out),]$logFC, x = fixed_out$logFC, 
     ylab = "mixed effects logFC", xlab = "fixed effects logFC", main = "Equal-Weighted Fixed Effects Estimates\nvs.\nMixed Effects Estimates", cex = 0.2)
abline(0, 1, col = "blue", lwd = 3)
text(-6, 1, "y=x", col = "blue", cex = 2, adj = c(0, 0))
text(-6, 1.5, 
     paste0("r=", round(cor(y = mixed_out[rownames(fixed_out),]$logFC, x = fixed_out$logFC), 3)), 
     col = "red", cex = 2, adj = c(0, 0))
#dev.off()

# And see the design for this? 
table(meta$patient_deid, meta$sub_types_v2)
#       Acinar_crib Ductal
# pt_10           0      2
# pt_11           1      0
# pt_14           0      1
# pt_15           0      1
# pt_23           0      1
# pt_26           1      0
# pt_28           1      0
# pt_29           0      1
# pt_30           1      1 <- ONLY ONE PATIENT SPANS!
# pt_31           1      0
# pt_33           0      1
# pt_34           1      0
# pt_35           0      1
# pt_4            1      0

# Furthermore, let's compare the model effects and the pt_30 differences.

#pdf(file = "~/Desktop/equal-weighted_vs_pt_30.pdf", width = 6, height = 6)
pt_30_diffs <- (norm$E[, meta$patient_deid == "pt_30"][,2] - norm$E[, meta$patient_deid == "pt_30"][,1])
model_effects <- fixed_out$logFC
plot(x = pt_30_diffs[rownames(fixed_out)], y = model_effects, 
     xlab = "pt_30 differences", ylab = "model effects", main = "pt_30 Differences\nvs.\nEqual-Weighted Fixed Effects Estimates")
abline(0, 1, col = "red", lty = 3, lwd = 3)
text(-6, 4, "y=x", col = "red", cex = 2, adj = c(0, 0))
#dev.off()

# They are literally the exact same values, which means that the entire result is 
# driven by one patient. 

# This is an extreme example, but the idea is evident. Using cell frequency 
# weighting protects against things like this happening when you are doing many
# comparisons and cannot diligently check the design for every single one. 

# Also note that if the design is balanced then the cell frequency weights will 
# be the same as the equal weights.

# I guess I see this method as having favorable trade-offs when the design is not 
# balanced and the sample size is small. Best case, you have a balanced design 
# and it is the same as using the more canonical equal weights. Worst case, it 
# is not balanced, but you approximate the ordinary marginal means. In either
# case, lack of spanning across the two levels you are comparing will not break
# down the analysis. 

# See here: 
# https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html#weights
# https://cran.r-project.org/web/packages/emmeans/vignettes/messy-data.html#weights
# https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html#exper
# https://www.tandfonline.com/doi/abs/10.1080/00031305.1980.10483031

