# Project Information

Analysis completed on behalf of the MGH KFCCR Tumor Cartography Core for the Wu Lab. 
- All steps were conducted in a `renv` environment that can be restored from the `renv.lock` file in this repo. 
- All raw and/or processed data can be requested by emailing the corresponding author.

## Notes About DE Analyses

One may notice in the DE scripts that the way I construct the contrasts is atypical.

Let's use the following comparison in the prostate project as an example: Ductal vs. Acinar_crib

We could construct the contrast with the following code:

```         
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr_mat <- limma::makeContrasts(
  DvsA = sub_types_v2Ductal,
  levels = mm
)
```

This leads to this contrast vector, which we will denote with X:

```         
         Intercept sub_types_v2Ductal  patient_deidpt_11  patient_deidpt_14  patient_deidpt_15  patient_deidpt_23 
                 0                  1                  0                  0                  0                  0 
 patient_deidpt_26  patient_deidpt_28  patient_deidpt_29  patient_deidpt_30  patient_deidpt_31  patient_deidpt_33 
                 0                  0                  0                  0                  0                  0 
 patient_deidpt_34  patient_deidpt_35   patient_deidpt_4 
                 0                  0                  0 
```

Alternatively, we could do this for contrast construction:

```         
mm <- model.matrix(~sub_types_v2+patient_deid, data = dds@colData)
contr <- (mm[dds@colData$sub_types == "Ductal", ] |> colMeans()) - 
  (mm[dds@colData$sub_types == "Acinar_crib", ] |> colMeans())
```

This leads to this contrast vector, which we will denote with Y:

```         
       (Intercept) sub_types_v2Ductal  patient_deidpt_11  patient_deidpt_14  patient_deidpt_15  patient_deidpt_23 
        0.00000000         1.00000000        -0.14285714         0.11111111         0.11111111         0.11111111 
 patient_deidpt_26  patient_deidpt_28  patient_deidpt_29  patient_deidpt_30  patient_deidpt_31  patient_deidpt_33 
       -0.14285714        -0.14285714         0.11111111        -0.03174603        -0.14285714         0.11111111 
 patient_deidpt_34  patient_deidpt_35   patient_deidpt_4 
       -0.14285714         0.11111111        -0.14285714
```

The design for this comparison is this:

```         
        Acinar_crib Ductal
  pt_10           0      2
  pt_11           1      0
  pt_14           0      1
  pt_15           0      1
  pt_23           0      1
  pt_26           1      0
  pt_28           1      0
  pt_29           0      1
  pt_30           1      1
  pt_31           1      0
  pt_33           0      1
  pt_34           1      0
  pt_35           0      1
  pt_4            1      0
```

One can see that in X, the patients are all weighted equally. In Y, the patients have weights.

The former way of testing (contrast X) is the more typical method and is the method shown in nearly all DE analysis vignettes. However, we had a difficult study design in this project (see table above, for example).

Simply providing the formula `~sub_types_v2+patient_deid` and then testing the coefficient of `sub_types_v2` is a "paired" test that is attempting to estimate and evaluate the within-patient effect size of the `sub_types_v2` variable. In other words, it would estimate the subtype effect size, while conditioning on the patient.

However, this is not exactly what we wanted in this study. We wanted to test the between-patient effect size, while accounting for the fact that some patients may have been sampled multiple times (i.e. pseudo-replicated).

Essentially, I felt that it made sense to try to account for the possible correlation between samples from the same patient, while also allowing these patients to contribute more to the comparison at hand, since they provided more data.

In practice, I tested **weighted** contrasts after including the patient identifier in the model. The modeling "regresses out" patient-specific effects. However, with the specific contrasts that I computed, if a patient has more samples in one direction of the comparison, it is then able to contribute more to the logFC estimate than the other patients. As a result, this method allows the final logFC to better represent the actual unbalanced study design at hand, instead of testing the contrasts under the assumption that the study design was balanced. **Really, we estimated a marginal subtype effect that is sample composition weighted.** I think it is important to note this distinction.

To provide one example, the PIGR gene has a huge logFC value for contrast X. This seems to be because the within-patient expression difference is huge for pt_30. Contrast Y allows the other patients to contribute more to the comparison, which can be protective, since we obtain an effect size that is more consistent with cross-patient patterns. We could see from the actual data points that PIGR really is not DE in any striking way across all patients.

Two questions may arise:

**1) Why not pseudo-bulk?**

This could be a good idea. However, notice in our design that some patients might have both levels of the variable of interest. Thus, the situation does not always not become much simpler after pseudo-bulking. Furthermore, for this study, the postdoc wanted to retain the individual AOIs for other analyses, like ssGSEA. This was a bit of a Catch-22, to be honest. I opted to keep the individual AOIs, account for pseudo-replication as best I could, and represent the unbalanced structure of the data in the final contrasts being tested.

**2) Why not use a mixed model?**

Honestly, a mixed model with a random intercept term for patient is probably better. However, this can be computationally expensive. At the time of the analysis, I did not know much about mixed models, so I opted for a familiar tool. However, I would probably opt for `variancePartition`'s `dream` function now.

**Conclusions**

I think that we can debate all day about the "correct" way to construct the contrasts for testing. The whole point of this analysis though was to identify genes that have differing expression patterns across the groups of interest. Looking at the heatmaps, I think that we were able to do that; one can see clear differences. Overall, my outlook is this: we tried hard to be statistically rigorous in our methods and thought through how we were representing the data, which allowed us to identify some real differences, evidenced by visualization. Our modeling and testing may not be perfect. The exact steps taken may be improved by using certain alternative modeling strategies. However, we nonetheless captured some reasonable differences in the data that are biologically interesting and supportable.
