# Methodological note: sensitivity of the PCA to the `NA` source column

The published PCA (the original version v5) was reproduced by retaining all nonconstant percentage columns, **including `NA` when nonconstant**. Independent inspection of the source XLSX files confirmed the presence of this class in DOF, SSF and DSF. This class is labeled simply `NA` in the supplied source sheets and its underlying MapBiomas classification meaning has not been independently verified for this dataset.

With the original class set, the two-PC cumulative variance is: Basin 68.4%; DOF 66.0%; SSF 74.7%; DSF 58.6%.

An independent *sensitivity calculation* excluding only the `NA` class produced approximate two-PC cumulative explained variances: Basin 68.4%; DOF 69.5%; SSF 78.9%; DSF 63.7%. These are **not** published estimates. The changes indicate that the `NA` variable is not numerically immaterial.

Before publication, the authors should decide explicitly, based on the MapBiomas class mapping, whether `NA` belongs in the environmental explanatory variables. If it is excluded, rerun the full PCA and revise Figure 5, Table 2, Results, and supporting methods. If retaining it, disclose what `NA` represents and why it was considered an analytical class. The main reproducibility script retains it to guarantee agreement with the submitted numerical results and to avoid silently changing the scientific analysis.

The optional script `scripts/03_optional_pca_sensitivity_NA.R` enables a transparent alternative calculation without modifying the archived reference values.
