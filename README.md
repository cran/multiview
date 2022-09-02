
# Cooperative Learning for Multi-view Analysis <img src="man/figures/logo.png" width="100" align="right" />

The package `multiview` is a new method for supervised learning with
multiple sets of features called *views*. The multi-view problem is
especially important in biology and medicine, where “-omics” data such
as genomics, proteomics and radiomics are measured on a common set of
samples. Cooperative learning combines the usual squared error loss of
predictions with an “agreement” penalty to encourage the predictions
from different data views to agree. By varying the weight of the
agreement penalty, we get a continuum of solutions that include the
well-known early and late fusion approaches. Cooperative learning
chooses the degree of agreement (or fusion) in an adaptive manner, using
a validation set or cross-validation to estimate test set prediction
error.

In the setting of cooperative regularized linear regression, the method
combines the lasso penalty with the agreement penalty, yielding feature
sparsity. The method can be especially powerful when the different data
views share some underlying relationship in their signals that can be
exploited to boost the signals.

As shown in Ding et al. ([2021](#ref-cooperative)), cooperative learning
achieves higher predictive accuracy on simulated data and real
multiomics examples of labor onset prediction and breast ductal
carcinoma in situ and invasive breast cancer classification. Leveraging
aligned signals and allowing flexible fitting mechanisms for different
modalities, cooperative learning offers a powerful approach to
multiomics data fusion.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-cooperative" class="csl-entry">

Ding, Daisy Yi, Shuangning Li, Balasubramanian Narasimhan, and Robert
Tibshirani. 2021. “Cooperative Learning for Multi-View Analysis.” *arXiv
Preprint arXiv:2112.12337*.

</div>

</div>
