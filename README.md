## Bayesian Hierarchical Methods for Non-standard Compositional Data 

### PhD Thesis
### Catherine Holland

Code and Supplementary Material for PhD Thesis: *LINK*


### Abstract

Compositional data takes the form of parts of some whole with non-negative vectors that is
subject to a constraint. Compositional data can appear as proportions, percentages, general
non-negative values or counts. The inherent characteristics of compositional data, i.e. non-
negativity and the constraint to some total, pose unique challenges for traditional statistical
techniques. Compositional data arises across many real-world applications such as health,
environmental, forensic, financial and sports science. Further challenges occur when com-
positional data also includes other advanced data challenges such as multilevel hierarchical
structure, non-smooth time series or a spatial structure.

The main technique in the literature to overcome the complexities of compositional data
is to transform the components from the simplex (the sample space of compositional data)
into Euclidean space (the standard statistical space) using a log-ratio transformation. Once
transformed, standard statistical models can be applied. However, while this transformation
is powerful, it is not always suitable in practice. There are many features commonly found
within compositional data that prohibit log-ratio transformations. For example, when com-
positional data contains zeros, the log-ratios become undefined. Similarly, when the compon-
ents contain missing values, some or all of the log-ratio transformations may not produce
sensible results. Lastly, when compositional data consists of counts, applying a log-ratio
transformation may discard information on how the total count may impact the variance
and the possible values the counts can take. Thus, there is a need for frameworks that can
handle compositional data containing these features, as well as addressing advanced data
challenges.

This thesis presents novel Bayesian hierarchical frameworks designed to overcome the limit-
ations of log-ratio transformations in these instances. We apply and evaluate our proposed
frameworks to three applications of compositional data containing both a feature which
prevents log-ratio transformations and an advanced data challenge. These include: com-
positional data containing many zeros and a multilevel hierarchical structure, applied to
forensic elemental glass data; non-smooth time series containing a count structure and zero
values, applied to COVID-19 variant counts; and compositional data with a spatial pattern
containing zeros, applied to tree species proportions across a spatial grid. We assess the
performance of our frameworks through both in-sample and out-of-sample predictive exper-
iments, comparing with commonly used models. The results from the predictive experiments
demonstrate the effectiveness of our approaches, highlighting their contribution to compos-
itional data analysis and offering a robust alternative for handling real-world compositional
data.
