## Chapter 5: Methods for Spatial Compositional Data

Code and Supplementary Material.

We do not have permissions to share the data within this chapter so the data file is not uploaded here. 
I have however created a script to simulate spatial compositional proportions which can be used with the uploaded code to fit the framework. The outputted data is located in the `simulated_spatial_compositional_data.rds` file with the R script to produce this in `simulate_data.R`.


### Abstract

In this chapter, we propose an approach to modelling compositional data arranged over
a spatial domain, such that we need to account for both the compositional and spatial
structures. Furthermore, we target the additional challenges of accounting for both zero and
missing values in the spatial compositions.

Here, we propose a framework combining the Generalised-Dirichlet-Multinomial (GDM)
family of distributions with two-dimensional penalised regression splines that capture spatial
structure. We evaluate our approach through two posterior predictive experiments, one to
assess a novel variance parameter specification and another to assess how well the framework
can predict missing compositional counts.
