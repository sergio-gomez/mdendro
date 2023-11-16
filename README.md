[![DOI](https://zenodo.org/badge/371634018.svg)](https://zenodo.org/badge/latestdoi/371634018)

# mdendro
Extended Agglomerative Hierarchical Clustering in R

## Description

[R](https://www.r-project.org) package [mdendro](https://github.com/sergio-gomez/mdendro) enables the calculation of **agglomerative hierarchical clustering** (AHC), extending the standard functionalities in several ways:

- Native handling of both **similarity** and **dissimilarity** (distances) matrices.

- Calculation of pair-group dendrograms and variable-group **multidendrograms** [[1](#references)].

- Implementation of the most common AHC methods in both **weighted** and **unweighted** forms: _single linkage_, _complete linkage_, _average linkage_ (UPGMA and WPGMA), _centroid_ (UPGMC and WPGMC), and _Ward_.

- Implementation of two additional parametric families of methods: **versatile linkage** [[2](#references)], and **beta flexible**. Versatile linkage leads naturally to the definition of two additional methods: _harmonic linkage_, and _geometric linkage_.

- Calculation of the **cophenetic** (or **ultrametric**) matrix.

- Calculation of five **descriptors** of the final dendrogram: _cophenetic correlation coefficient_, _space distortion ratio_, _agglomerative coefficient_, _chaining coefficient_, and _tree balance_.

- Calculation and plots of the descriptors for the parametric methods.

All this functionality is obtained with three functions: `linkage`, `descval` and `descplot`. Function `linkage` may be considered as a replacement for functions `hclust` (in package [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html)) and `agnes` (in package [cluster](https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/00Index.html)). To enhance usability and interoperability, the `linkage` class includes several methods for plotting, summarizing information, and class conversion.


## References

1. A. Fernández, S. Gómez. Solving non-uniqueness in agglomerative hierarchical clustering using multidendrograms. _Journal of Classification_ **25**, 43-65 (2008). DOI:[10.1007/s00357-008-9004-x](https:/doi.org/10.1007/s00357-008-9004-x).
2. A. Fernández, S. Gómez. Versatile linkage: A family of space-conserving strategies for agglomerative hierarchical clustering. _Journal of Classification_ **37**, 584-597 (2020). DOI:[10.1007/s00357-019-09339-z](https:/doi.org/10.1007/s00357-019-09339-z).


## Authors

- **Alberto Fernández**: Dept. Enginyeria Química, Universitat Rovira i Virgili, Tarragona (Spain). ([email](mailto:alberto.fernandez@urv.cat?subject=[mdendro])) ([ORCID](https://orcid.org/0000-0002-1241-1646)) ([Google Scholar](https://scholar.google.es/citations?user=AbH4r0IAAAAJ)) ([GitHub](https://github.com/albyfs))

- **Sergio Gómez**: Dept. Enginyeria Informàtica i Matemàtiques, Universitat Rovira i Virgili, Tarragona (Spain). ([web](https://webs-deim.urv.cat/~sergio.gomez/)) ([email](mailto:sergio.gomez@urv.cat?subject=[mdendro])) ([ORCID](http://orcid.org/0000-0003-1820-0062)) ([Google Scholar](https://scholar.google.es/citations?user=ETrjkSIAAAAJ)) ([GitHub](https://github.com/sergio-gomez)) ([Twitter](https://twitter.com/SergioGomezJ))


## Documentation

The full documentation of [mdendro](https://github.com/sergio-gomez/mdendro), including description, installation, tutorial, rationale and reference manual, can be found [here](https://sergio-gomez.github.io/mdendro/).
