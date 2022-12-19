Normal distribution testing. If pvalue is significant, data is not normal
```
res140 = df_resolution3.query('Hidef_resolution_maximum == 140')
stats.shapiro(res140['frequency'])

```

A highly cited [paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y) from Rob Knight's group on differential abundace testing. Basically, they discuss fancier methods to test compositional relative abundance changes
>"For OTU differential abundance testing between groups (e.g., case vs. control), a common approach is to first rarify the count matrix to a fixed depth and then apply a nonparametric test (e.g., the Mann-Whitney/Wilcoxon rank-sum test for tests of two groups; the Kruskal-Wallis test for tests of multiple groups). Nonparametric tests are often preferred because OTU counts are not exactly normally distributed [40]. However, when analyzing relative abundance data, this approach does not account for the fact that the relative abundances are compositional". In other words, groups/categories, as opposed to just raw counts.

Kruskal wallis
```
stats.kruskal(df['frequency'], df['frequency'])
```

Wilcoxon rank-sum/Mann Whitney
```
stats.wilcoxon(df['frequency'], df['frequency'])
```
