layout: page
title: "Stats"
permalink: /Stats

## <a name='TOC'>Stats</a>

1. [Normal distribution testing](#Normal)
2. [Differential abundance testing](#Differential)

**<a name="Normal">Test if a distribution is normal</a>**

Run a shapiro-wilk test for normality [source](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.shapiro.html).

The example uses a pandas dataframe with a column 'frequency'. If pvalue is significant, data is not normal
```python
from scipy import stats
stats.shapiro(df['frequency'])

```

**<a name="Differential">Differential abundance testing</a>**

Some background: A highly cited [paper](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y) from Rob Knight's group on differential abundace testing. Basically, they discuss fancier methods to test compositional relative abundance changes
>"For OTU differential abundance testing between groups (e.g., case vs. control), a common approach is to first rarify the count matrix to a fixed depth and then apply a nonparametric test (e.g., the Mann-Whitney/Wilcoxon rank-sum test for tests of two groups; the Kruskal-Wallis test for tests of multiple groups). Nonparametric tests are often preferred because OTU counts are not exactly normally distributed [40]. However, when analyzing relative abundance data, this approach does not account for the fact that the relative abundances are compositional".

Compositional, in other words, means groups/categories (e.g., Family, Order...), as opposed to just raw counts.

Regardless, here are two examples:

Kruskal wallis
```python
stats.kruskal(df['frequency1'], df['frequency2'])
```

Wilcoxon rank-sum/Mann Whitney
```python
stats.wilcoxon(df['frequency1'], df['frequency2'])
```
