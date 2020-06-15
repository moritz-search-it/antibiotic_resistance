# Antibiotic Resistance genes
Mass screening of contigs for antimicrobial resistance or virulence genes with bioinformatic tool **abricate** by Seemann T. In total 4 different databases for ab resistance genes were screened.

* Abricate, Github <https://github.com/tseemann/abricate>
  + NCBI AMRFinderPlus - doi: 10.1128/AAC.00483-19
  + CARD - doi:10.1093/nar/gkw1004
  + Resfinder - doi:10.1093/jac/dks261
  + ARG-ANNOT - doi:10.1128/AAC.01310-13

Minimum DNA %identity 90, Minimum DNA %coverage 90

```{bash, }
conda activate abricate
abricate --minid [90] --mincov [90] --db ncbi --fofn fofn.txt >results_ncbi.tab
abricate --summary results_ncbi.tab > summary.txt
```


Manually screened in excel file, grouped by antibiotic classes and finally converted to csv and imported into R studio.

* Glycopeptide
* Tetracycline
* Aminoglycoside
* Aminocoumarin
* Diaminopyrimidine
* Fosfomycin
* Nitroimidazole
* Polypeptide antibiotic
* Sulfomanide
* Fluoroquinolone
* $\beta$-lactam antibiotic
* Broad-spectrum antimicrobial activity
* Protein-synthese inhibitor antibiotic
* multiple antibiotic resistances

## Database comparison
Database comparison showed different number of ab genes present 

![](https://github.com/moritz-search-it/antibiotic_resistance/blob/master/result_databases_faecis.png)

In _enterococcus faecis_, similar results are obtained from different databases. Manual comparasion revealed no new genes between CARD, ARGANNOT, NCBI and Resfinder. Only different names (_eat(A)_, _ant6_ in NCBI and _efm(A)_, _aad6_ in CARD)

![](https://github.com/moritz-search-it/antibiotic_resistance/blob/master/result_databases_coli.png)

Since CARD has also found more entries in e.coli, CARD resistance database is used for further exploration

CARD is a curated collection of characterized, peer-reviewed resistance genes and associated antibiotics and is updated monthly.

## Screening of antibiotic genes

![](https://github.com/moritz-search-it/antibiotic_resistance/blob/master/result_abscreening_coli.png)

![](https://github.com/moritz-search-it/antibiotic_resistance/blob/master/result_abscreening_faecis.png)

