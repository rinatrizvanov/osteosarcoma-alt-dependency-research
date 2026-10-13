# Axis G — Expression vs structural proxies: when expression classifiers miss SV-driven phenotypes

## Note on this axis

The axis as literally posed has **no dedicated corpus**. The pre-registered query below
returned only 6 works, and on inspection none were on topic (PDX
analysis workflows, breast organoid biobanks, DCIS multi-omics). Four further broadenings
were tried and drifted into maize, rice, yeast and Parkinson's genomics.

This is itself a finding: "when does an expression proxy fail to capture an SV-driven
phenotype" is not a named research area with its own literature. The nearest real corpus
inverts the question — how much of gene expression variance structural and copy-number
variants actually explain — which bounds how well expression can ever proxy SV status.
Both queries are recorded.

## Queries

Pre-registered (returned nothing on topic, retained for the record):
```
https://api.openalex.org/works?filter=title_and_abstract.search%3A%28%22gene+expression+signature%22+OR+%22transcriptomic+classifier%22+OR+%22expression-based+classifier%22%29+AND+%28%22structural+variant%22+OR+%22copy+number+alteration%22%29+AND+%28%22surrogate%22+OR+%22proxy%22+OR+%22concordance%22+OR+%22misclassification%22%29&sort=cited_by_count%3Adesc&per_page=15&mailto=<CONTACT_EMAIL>&api_key=<API_KEY>
```

Reframed, sorted by citations:
```
https://api.openalex.org/works?filter=title_and_abstract.search%3A%28%22structural+variant%22+OR+%22structural+variation%22+OR+%22copy+number+alteration%22%29+AND+%28%22gene+expression%22%29+AND+%28%22cancer%22+OR+%22tumour%22+OR+%22tumor%22%29+AND+%28%22impact%22+OR+%22consequence%22+OR+%22cis%22+OR+%22dosage%22%29&sort=cited_by_count%3Adesc&per_page=15&mailto=<CONTACT_EMAIL>&api_key=<API_KEY>
```

Reframed, sorted by recency:
```
https://api.openalex.org/works?filter=title_and_abstract.search%3A%28%22structural+variant%22+OR+%22structural+variation%22+OR+%22copy+number+alteration%22%29+AND+%28%22gene+expression%22%29+AND+%28%22cancer%22+OR+%22tumour%22+OR+%22tumor%22%29+AND+%28%22impact%22+OR+%22consequence%22+OR+%22cis%22+OR+%22dosage%22%29%2Cpublication_year%3A%3E2023&sort=publication_date%3Adesc&per_page=8&mailto=<CONTACT_EMAIL>&api_key=<API_KEY>
```

Total matching works (reframed citation query): 548

---

## Top by citation count

https://doi.org/10.1038/s41576-018-0007-0
Structural variation in the 3D genome — Malte Spielmann (2018, Nature Reviews Genetics; 737 citations)

https://doi.org/10.1136/gutjnl-2011-301839
A comprehensive survey of genomic alterations in gastric cancer reveals systematic patterns of molecular exclusivity and co-occurrence among distinct therapeutic targets — Niantao Deng (2012, Gut; 645 citations)

https://doi.org/10.1101/gr.145144.112
The transcriptional landscape and mutational profile of lung adenocarcinoma — Jeong‐Sun Seo (2012, Genome Research; 585 citations)

no DOI available — OpenAlex W2135185276
Impact of DNA amplification on gene expression patterns in breast cancer. — Elizabeth Hyman (2002, PubMed; 494 citations)

https://doi.org/10.1038/ng.3722
Pan-cancer analysis of somatic copy-number alterations implicates IRS4 and IGF2 in enhancer hijacking — Joachim Weischenfeldt (2016, Nature Genetics; 406 citations)

https://doi.org/10.1038/s41586-020-1970-0
Genomic basis for RNA alterations in cancer — PCAWG Transcriptome Core Group (2020, Nature; 393 citations)

https://doi.org/10.1158/1078-0432.ccr-19-3231
Single-Cell Genomic Characterization Reveals the Cellular Reprogramming of the Gastric Tumor Microenvironment — Anuja Sathe (2020, Clinical Cancer Research; 390 citations)

https://doi.org/10.1038/ng.3173
Gene expression analysis identifies global gene dosage sensitivity in cancer — Rudolf S.N. Fehrmann (2015, Nature Genetics; 349 citations)

https://doi.org/10.1016/j.cell.2012.12.034
Integrative eQTL-Based Analyses Reveal the Biology of Breast Cancer Risk Loci — Qiyuan Li (2013, Cell; 308 citations)

https://doi.org/10.1016/j.ajpath.2017.04.022
Breast Cancer Molecular Stratification — Hege Giercksky Russnes (2017, American Journal Of Pathology; 280 citations)

https://doi.org/10.3390/cancers12030764
Global DNA Hypomethylation in Epithelial Ovarian Cancer: Passive Demethylation and Association with Genomic Instability — Wa Zhang (2020, Cancers; 278 citations)

https://doi.org/10.1038/ng.3849
Genomic consequences of aberrant DNA repair mechanisms stratify ovarian cancer histotypes — Yi Kan Wang (2017, Nature Genetics; 274 citations)

https://doi.org/10.1038/ng.3073
An integrated genomics approach identifies drivers of proliferation in luminal-subtype human breast cancer — Michael L. Gatza (2014, Nature Genetics; 256 citations)

https://doi.org/10.1016/j.semcancer.2010.03.001
A mobile threat to genome stability: The impact of non-LTR retrotransposons upon the human genome — Miriam K. Konkel (2010, Seminars in Cancer Biology; 191 citations)

https://doi.org/10.1002/path.2424
Genome‐wide analysis of DNA copy number alterations and gene expression in gastric cancer — Yoshiyuki Tsukamoto (2008, The Journal of Pathology; 130 citations)

---

## Most recent (post-2024, deduplicated against the list above)

https://doi.org/10.1200/jco-26-00742
Prognostic Impact of Telomere Maintenance Mechanisms in Neuroblastoma — Carolina Rosswog (2026, Journal of Clinical Oncology; 0 citations)

https://doi.org/10.3389/fbinf.2026.1811813
The impact of deep learning and omics data in transforming precision therapy for brain cancer — Sana Munquad (2026, Frontiers in Bioinformatics; 0 citations)

https://doi.org/10.82444/warw.33820447
Particle Monte Carlo methods for the integrative cluster analysis of multiple genomic datasets — Nathan J. Cunningham (2026, Warwick Research Archive Portal (University of Warwick); 0 citations)

https://doi.org/10.1007/s10142-026-02012-w
Exploring the transcriptional crosstalk between adipose tissue and locoregional recurrence in breast cancer using independent component analysis — Marlous Arjaans (2026, Functional & Integrative Genomics; 0 citations)

https://doi.org/10.1038/s41467-026-77331-4
Persistent chromatin alterations and gene expression reprogramming follow widespread DNA damage in glioblastoma — Aram Sandaldjian Modrek (2026, Nature Communications; 0 citations)

