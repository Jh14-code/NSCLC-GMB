# NSCLC-GMB (R Code of GMB VS TMB)
This repository contains the code for the paper "Gene Mutational Burden as A Potential Better Predictor of Immunotherapy Benefits than Tumor Mutational Burden in Patients With NSCLC" by Jing Huang et al.
# The input data source
All clinical cohort and patient data were obtained from the cBioPortal database (https://www.cbioportal.org) and previous studies[1-9]. The discovery and validation cohort included four ICI-treated cohorts (1661 patients received ICB therapy from the cohort of Samstein et al.[1], 35 NSCLC patients treated with anti-PD-1 from the cohort of Rizvi et al. (2015)[2], 240 NSCLC patients treated with anti-PD-(L)1 from the cohort of Rizvi et al. (2018)[3], and 75 NSCLC patients treated with anti-PD-1 plus anticytotoxic T-lymphocyte antigen 4 (anti-CTLA-4) from the cohort of Hellmann et al.[4]). We also collected three non-ICI-treated cohorts (566 LUAD patients from TCGA[5], 487 LUSC patients from TCGA[5], and 10336 advanced tumor patients from zehir et al.[6]), and the blood-based data from POPLAR and OAK trials[7-9].
# The content of the code
All the code were performed by R software (version 4.2.0). The repository is organized as follows:
 1. Gene lists and the predictive value of GMB in the discovery cohort.R
 2. The predictive value of GMB in the validation cohorts.R
 3. GMB was not a prognostic factor in the non-ICB-treated cohort.R
 4. GMB was also useful in multiple cancers.R
 5. Blood-based GMB (bGMB) was a potential alternative to blood-based TMB (bTMB) in predicting the clinical benefits of ICB treatment in NSCLC.R
 6. Potential mechanisms associated with GMB in predicting the ICB treatment efficacy.R
# Detailed description of the code
1. Gene lists and the predictive value of GMB in the discovery cohort.R

 1.1 select candidate genes which are associated with OS after ICB therapy.

 1.2 The predictive value of GMB in the discovery cohort. We explore the association between GMB and ICB treatment outcomes: OS and PFS. Then we compare the predictive power of GMB VS. TMB by kaplan-Meier curves.

2. The predictive value of GMB in the validation cohorts.R

 2.1 Validation of the predictive value of GMB in the Rizvi cohort(2018).

 2.2 Validation of the predictive value of GMB in the Hellmann cohort.

 2.3 Validation of the predictive value of GMB in the Rizvi cohort (2015).

 2.4 The forest of the predictive value of GMB in validation cohort.

 2.5 Validation of the predictive value of GMB in the validation cohort.

 2.6 Compare of the C index of GMB,TMB and other single mutatin gene.

3. GMB was not a prognostic factor in the non-ICB-treated cohort.R

 3.1 GMB was not a prognostic factor in the Zehir Cohort.

 3.2 GMB was not a prognostic factor in the TCGA cohort.

4. GMB was also useful in multiple cancers.R

 4.1 GMB was also applicable to multiple cancers with proper gene mutation sets.

 4.2 GMB was better than risk model in NSCLC.

5. Blood-based GMB (bGMB) was a potential alternative to blood-based TMB (bTMB) in predicting the clinical benefits of ICB treatment in NSCLC.R

 5.1 Blood-based GMB (bGMB) predicting the clinical benefits of ICB treatment in NSCLC.

 5.2 Blood-based TMB (bTMB) predicting the clinical benefits of ICB treatment in NSCLC.
 
 5.3 clinical factors associate with OS in the POPLAR and OAK cohort.

6. Potential mechanisms associated with GMB in predicting the ICB treatment efficacy.R

 6.1 Neoantigen in GMB-H vs. GMB-L.

 6.2 Pathway analysis by GSEA.
 
 6.3 tumor immune microenvironment in GMB-H vs. GMB-L

# References

1.Samstein RM, Lee CH, Shoushtari AN, et al. Tumor mutational load predicts survival after immunotherapy across multiple cancer types. Nature genetics 2019;51:202-206.

2.Rizvi NA, Hellmann MD, Snyder A, et al. Cancer immunology. Mutational landscape determines sensitivity to PD-1 blockade in non-small cell lung cancer. Science (New York, NY) 2015;348:124-128.

3.Rizvi H, Sanchez-Vega F, La K, et al. Molecular Determinants of Response to Anti-Programmed Cell Death (PD)-1 and Anti-Programmed Death-Ligand 1 (PD-L1) Blockade in Patients With Non-Small-Cell Lung Cancer Profiled With Targeted Next-Generation Sequencing. Journal of clinical oncology : official journal of the American Society of Clinical Oncology 2018;36:633-641.

4.Hellmann MD, Nathanson T, Rizvi H, et al. Genomic Features of Response to Combination Immunotherapy in Patients with Advanced Non-Small-Cell Lung Cancer. Cancer cell 2018;33:843-852.e844.

5.Liu J, Lichtenberg T, Hoadley KA, et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell 2018;173:400-416.e411.

6.Zehir A, Benayed R, Shah RH, et al. Mutational landscape of metastatic cancer revealed from prospective clinical sequencing of 10,000 patients. Nat Med 2017;23:703-713.

7.Fehrenbacher L, Spira A, Ballinger M, et al. Atezolizumab versus docetaxel for patients with previously treated non-small-cell lung cancer (POPLAR): a multicentre, open-label, phase 2 randomised controlled trial. Lancet (London, England) 2016;387:1837-1846.

8.Rittmeyer A, Barlesi F, Waterkamp D, et al. Atezolizumab versus docetaxel in patients with previously treated non-small-cell lung cancer (OAK): a phase 3, open-label, multicentre randomised controlled trial. Lancet (London, England) 2017;389:255-265.

9.Gandara DR, Paul SM, Kowanetz M, et al. Blood-based tumor mutational burden as a predictor of clinical benefit in non-small-cell lung cancer patients treated with atezolizumab. Nat Med 2018;24:1441-1448.
