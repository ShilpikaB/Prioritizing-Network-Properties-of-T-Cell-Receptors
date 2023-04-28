# 1 Introduction
Tcells are crucial components of the adaptive immune system, mediating anti-tumoral immunity and immune response to infections. They are necessary for effective host-response to a wide range of pathogens. Tcells are defned by their Tcell receptors (TCRs), which are protein complexes on the Tcell surface. TCRs mount a response to harmful foreign invaders by targeting specifc antigens based on nucleotide sequence. TCRs act as the arms of the Tcells with memory and can remember harmful pathogens they have seen before, thereby providing a life-long protection which enables a swift response in case of a similar future encounter. Thus, understanding the TCR repertoire could lead to insights regarding immune response pathology while also discovering indicative bio-markers and lead to therapeutic strategies.

### 1.1 Motivation & Objective
The TCR repertoire network data of patients collected as a part of a drug trail showed significant difference in the network structures. It was observed that patients exhibiting focused and larger clusters of the TCR repertoire attained significantly longer overall survival (OS) than those with smaller clusters (Naidus et al., 2021). Therefore, drawing quantitative analysis on the TCR repertoire network properties in the ‘longer overall survival’ and the ‘shorter overall survival’ cohorts may provide a better understanding of the immune landscape involving Tcell response.

<div align="center">
  <img src="https://user-images.githubusercontent.com/82466266/235205206-aaf85ad8-29a0-45e3-9722-0940c51b2243.JPG" width=40% height=30%>
</div>

### 1.2 Challenges
___Data Heterogeneity___
- Some of the TCR repertoire network properties have single values whereas, other properties are of varying lengths. This makes the TCR network data heterogeneous and unsuitable for direct statistical inferences across any two patients.
- Since the TCR repertoire is constantly adapting to the health and the environmental factors of the patient, the network properties are continually shaping. Given any two patients the TCR repertoire is never the same. 
- Less than 20% overlap is observed in the TCR repertoires for the same subject. 
- Heterogeneity also complicates the data simulation process required to perform the simulation study.

### 1.3 Contribution
- Strategy to extract features from heterogeneous TCR repertoire network data.
- A novel statistical method to prioritize the network properties - Group Lasso with Permutation Assisted Tuning.
- Mimicked real property distributions and correlation structure to simulate network properties.
- Demonstrated proposed methods and schemes using simulation study.

# Data Preprocessing
 The overall survival months for each of the enrolled patients was also provided. The median overall survival months is 20.3. For this analysis, the network properties are used as the explanatory variables and the overall survival month (OS mon) as the response variable. Patients with OS mon ≥ 20.3 are categorized into ‘longer overall survival’ group and patients with OS mon < 20.3 are categorized into ‘shorter overall survival’ group.
# Methods & Implementation
<div align="center">
  <img src="" width=80% height=50%>
</div>

# Code (in R)
R code: https://github.com/ShilpikaB/Prioritizing-Network-Properties-of-T-Cell-Receptors/blob/main/R_script.R


# Results
- 


# Discussion



# References
- Elliot Naidus, Jerome Bouquet, David Y. Oh, Timothy J. Looney, Hai Yang, Lawrence Fong, Nathan E. Standifer, Li Zhang. “Early changes in the circulating T cells are associated with clinical outcomes after PD-L1 blockade by durvalumab in advanced NSCLC patients”. In: Cancer Immunology, Immunotherapy 70:2095–2102 (2021).
- Enkelejda Miho, Rok Roskar, Victor Greif and Sai T.Reddy. “Large-scale network analysis reveals the sequence space architecture of antibody repertoires”. In: Nature Communications 10:1321 (2019).
- Ming Yuan and Yi Lin. “Model selection and estimation in regression with grouped variables”. In: Journal of the Royal Statistical Society. Series B 68.Part 1 (2006), pp. 49–67.
- Songshan Yang, Jiawei Wen, Scott T. Eckert, Yaqun Wang, Dajiang J. Liu, Rongling Wu, Runze Li1 and Xiang Zhan. “Prioritizing genetic variants in GWAS with lasso using permutation-assisted tuning”. In: Bioinformatics 36:3811-7 (2020).
- Robert Tibshirani. “Regression Shrinkage and Selection via the Lasso”. In: Journal of the Royal Statistical Society. Series B (Methodological) 58 (1996), pp. 267–288.
- Yang Zhou, Rong Jin, Steven C. H. Hoi. “Exclusive Lasso for Multi-task Feature Selection”. In: JMLR Workshop and Conference Proceedings: 13th International Conference on Artifcial Intelligence and Statistics (AISTATS) 9 (2010), pp. 988–995.
