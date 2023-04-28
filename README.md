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
- A novel statistical method to prioritize the network properties - _Group Lasso with Permutation Assisted Tuning_.
- Mimicked real property distributions and correlation structure to simulate network properties.
- Demonstrated proposed methods and schemes using simulation study.

# 2 Methods & Implementation
___2.1 Method outline___
<div align="center">
  <img src="https://user-images.githubusercontent.com/82466266/235205807-1a306922-ae92-472f-945f-87f0027466c8.JPG" width=60% height=40%>
  <img src="https://user-images.githubusercontent.com/82466266/235209064-86c9d408-8dc9-420b-aecb-c4418eb45729.JPG" width=60% height=40%>
 </div>
 
___2.2 Implementation___: So far Permutation Assisted Tuning (PAT) has been used only on Lasso models and is known to have lower false positives than Cross-validation technique. We add novelty by extending the application of PAT to Group Lasso model and formulating the required setup.
 
 <div align="center">
  <img src="https://user-images.githubusercontent.com/82466266/235208784-2a8ce7d5-a2f0-424c-8488-36f8b8d4711b.JPG" width=50% height=50%>
  <img src="https://user-images.githubusercontent.com/82466266/235208813-a25b034a-a18d-4ddc-951f-44b507a96c9b.JPG" width=40% height=40%>
</div>

# 3 Code (in R)
R code: https://github.com/ShilpikaB/Prioritizing-Network-Properties-of-T-Cell-Receptors/blob/main/R_script.R


# 4 Results
### 4.1 Real data analysis
 <div align="center">
  <img src="https://user-images.githubusercontent.com/82466266/235210760-c38026db-a2d4-4e8e-9ae3-23d455de5707.JPG" width=60% height=50%>
</div>

### 4.2 Simulated data analysis
<div align="center">
  <img src="https://user-images.githubusercontent.com/82466266/235210922-7cfd7366-cd16-47fb-b215-8e03076675fb.JPG" width=60% height=50%>
</div>

# 5 Discussion
This study provides evidence that Permutation Assisted Tuning performs better than Cross-validation technique for both grouped and ungrouped feature selection scenarios.

# 6 References
- Elliot Naidus, Jerome Bouquet, David Y. Oh, Timothy J. Looney, Hai Yang, Lawrence Fong, Nathan E. Standifer, Li Zhang. “Early changes in the circulating T cells are associated with clinical outcomes after PD-L1 blockade by durvalumab in advanced NSCLC patients”. In: Cancer Immunology, Immunotherapy 70:2095–2102 (2021).
- Enkelejda Miho, Rok Roskar, Victor Greif and Sai T.Reddy. “Large-scale network analysis reveals the sequence space architecture of antibody repertoires”. In: Nature Communications 10:1321 (2019).
- Ming Yuan and Yi Lin. “Model selection and estimation in regression with grouped variables”. In: Journal of the Royal Statistical Society. Series B 68.Part 1 (2006), pp. 49–67.
- Songshan Yang, Jiawei Wen, Scott T. Eckert, Yaqun Wang, Dajiang J. Liu, Rongling Wu, Runze Li1 and Xiang Zhan. “Prioritizing genetic variants in GWAS with lasso using permutation-assisted tuning”. In: Bioinformatics 36:3811-7 (2020).
- Robert Tibshirani. “Regression Shrinkage and Selection via the Lasso”. In: Journal of the Royal Statistical Society. Series B (Methodological) 58 (1996), pp. 267–288.
- Yang Zhou, Rong Jin, Steven C. H. Hoi. “Exclusive Lasso for Multi-task Feature Selection”. In: JMLR Workshop and Conference Proceedings: 13th International Conference on Artifcial Intelligence and Statistics (AISTATS) 9 (2010), pp. 988–995.

# 7 Acknowledgements
<div align="center">
  <img src="https://user-images.githubusercontent.com/82466266/235212142-6c1b55b1-8328-4c92-885e-82d5a60e735f.JPG" width=60% height=50%>
</div>
