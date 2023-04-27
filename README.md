# 1 Introduction
T cells are crucial components of the adaptive immune system, mediating anti-tumoral immunity and immune response to infections. They are necessary for efective host-response to a wide range of pathogens. T cells are defned by their T cell receptors (TCRs), which are protein complexes on the T-cell surface. TCRs mount a response to harmful foreign invaders by targeting specifc antigens based on nucleotide sequence. TCRs act as the arms of the T cells with memory and can remember harmful pathogens they have seen before, thereby providing a life-long protection which enables a swift response in case of a similar future encounter. Thus, understanding the TCR repertoire could lead to insights regarding immune response pathology while also discovering indicative bio-markers and lead to therapeutic strategies.
As the immune repertoire ages, it is shaped based on the environmental exposure of an individual throughout their lifetime. Therefore, performing statistical inferences directly on TCR data between subjects is challenging due to its heterogeneous nature. In fact, there is less than 20% overlap across repertoire, even for the same subject. However, it is observed that the similarity among TCRs sequence directly infuences the antigen recognition breadth. Therefore, interrogation of TCR sequence similarity can add an important layer of information. This can be achieved through network property analysis of TCR repertoire. A clonal network is constructed where each clone is defned as a node, and then based on the sequence distance (Levenshtein distance), an edge is drawn based on a certain similarity condition (e.g., one letter diference in sequence).

### 1.1 Background
In this work we analysis the data of 65 patients from the Phase I trial (NCT01693562, 14 September 2012) of durvalumab, an immune checkpoint inhibitor (ICI) designed to activate exhausted tumor-reactive T cells. Durvalumab consolidation therapy is administered to patients with stage III, non-small cell lung cancer (NSCLC) and their immunophenotypic responses are observed. The patients exhibiting increased TCR repertoire diversity on day 15 attained signifcantly longer overall survival (OS) than those with decreased diversity. Patients with larger TCR clusters showed improved OS than patients with smaller TCR clusters ([1]Elliot Naidus et al., 2021). It was inferred that early TCR repertoire diversifcation after durvalumab therapy for NSCLC may be predictive of increased survival. Therefore, drawing quantitative analysis of the TCR repertoire in ‘longer overall survival’ and ‘shorter overall survival’ cohorts may provide a better understanding of the immune landscape involving T cell response. This information can then be used to develop tools to improve patient stratifcation, prediction of disease outcome, and patient response to treatments.

### 1.2 Motivation & Objective
The immunophenotypic response data captures the TCR repertoire details for the 65 patients. The TCR network is continually shaping over a patient’s lifetime and is also impacted as a response to the immunotherapy administered to the patient, thereby making the data heterogeneous in nature. The network data captured for our analysis is shown in the Figure 1.1. A total of ffteen network and non-network properties of the TCR repertoire data are used for this work which are referenced as the TCR network properties collectively. Some of these network properties are global and some are clonal (local) network properties ([2]Miho et al., 2019). Another data set representing the overall survival stats for these 65 patients was also referenced. It was also made available that the patients with overall survival months (OS mon) ≥ 20.3 have a higher survival chance than the other patients. For the analysis, the network properties are used as the explanatory variables and the overall survival month (OS mon) as the response variable. Patients with OS mon ≥ 20.3 are categorized into ‘longer overall survival’ group and patients with OS mon < 20.3 are categorized into ‘shorter overall survival’ group. The objective here is to investigate the TCR repertoire network properties and develop novel statistical method to prioritize the important network properties that are associated with the clinical outcome of increased overall survival.

### Challenges
The response variable, OS mon, has a defnite value for each of the 65 patients. The TCR network data (the explanatory variables) consists of a mix of global and local variables. The global variables are described by a single set of values, while the local variables are vectors of varying lengths. Since the TCR repertoire is constantly adapting to the health and the environmental factors of the patient, the network properties are continually shaping. Given any two patients the TCR repertoire is never the same. Less than 20% overlap is observed in the TCR repertoires for the same subject. This heterogeneous nature of the TCR repertoire and network properties makes it difcult to perform statistical inference or machine learning directly between subjects. The heterogeneity issue also complicates the data simulation process required to perform the simulation study. Therefore, we require to develop ingenious ways to handle the TCR network data throughout this work and derive meaningful inferences.

### 1.3 Contribution
- Strategy to extract features from heterogeneous TCR repertoire network data.
- A novel statistical method to prioritize the network properties - Group Lasso with Permutation Assisted Tuning.
- Mimicked real property distributions and correlation structure to simulate network properties.
- Demonstrated proposed methods and schemes using simulation study.

# Dataset Description
- 

# Data Visualization
<div align="center">
  <img src="" width=80% height=50%>
</div>

# Methods & Implementation


# Code (in R)
R code: 


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
