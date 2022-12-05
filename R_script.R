library(ggplot2)
library(ggtext)
library(plyr)
library(dplyr)
library(glmnet) 
library(MASS)
library(gglasso)
library(aod)
library(MLmetrics)
library(caret)
library(Rlab)

library(devtools)
#install_github("DataSlingers/ExclusiveLasso")
library(ExclusiveLasso)


main_function = function()
{
  
  ds_os = read.csv("C:\\Users\\sengu\\Dropbox\\PC\\Desktop\\SFSU_MasterThesis\\Load_dt_analysis.csv", 
                   header=TRUE, sep=",")
  
  #---------   
  # Dataframe with all CLUSTER INFO csv files and add patient_id column 
  #----------
  path = "C:\\Users\\sengu\\Dropbox\\PC\\Desktop\\SFSU_MasterThesis\\ClusterFiles\\"
  files = list.files(path, pattern = '\\.csv$', full.names = TRUE)
  
  all_data = do.call(rbind, lapply(files, function(x) transform(read.csv(x), File = basename(x))))
  
  all_data$File = gsub("([0-9]+)_.*", "\\1", all_data$File)
  names(all_data)[names(all_data) == "File"] = "Patient_id"
  
  
  #----------
  # Add column OS_MON to the cluster level dataframe    
  #----------
  
  lookup = ds_os %>%
    dplyr::select("Patient_id", "OS_mon")
  
  lookup = dplyr::select(ds_os, "Patient_id", "OS_mon")
  ds_data = join(all_data, lookup, by = "Patient_id")
  
  #----------
  # Add column RESPONSE #######
  # if OS_mon >= 20.3 Response = 1; else Response = 0    
  #----------
  ds_data$Response = -1
  ds_data$Response[ds_data$OS_mon >= 20.3] = 1
  

  
  #----------
  # Create response and feature matrices    
  #---------- 
  
  ##########  Response Variable Y  ###########
  mat_response = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(Response = mean(Response))
  #head(aggr_response)
  
    
  #######  Membership matrix ########
  aggr_cluster = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(count_cluster = n())
  #head(aggr_cluster)
  mat_cluster = as.matrix(aggr_cluster)
  
  #######  Node Count matrix ########
  
  summary(ds_data$node_count)
  
  aggr_node_count = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_node_count = min(node_count),
                     q1_node_count = quantile(node_count,0.25),
                     med_node_count = median(node_count),
                     mean_node_count = mean(node_count),
                     q3_node_count = quantile(node_count,0.75),
                     max_node_count = max(node_count))
  
  #head(aggr_node_count)
  mat_node_count = as.matrix(aggr_node_count)
  
  
  #######  Degree matrix ########
  aggr_deg = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_deg = min(deg),
                     q1_deg = quantile(deg,0.25),
                     med_deg = median(deg),
                     mean_deg = mean(deg),
                     q3_deg = quantile(deg,0.75),
                     max_deg = max(deg))
  #head(aggr_deg)
  mat_deg = as.matrix(aggr_deg)
  
  
  #######  AA_length matrix ########
  aggr_AA_len = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_AA_len = min(AA_length),
                     q1_AA_len = quantile(AA_length, 0.25),
                     med_AA_len = median(AA_length),
                     mean_AA_len = mean(AA_length),
                     q3_AA_len = quantile(AA_length, 0.75),
                     max_AA_len = max(AA_length))
  #head(aggr_AA_len)
  mat_AA_len = as.matrix(aggr_AA_len)
  
  
  #######  count_pre_infusion matrix ########
  aggr_cnt_pre_infu = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_cnt_pre_infu = min(Count_PRE_INFUSION),
                     q1_cnt_pre_infu = quantile(Count_PRE_INFUSION, 0.25),
                     med_cnt_pre_infu = median(Count_PRE_INFUSION),
                     mean_cnt_pre_infu = mean(Count_PRE_INFUSION),
                     q3_cnt_pre_infu = quantile(Count_PRE_INFUSION, 0.75),
                     max_cnt_pre_infu = max(Count_PRE_INFUSION))
  #head(aggr_cnt_pre_infu)
  mat_cnt_pre_infu = as.matrix(aggr_cnt_pre_infu)
  
  #######  count_dose_2 matrix ########
  aggr_cnt_dose2 = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_cnt_dose2 = min(Count_DOSE_2),
                     q1_cnt_dose2 = quantile(Count_DOSE_2, 0.25),
                     med_cnt_dose2 = median(Count_DOSE_2),
                     mean_cnt_dose2 = mean(Count_DOSE_2),
                     q3_cnt_dose2 = quantile(Count_DOSE_2, 0.75),
                     max_cnt_dose2 = max(Count_DOSE_2))
  #head(aggr_cnt_dose2)
  mat_cnt_dose2 = as.matrix(aggr_cnt_dose2)
  
  #######  deg_avg matrix ########
  aggr_deg_avg = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_deg_avg = min(deg_avg),
                     q1_deg_avg = quantile(deg_avg, 0.25),
                     med_deg_avg = median(deg_avg),
                     mean_deg_avg = mean(deg_avg),
                     q3_deg_avg = quantile(deg_avg, 0.75),
                     max_deg_avg = max(deg_avg))
  #head(aggr_deg_avg)
  mat_deg_avg = as.matrix(aggr_deg_avg)
  
  #######  Diameter Length matrix ########
  aggr_dia_len = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_dia_len = min(diam_length, na.rm=TRUE),
                     q1_dia_len = quantile(diam_length, na.rm=TRUE, 0.25),
                     med_dia_len = median(diam_length, na.rm=TRUE),
                     mean_dia_len = mean(diam_length, na.rm=TRUE),
                     q3_dia_len = quantile(diam_length, na.rm=TRUE, 0.75),
                     max_dia_len = max(diam_length, na.rm=TRUE))
  #head(aggr_dia_len)
  mat_dia_len = as.matrix(aggr_dia_len)
  
  #######  Assortativity matrix ########
  summary(all_data$assortativity)
  par(mfrow=c(1,1))
  hist(all_data$assortativity[all_data$Patient_id=="1093501642"], freq=FALSE,
       col = "#4371A6",
       xlab = "Assortativity", 
       ylab = "Density", main = "Assortativity of a Subject",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.8, cex.sub=1.5)
  hist(log(all_data$assortativity[all_data$Patient_id=="1093501642"]),col = "#4371A6",
       xlab = "Assortativity", 
       ylab = "Density", main = "Log-normal form",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.8, cex.sub=1.5)
  
  aggr_assort = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(prob_NA_assort = sum(is.na(assortativity))/n(),
                     min_assort = min(assortativity, na.rm=TRUE),
                     q1_assort = quantile(assortativity, na.rm=TRUE, 0.25),
                     med_assort = median(assortativity, na.rm=TRUE),
                     mean_assort = mean(assortativity, na.rm=TRUE),
                     q3_assort = quantile(assortativity, na.rm=TRUE, 0.75),
                     max_assort = max(assortativity, na.rm=TRUE))
  #head(aggr_assort)
  mat_assort = as.matrix(aggr_assort)
  
  
  #######  Transitivity matrix ########
  summary(all_data$transitivity)
  aggr_trans = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(prob_NA_trans = sum(is.na(transitivity))/n(),
                     min_trans = min(transitivity, na.rm=TRUE),
                     q1_trans = quantile(transitivity, na.rm=TRUE, 0.25),
                     med_trans = median(transitivity, na.rm=TRUE),
                     mean_trans = mean(transitivity, na.rm=TRUE),
                     q3_trans = quantile(transitivity, na.rm=TRUE, 0.75),
                     max_trans = max(transitivity, na.rm=TRUE))
  #head(aggr_trans)
  
  colnames(aggr_trans)
  aggr_trans%>%filter(Patient_id=="1093501642")
  
  par(mfrow=c(1,2))
  hist(all_data$transitivity[all_data$Patient_id=="1093501642"], freq=FALSE,
       col = "#4371A6",
       xlab = "Transitivity", 
       ylab = "Density", main = "Transitivity Distribution for a Subject",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.8, cex.sub=1.5)
  hist(log(all_data$transitivity[all_data$Patient_id=="1093501642"]),col = "#4371A6",
       xlab = "Log(Transitivity)", 
       ylab = "Density", main = "Log-normal form",
       cex.lab=1.5, cex.axis=1.5, cex.main=1.8, cex.sub=1.5)
  
  mat_trans = as.matrix(aggr_trans)
  
  #######  Edge_density matrix ########
  aggr_edg_den = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_edg_den = min(edge_density),
                     q1_edg_den = quantile(edge_density, 0.25),
                     med_edg_den = median(edge_density),
                     mean_edg_den= mean(edge_density),
                     q3_edg_den = quantile(edge_density, 0.75),
                     max_edg_den = max(edge_density))
  #head(aggr_edg_den)
  mat_edg_den = as.matrix(aggr_edg_den)
  
  ####### Centr_degreematrix ########
  aggr_cen_deg = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_cen_deg = min(centr_degree),
                     q1_cen_deg = quantile(centr_degree, 0.25),
                     med_cen_deg = median(centr_degree),
                     mean_cen_deg = mean(centr_degree),
                     q3_cen_deg = quantile(centr_degree, 0.75),
                     max_cen_deg = max(centr_degree))
  #head(aggr_cen_deg)
  mat_cen_deg = as.matrix(aggr_cen_deg)
  
  #######  Centr_clo matrix ########
  aggr_cen_clo = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(prob_NA_cen_clo = sum(is.na(centr_clo))/n(),
                     min_cen_clo = min(centr_clo, na.rm=TRUE),
                     q1_cen_clo = quantile(centr_clo, na.rm=TRUE, 0.25),
                     med_cen_clo = median(centr_clo, na.rm=TRUE),
                     mean_cen_clo = mean(centr_clo, na.rm=TRUE),
                     q3_cen_clo = quantile(centr_clo, na.rm=TRUE, 0.75),
                     max_cen_clo = max(centr_clo, na.rm=TRUE))
  #head(aggr_cen_clo)
  mat_cen_clo = as.matrix(aggr_cen_clo)
  
  ####### Eigen_centrality  matrix ########
  aggr_egn_cen = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(min_egn_cen = min(eigen_centrality),
                     q1_egn_cen = quantile(eigen_centrality, 0.25),
                     med_egn_cen = median(eigen_centrality),
                     mean_egn_cen = mean(eigen_centrality),
                     q3_egn_cen = quantile(eigen_centrality, 0.75),
                     max_egn_cen = max(eigen_centrality))
  #head(aggr_egn_cen)
  mat_egn_cen = as.matrix(aggr_egn_cen)
  
  #######  Centr_eigen matrix ########
  aggr_cen_egn = ds_data %>%
    group_by(Patient_id) %>%
    dplyr::summarise(prob_NA_cen_egn = sum(is.na(centr_eigen))/n(),
                     min_cen_egn = min(centr_eigen, na.rm=TRUE),
                     q1_cen_egn = quantile(centr_eigen, na.rm=TRUE, 0.25),
                     med_cen_egn = median(centr_eigen, na.rm=TRUE),
                     mean_cen_egn = mean(centr_eigen, na.rm=TRUE),
                     q3_cen_egn = quantile(centr_eigen, na.rm=TRUE, 0.75),
                     max_cen_egn = max(centr_eigen, na.rm=TRUE))
  #head(aggr_cen_egn)
  mat_cen_egn = as.matrix(aggr_cen_egn)
  
  
################################################################################  
  
  y_obsvd = mat_response$Response
   
  X = cbind(aggr_cluster[,-1], aggr_node_count[,-1], aggr_deg[,-1], aggr_AA_len[,-1], 
            aggr_cnt_pre_infu[,-1], aggr_cnt_dose2[,-1], aggr_deg_avg[,-1], 
            aggr_dia_len[,-1], aggr_assort[,-1], aggr_trans[,-1], aggr_edg_den[,-1], 
            aggr_cen_deg[,-1], aggr_cen_clo[,-1], aggr_egn_cen[,-1], aggr_cen_egn[,-1])
  
  #-----------
  # Standardized features 
  #-----------  
  X_scaled = apply(X, 2, function(y_obsvd) (y_obsvd - mean(y_obsvd)) / sd(y_obsvd) ^ as.logical(sd(y_obsvd)))
  
  #-----------
  # GROUP LASSO 
  #-----------
  set.seed(25)
  X = X_scaled
  
  # Cross Validation
  set.seed(25)
  grp_lasso_features = grp_lasso_cv(X, y_obsvd, 0.065) 
  grp_lasso_features
  grp_lasso_grp_idx = select_grps(grp_lasso_features)
  grp_lasso_grp_idx
  
  
  # Permutation Tuning
  set.seed(25)
  grp_plasso_features = grp_plassob(X, y=y_obsvd, pB=10, SS=0.25)
  grp_plasso_grp_idx = select_grps(grp_plasso_features)
  grp_plasso_grp_idx
  
  
  
  #-----------
  # LASSO 
  #-----------
  # Cross Validation
  set.seed(25)
  la_cv = cv.glmnet(X, y=y_obsvd, family="binomial", nfolds=5, 
                    alpha = 1)
  paste(round(la_cv$lambda.min,6), round(la_cv$lambda.1se,6))
  plot(la_cv)
  
  set.seed(25)
  la = glmnet(X, y=y_obsvd, family="binomial", lambda = la_cv$lambda.min, 
              alpha = 1)
  
  round(coef(la)[abs(coef(la)[,1]) > 0,],6)
  
  lasso_cv_feature_idx = NULL
  # extract non-zero coefficients 
  for(i in 2:length(coef(la))){
    if(abs(coef(la)[i,1])>0){
      lasso_cv_feature_idx = c(lasso_cv_feature_idx,i-1)
    }
  }
  lasso_cv_feature_idx
  
  
  # LASSO with permutation tuning
  set.seed(25)
  plasso_feature_idx = plassob(X_scaled, y_obsvd, pB=10, SS=0.1)
  plasso_feature_idx
  
  
  
  #------------
  # Exclusive Lasso with Cross Validation
  #------------
  exclsv_lasso_features_idx = exclsv_lasso(X, y_obsvd, nlambda = 50)
  exclsv_lasso_features_idx
  
#-------------------------  
# Simulate # of clusters in each sample
#-------------------------     
  set.seed(81)
  n = 1000
  
  membership_sample = ds_data%>%
    dplyr::group_by(Patient_id)%>%
    dplyr::summarise(max_mem = max(membership))
  
  summary(membership_sample$max_mem)
  
  
  par(mfrow=c(1,2))
  hist(membership_sample$max_mem, col = "#90C0AF", freq=FALSE,
       xlab = "# of Clusters", 
       ylab = "Density", 
       main = "Distribution of # of Clusters",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  
  hist(log(membership_sample$max_mem), col = "#D1BED1", freq=FALSE,
       xlab = "# (Log count) Clusters", 
       ylab = "Density", 
       main = "Distribution of (Log count) of Clusters",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)# approx normal
  
  
  mean_cluster_count = mean(log(membership_sample$max_mem))
  sd_cluster_count = sd(log(membership_sample$max_mem))
  
  cluster_sim = vector(length=n)
  i = 1
  while(i<=n){
    val = round(rlnorm(1, meanlog=mean_cluster_count, sdlog=sd_cluster_count),0)
    if(val>=15 && val<=900){
      cluster_sim[i]=val
      i = i+1
    }
  }
  
  summary(cluster_sim)
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_cluster[,2]), freq=FALSE, col = "#FFFAC2",
       xlab = "Observed Cluster size per patient/sample from Real Data", 
       ylab = "Density", 
       main = "Density of Observed Cluster size",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  
  hist(as.numeric(cluster_sim), freq=FALSE, col = "#BDE496",
       xlab = "Simulated Cluster size per patient/sample", 
       ylab = "Density", 
       main = "Density of Simulated Cluster size",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  
#-------------------------  
# Simulate Node count data
#-------------------------  
  summary(ds_data$node_count)
  
  # Overall Node_count distribution 
  aggr_nodes = ds_data%>%
    dplyr::summarize(prob_2 = sum(node_count == 2)/n(), 
                     prob_num = sum(node_count > 2)/n(),
                     prob_3 = sum(node_count == 3)/n())
  sum(aggr_nodes)
  aggr_nodes
  
  hist(ds_data$node_count[ds_data$node_count>3], freq=FALSE)
  hist(log(ds_data$node_count[ds_data$node_count>3]), freq=FALSE)
  
  mean_node_grt_2 = mean(ds_data$node_count[ds_data$node_count>2])
  sd_node_grt_2 = sd((ds_data$node_count[ds_data$node_count>2]))

  nodes_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    nodes_per_cluster = vector(length=node_len)
    num_2 = round(cluster_sim[i]*aggr_nodes[1],0)
    num_3 = round(cluster_sim[1]*aggr_nodes[3],0)
    
    j = 1
    while(j <= num_2){
      nodes_per_cluster[j] = 2
      j=j+1
    }
    
    while(j <= sum(num_2,num_3)){
      nodes_per_cluster[j] = 3
      j=j+1
    }
    
    while(j <= node_len){
      val = round(rlnorm(1, meanlog=log(mean_node_grt_2), sdlog=log(sd_node_grt_2)/3),0)
      if(val<= 404 && val>=3){
        nodes_per_cluster[j] = val
        j=j+1  
      }
      
    }
    
    for(k in 1:6){
      nodes_sim[i,k] = summary(nodes_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_node_count[, 2:7]), freq=FALSE, col = "#FFFAC2", 
       main="Density of Observed Node Count data", xlab="Observed ")
  hist(as.numeric(nodes_sim), freq=FALSE, 
       main="Density of Simulated Node Count data", 
       col = "#BDE496", xlab="SimulatedNode Count")
  

#-------------------------  
# Simulate Degree data
#------------------------- 
  set.seed(81)
  summary(ds_data$deg)
  cor(ds_data$node_count, ds_data$deg)
  unique(ds_data$node_count[ds_data$deg == 1])
  
  deg_prob = ds_data%>%
    dplyr::summarize(prob_1 = sum(deg == 1)/n(), 
                     prob_num = sum(deg > 1)/n())
  sum(deg_prob)
  
  hist(ds_data$deg[ds_data$deg > 1])
  hist(log(ds_data$deg[ds_data$deg > 1]))
  
  deg_mean = mean(log(ds_data$deg[ds_data$deg > 1]))
  deg_sd = sd(log(ds_data$deg[ds_data$deg > 1]))
  
  par(mfrow=c(1,1))
  hist(ds_data$deg[ds_data$deg>1], freq = FALSE,
       xlim=c(1,5), ylim=c(0,2.5), main="# of records with Degree > 1",
       xlab="Degree value", col = "lightblue")
  
  curve(dlnorm(x, meanlog=deg_mean, sdlog=deg_sd), from=1, to=5, add=TRUE, col="red", lwd=2)
  
  
  # Data simulation for Degree
  deg_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    deg_per_cluster = vector(length=node_len)
    node_2 = round(cluster_sim[i]*aggr_nodes[1],0)
    
    j = 1
    while(j <= node_2){
      deg_per_cluster[j] = 1
      j=j+1
    }
    
    while(j <= cluster_sim[i]){
      #val = rlnorm(1, meanlog=deg_mean-1.7, sdlog=deg_sd)
      val = rlnorm(1, meanlog=deg_mean, sdlog=deg_sd)
      if(val>1 && val<5){
        deg_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    for(k in 1:6){
      deg_sim[i,k] = summary(deg_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_deg[,2:7]), freq=FALSE, col="lightyellow", 
       main="Density of Observed Degree data", xlab="Observed Degree")
  hist(as.numeric(deg_sim), freq=FALSE, main="Density of Simulated Degree data", 
       col="lightgreen", xlab="Simulated Degree")
  
 
#-------------------------  
# Simulate AA Length data
#------------------------- 
  set.seed(81)
  summary(ds_data$AA_length)
  cor(ds_data$node_count, ds_data$AA_length)
  #unique(ds_data$node_count[ds_data$AA_length == 10])
  #unique(ds_data$AA_length)
  
  hist(ds_data$AA_length)
  
  aa_prob = ds_data%>%
    dplyr::summarize(prob_lessthan_8 = sum(AA_length < 8)/n(), 
                     prob_norm = sum(AA_length >= 8 & AA_length <= 15)/n() ,
                     prob_morethan_15 = sum(AA_length > 15)/n())
  sum(aa_prob)
  
  
  hist(ds_data$AA_length[ds_data$AA_length >= 8 && ds_data$AA_length <= 15])
  
  mean_aa_norm = mean(ds_data$AA_length[ds_data$AA_length >= 8 && ds_data$AA_length <= 15])
  sd_aa_norm = sd(ds_data$AA_length[ds_data$AA_length >= 8 && ds_data$AA_length <= 15])
  
  
  # Data simulation for AA Length
  aa_len_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    aa_len_per_cluster = vector(length=node_len)
    aa_norm = round(aa_prob[2]*cluster_sim[i],0)
    
    j = 1
    while(j <= node_len){
      val = rnorm(1, mean=mean_aa_norm, sd=sd_aa_norm+1)
      if(val>=4 && val<=17){
        aa_len_per_cluster[j] = val
        j=j+1
      }
      
    }
    
    for(k in 1:6){
      aa_len_sim[i,k] = summary(aa_len_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_AA_len[,2:7]), freq=FALSE, col="lightyellow", 
       main="Density of Observed AA Length data", xlab="Observed AA Length")
  hist(as.numeric(aa_len_sim), freq=FALSE, main="Density of Simulated AA Length data", 
       col="lightgreen", xlab="Simulated AA Length")
  
  
#-------------------------  
# Simulate Count Pre Infusion data
#------------------------- 
  set.seed(81)
  summary(ds_data$Count_PRE_INFUSION)
  cor(ds_data$node_count, ds_data$Count_PRE_INFUSION)
  #unique(ds_data$Count_PRE_INFUSION)
  
  hist(ds_data$Count_PRE_INFUSION)
  
  pre_infu_prob = ds_data%>%
    dplyr::summarize(prob_0 = sum(Count_PRE_INFUSION == 0)/n(), 
                     prob_num = sum(Count_PRE_INFUSION > 100)/n())
  sum(pre_infu_prob)
  pre_infu_prob
  
  hist(log(ds_data$Count_PRE_INFUSION), freq = FALSE)
  
  mean_pre_inf = mean(log(ds_data$Count_PRE_INFUSION[ds_data$Count_PRE_INFUSION>0]))
  sd_pre_inf = sd(log(ds_data$Count_PRE_INFUSION[ds_data$Count_PRE_INFUSION>0]))
  
  
  # Data simulation for Count Pre Infusion
  pre_infu_sim = matrix(, nrow = n, ncol = 6)
  
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    pre_infu_per_cluster = vector(length=node_len)
    num_0 = round(cluster_sim[i]*pre_infu_prob[1],0)
    
    j = 1
    while(j <= num_0){
      pre_infu_per_cluster[j] = 0
      j=j+1
    }
    
    while(j <= node_len){
      val = rlnorm(1, meanlog = mean_pre_inf+2.5, sdlog = sd_pre_inf-0.09)
      if(val>0 && val<=410000){
        pre_infu_per_cluster[j] = val
        j=j+1
      }
    }
    
    for(k in 1:6){
      pre_infu_sim[i,k] = summary(pre_infu_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_cnt_pre_infu[,2:7]), freq=FALSE, col="lightyellow", 
       main="Density of Observed Count Pre Infusion data", xlab="Observed Count Pre Infusion")
  hist(as.numeric(pre_infu_sim), freq=FALSE, main="Density of Simulated Count Pre Infusion data", 
       col="lightgreen", xlab="Simulated Count Pre Infusion")
  
 
#-------------------------  
# Simulate Count dose 2 data
#------------------------- 
  set.seed(81)
  summary(ds_data$Count_DOSE_2)
  cor(ds_data$node_count, ds_data$Count_DOSE_2)
  #unique(ds_data$Count_DOSE_2)
  
  hist(log(ds_data$Count_DOSE_2[ds_data$Count_DOSE_2 > 80000]))
  boxplot(log(ds_data$Count_DOSE_2))
  
  dose2_prob = ds_data%>%
    dplyr::summarize(prob_0 = sum(Count_DOSE_2 == 0)/n(), 
                     prob_num = sum(Count_DOSE_2 > 0)/n())
  sum(dose2_prob)
  dose2_prob
  
  mean_dose = mean(log(ds_data$Count_DOSE_2[ds_data$Count_DOSE_2 > 0]))
  sd_dose = sd(log(ds_data$Count_DOSE_2[ds_data$Count_DOSE_2 > 0]))
  
  
  # Data simulation for Count Dose 2
  dose2_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    dose2_per_cluster = vector(length=node_len)
    num_0 = round(cluster_sim[i]*dose2_prob[1],0)
    
    j = 1
    while(j <= num_0){
      dose2_per_cluster[j] = 0
      j=j+1
    }
    
    while(j <= node_len){
      val = rlnorm(1, meanlog = mean_dose+4, sdlog = sd_dose-0.7)
        if(val>=0 && val<=300000){
          dose2_per_cluster[j] = val
          j=j+1
        }
     }
    
    for(k in 1:6){
      dose2_sim[i,k] = summary(dose2_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_cnt_dose2[,2:7]), freq=FALSE, col="lightyellow", 
       main="Density of Observed Count Dose 2 data", xlab="Observed Count Dose 2")
  hist(as.numeric(dose2_sim), freq=FALSE, main="Density of Simulated Count Dose 2 data", 
       col="lightgreen", xlab="Simulated Count Dose 2")
  
  
#-------------------------  
# Simulate Degree Average data
#------------------------- 
  set.seed(81)
  summary(ds_data$deg_avg)
  cor(ds_data$node_count, ds_data$deg_avg)
  #unique(ds_data$node_count[ds_data$deg_avg == 1])
  
  davg_prob = ds_data%>%
    dplyr::summarize(prob_1 = sum(deg_avg == 1)/n(), 
                     prob_num = sum(deg_avg > 1)/n())
  sum(davg_prob)
  davg_prob
  
  hist(ds_data$deg_avg[ds_data$deg_avg > 1])
  hist(log(ds_data$deg_avg[ds_data$deg_avg > 1]))
  
  mean_davg_grt1 = mean(log(ds_data$deg_avg[ds_data$deg_avg > 1]))
  sd_davg_grt1 = sd(log(ds_data$deg_avg[ds_data$deg_avg > 1]))
  
  par(mfrow=c(1,1))
  hist(ds_data$deg_avg[ds_data$deg_avg>1.5], freq = FALSE)
  
  
  # Data simulation for Degree Average
  davg_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    davg_per_cluster = vector(length=node_len)
    num_2 = round(cluster_sim[i]*aggr_nodes[1],0)
    
    j = 1
    while(j <= num_2){
      davg_per_cluster[j] = 1
      j=j+1
    }
    
    while(j <= cluster_sim[i]){
      val = rlnorm(1, meanlog=mean_davg_grt1, sdlog=sd_davg_grt1)
      if(val>1 && val<5){
        davg_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    for(k in 1:6){
      davg_sim[i,k] = summary(davg_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_deg_avg[,2:7]), freq=FALSE, col="lightyellow", 
       main="Density of Observed Degree Avg. data", xlab="Observed Degree Avg.")
  hist(as.numeric(davg_sim), freq=FALSE, main="Density of Simulated Degree Avg. data", 
       col="lightgreen", xlab="Simulated Degree Avg.")
  
 
  
  
#-------------------------  
# Simulate Diameter Length data
#------------------------- 
  set.seed(81)
  summary(ds_data$diam_length)
  cor(ds_data$node_count, ds_data$diam_length)
  #unique(ds_data$node_count[ds_data$deg_avg == 1])
  
  dialen_prob = ds_data%>%
    dplyr::summarize(prob_1 = sum(diam_length == 2)/n(), 
                     prob_num = sum(diam_length > 2)/n())
  sum(dialen_prob)
  dialen_prob
  
  hist(ds_data$diam_length[ds_data$diam_length > 2])
  hist(log(ds_data$diam_length[ds_data$diam_length > 2]))
  
  mean_dialen = mean(log(ds_data$diam_length[ds_data$diam_length > 2]))
  sd_dialen = sd(log(ds_data$diam_length[ds_data$diam_length > 2]))
  
  par(mfrow=c(1,1))
  
  
# Data simulation for Diameter Length
  dialen_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    #i = 1
    node_len = cluster_sim[i]
    dialen_per_cluster = vector(length=node_len)
    num_2 = round(cluster_sim[i]*aggr_nodes[1],0)
    
    j = 1
    while(j <= num_2){
      dialen_per_cluster[j] = 2
      j=j+1
    }
    
    while(j <= cluster_sim[i]){
      val = rlnorm(1, meanlog=mean_dialen+0.25, sdlog=sd_dialen+.01)
      if(val>2 && val<22){
        dialen_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    for(k in 1:6){
      dialen_sim[i,k] = summary(dialen_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_dia_len[,2:7]), freq=FALSE, col="lightyellow", 
       main="Density of Observed Diameter Length data", xlab="Observed Diameter Length")
  hist(as.numeric(dialen_sim), freq=FALSE, main="Density of Simulated Diameter Length data", 
       col="lightgreen", xlab="Simulated Diameter Length")
  
  
  
  
#--------------
# Diameter Length Data Simulation
#--------------  
  summary(ds_data$diam_length)
  par(mfrow=c(1,1))
  hist(ds_data$diam_length, freq=FALSE, col="lightyellow")
  unique(ds_data$diam_length)
  
  dia_prob = ds_data%>%
    dplyr::summarize(prob_2 = sum(diam_length == 2)/n(), 
                     prob_3 = sum(diam_length == 3)/n(),
                     prob_4 = sum(diam_length == 4)/n(),
                     prob_num = sum(diam_length > 4)/n())
  sum(dia_prob)
  dia_prob
  
  hist(log(ds_data$diam_length[ds_data$diam_length > 4]))
  
  mean_dia = mean(log(ds_data$diam_length[ds_data$diam_length > 4]))
  sd_dia = sd(log(ds_data$diam_length[ds_data$diam_length > 4]))
  
  
  # Data simulation for Diameter Length
  dia_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    dia_per_cluster = vector(length=node_len)
    num_2 = round(cluster_sim[i]*dia_prob[1],0)
    num_3 = round(cluster_sim[i]*dia_prob[2],0)
    num_4 = round(cluster_sim[i]*dia_prob[3],0)
    
    j = 1
    while(j <= num_2){
      dia_per_cluster[j] = 2
      j=j+1
    }
    while(j <= sum(num_2,num_3)){
      dia_per_cluster[j] = 3
      j=j+1
    }
    while(j <= sum(num_2,num_3,num_4)){
      dia_per_cluster[j] = 4
      j=j+1
    }
    
    while(j <= cluster_sim[i]){
      val = round(rlnorm(1, meanlog=mean_dia-0.25, sdlog=sd_dia+0.01),0)
      if(val>4 && val<=21){
        dia_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    for(k in 1:6){
      dia_sim[i,k] = summary(dia_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_dia_len[,2:7]), freq=FALSE, col="lightyellow", 
       main="Density of Observed Diameter Length data", xlab="Observed Diameter Length")
  hist(as.numeric(dia_sim), freq=FALSE, main="Density of Simulated Diameter Length data", 
       col="lightgreen", xlab="Simulated Diameter Length")
  
  
#-------------------------  
# Simulate Assortativity data
#------------------------- 
  summary(ds_data$assortativity)
  
  assort_prob = ds_data%>%
    dplyr::summarize(prob_NA = sum(is.na(assortativity))/n(), 
                     prob_neg1 = sum(assortativity == -1, na.rm = T)/n(),
                     prob_num = sum(assortativity > -1, na.rm = T)/n())
  sum(assort_prob)
  assort_prob
  
  hist(ds_data$assortativity[ds_data$assortativity > -1])
  hist(log(ds_data$assortativity[ds_data$assortativity > -1]))
  
  par(mfrow=c(1,2))
  hist(ds_data$assortativity[ds_data$assortativity > -1], col = "#90C0AF",
       xlab = "Assortativity", 
       ylab = "Frequency", 
       main = "Assortativity > -1",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  
  hist(log(ds_data$assortativity[ds_data$assortativity > -1]), col = "#90C0AF",
       xlab = "# (Log) Assortativity", 
       ylab = "Frequency", 
       main = "Frequency of clusters per subject",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)# approx normal
  
  assort_mean = mean(ds_data$assortativity[ds_data$assortativity > -1], na.rm=T)
  assort_sd = sd(ds_data$assortativity[ds_data$assortativity > -1], na.rm=T)
  
  
  # Data simulation for Assortativity
  set.seed(81)
  assort_sim = matrix(, nrow = n, ncol = 7)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    assort_per_cluster = vector(length=node_len)
    num_2 = round(cluster_sim[i]*aggr_nodes[1],0)
    num_3 = round(cluster_sim[i]*aggr_nodes[3],0)
    
    j = 1
    while(j <= num_2){
      assort_per_cluster[j] = NA
      j=j+1
    }
    
    while(j <= sum(num_2,num_3)){
      assort_per_cluster[j] = -1
      j=j+1
    }
    
    while(j <= node_len){
      val = rlnorm(1, meanlog=assort_mean, sdlog=assort_sd)-1
      if(val> -1 && val<= 1){
        assort_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    assort_sim[i,1] = summary(assort_per_cluster)[7]/node_len
    for(k in 2:7){
      assort_sim[i,k] = summary(assort_per_cluster)[k-1]
    }
  }
  
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_assort[,2:8]), col = "#FFFAC2",freq=FALSE,
       xlab = "Observed Assortativity from Real Data", 
       ylab = "Density", 
       main = "Density of Observed Assortativity data",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  
  hist(as.numeric(assort_sim), freq=FALSE, col = "#BDE496",
       xlab = "Simulated Assortativity", 
       ylab = "Density", 
       main = "Density of Simulated Assortativity data",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  
  
  
  
#-------------------------  
# Simulate Transitivity data
#------------------------- 
  set.seed(81)
  summary(ds_data$transitivity)
  unique(ds_data$node_count[ds_data$transitivity == 0])
  unique(ds_data$transitivity[ds_data$node_count == 3])
  
  
  transtv_prob = ds_data%>%
    dplyr::summarize(prob_NA = sum(is.na(transitivity))/n(), 
                     prob_0 = sum(transitivity == 0, na.rm = T)/n(), 
                     prob_1 = sum(transitivity == 1, na.rm = T)/n(),
                     prob_num = sum(transitivity > 0, na.rm = T)/n() 
                     - sum(transitivity == 1, na.rm = T)/n())
  sum(transtv_prob)
  transtv_prob
  
  hist(ds_data$transitivity[ds_data$transitivity > 0 & ds_data$transitivity != 1])
  hist(log(ds_data$transitivity[ds_data$transitivity > 0 & ds_data$transitivity != 1]))
  
  trans_mean = mean(log(ds_data$transitivity[ds_data$transitivity > 0 
                                             & ds_data$transitivity != 1]), 
                    na.rm=T)
  trans_sd = sd(log(ds_data$transitivity[ds_data$transitivity > 0 
                                         & ds_data$transitivity != 1]), 
                na.rm=T)
  
  
  # Data simulation for Transitivity
  transtv_sim = matrix(, nrow = n, ncol = 7)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    transtv_per_cluster = vector(length=node_len)
    node_2 = round(cluster_sim[i]*aggr_nodes[1],0)
    
    j = 1
    while(j <= node_2){
      transtv_per_cluster[j] = NA
      j=j+1
    }
    
    num_0 = round(cluster_sim[i]*transtv_prob[2],0)
    num_1 = round(cluster_sim[i]*transtv_prob[3],0)
    
    
    while(j <= sum(node_2+num_0)){
      transtv_per_cluster[j] = 0
      j=j+1
    }
    
    while(j <= sum(node_2 + num_0 + num_1)){
      transtv_per_cluster[j] = 1
      j=j+1
    }
    
    while(j <= cluster_sim[i]){
      #val = rlnorm(1, meanlog=trans_mean+0.1, sdlog=trans_sd)-1
      val = rlnorm(1, meanlog=trans_mean+0.5, sdlog=trans_sd+0.3)
      if(val>0 && val<1){
        transtv_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    transtv_sim[i,1] = summary(transtv_per_cluster)[7]/node_len
    for(k in 2:7){
      transtv_sim[i,k] = summary(transtv_per_cluster)[k-1]
    }
  }
  
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_trans[,2:8]), col = "#FFFAC2",freq=FALSE,
       xlab = "Observed Transitivity from Real Data", 
       ylab = "Density", 
       main = "Density of Observed Transitivity data",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  
  hist(as.numeric(transtv_sim), freq=FALSE, col = "#BDE496",
       xlab = "Simulated Transitivity", 
       ylab = "Density", 
       main = "Density of Simulated Transitivity data",
       cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2)
  
  
  
#--------------
# Edge Density Data Simulation
#--------------  
  set.seed(81)
  summary(ds_data$edge_density)
  unique(ds_data$edge_density)
  par(mfrow=c(1,1))
  hist(ds_data$edge_density, freq=FALSE, col="lightyellow")
  
  distinct(ds_data, edge_density)
  
  edge_den_prob = ds_data%>%
    dplyr::summarize(prob_1 = sum(edge_density == 1)/n(), 
                     prob_num = sum(edge_density < 1)/n())
  sum(edge_den_prob)
  edge_den_prob
  
  hist(log(ds_data$edge_density[ds_data$edge_density < 1]))
  
  
  mean_edg_den = mean(log(ds_data$edge_density[ds_data$edge_density < 1]))
  sd_edg_den = sd(log(ds_data$edge_density[ds_data$edge_density < 1]))
  
  
  # Edge Density Data simulation for 
  edge_den_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    edg_den_per_cluster = vector(length=node_len)
    num_1 = round(cluster_sim[i]*edge_den_prob[1],0)
    
    j = 1
    while(j <= num_1){
      edg_den_per_cluster[j] = 1
      j=j+1
    }
   
    while(j <= cluster_sim[i]){
      val = rlnorm(1, meanlog=mean_edg_den, sdlog=sd_edg_den-0.1)
      if(val>0 && val<1){
        edg_den_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    for(k in 1:6){
      edge_den_sim[i,k] = summary(edg_den_per_cluster)[k]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_edg_den[,2:7]), freq=FALSE, col="lightyellow",
       main="Density of Observed Edge Density data", xlab="Observed Edge Density")
  hist(as.numeric(edge_den_sim), freq=FALSE, main="Density of Simulated Edge Density data", col="lightgreen"
       , xlab="Simulated Edge Density")
  
  
#--------------
# Center Degree Data Simulation
#--------------  
  set.seed(81)
  summary(ds_data$centr_degree)
  par(mfrow=c(1,1))
  hist(ds_data$centr_degree, freq=FALSE, col="lightyellow")
  
  cen_deg_prob = ds_data%>%
    dplyr::summarize(prob_0 = sum(centr_degree == 0, na.rm = T)/n(), 
                     prob_num = sum(centr_degree > 0, na.rm = T)/n())
  sum(cen_deg_prob)
  cen_deg_prob
  
  hist(log(ds_data$centr_degree[ds_data$centr_degree > 0],))
  
  mean_cen_deg = mean(log(ds_data$centr_degree[ds_data$centr_degree > 0]), na.rm=T)
  sd_cen_deg = sd(log(ds_data$centr_degree[ds_data$centr_degree > 0]), na.rm=T)
  
  
  #  Data simulation for Central Degree
  cen_deg_sim = matrix(, nrow = n, ncol = 6)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    cen_deg_per_cluster = vector(length=node_len)
    num_10 = round(cluster_sim[i]*cen_deg_prob[1],0)
    
    j = 1
    while(j <= num_0){
      cen_deg_per_cluster[j] = 0
      j=j+1
    }
    
    while(j <= cluster_sim[i]){
      val = rlnorm(1, meanlog=mean_cen_deg, sdlog=sd_cen_deg+0.7)
      if(val>0 && val<0.7){
        cen_deg_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    for(k in 1:6){
      cen_deg_sim[i,k] = summary(cen_deg_per_cluster)[k]
    }
  }
  
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_cen_deg[,2:7]), freq=FALSE, col="lightyellow",
       main="Density of Observed Central Degree data", xlab="Observed Central Degree")
  hist(as.numeric(cen_deg_sim), freq=FALSE, main="Density of Simulated Central Degree data",
       col="lightgreen", xlab="Simulated Central Degree")
  
  
  
#--------------
# Centr Closeness Data Simulation
#--------------  
  set.seed(81)
  summary(ds_data$centr_clo)
  par(mfrow=c(1,1))
  hist(ds_data$centr_clo, freq=FALSE, col="lightyellow")
  
  clo_prob = ds_data%>%
    dplyr::summarize(prob_NA = sum(is.na(centr_clo))/n(), 
                     prob_num = sum(centr_clo < 1, na.rm = T)/n(),
                     prob_1 = sum(centr_clo == 1, na.rm = T)/n())
  sum(clo_prob)
  clo_prob
  
  
  hist(ds_data$centr_clo[ds_data$centr_clo <1])
  
  clo_mean = mean(ds_data$centr_clo[ds_data$centr_clo <1], na.rm=TRUE)
  clo_sd = sd(ds_data$centr_clo[ds_data$centr_clo < 1], na.rm=TRUE)
  
  
  #Data simulation for Centr_clo
  clo_sim = matrix(, nrow = n, ncol = 7)
  
  for(i in 1:n){
    node_len = cluster_sim[i]
    clo_per_cluster = vector(length=node_len)
    num_na = round(cluster_sim[i]*clo_prob[1],0)
    num_1 = round(cluster_sim[i]*clo_prob[3],0)
    
    j = 1
    while(j <= num_na){
      clo_per_cluster[j] = NA
      j=j+1
    }
    
    while(j <= sum(num_na,num_1)){
      clo_per_cluster[j] = 1
      j=j+1
    }
    
    while(j <= cluster_sim[i]){
      val = rnorm(1, mean=clo_mean, sd=clo_sd+0.2)
      if(val>=0 && val<1){
        clo_per_cluster[j] = val
        j=j+1 
      } 
    }
    
    clo_sim[i,1] = summary(clo_per_cluster)[7]/node_len
    for(k in 2:7){
      clo_sim[i,k] = summary(clo_per_cluster)[k-1]
    }
  }
  
  par(mfrow=c(2,1))
  hist(as.numeric(mat_cen_clo[,2:8]), freq=FALSE, col="lightyellow", xlim=c(0,1),
       main="Density of Observed Centr_clo data", xlab="Observed Centr_clo")
  hist(as.numeric(clo_sim), freq=FALSE, main="Density of Simulated Centr_clo data", xlim=c(0,1), 
       col="lightgreen", xlab="Simulated Centr_clo")
  
  
  
#--------------
# Eigen Centrality Data Simulation
#--------------  
  set.seed(81)
  summary(ds_data$eigen_centrality)
  par(mfrow=c(1,1))
  hist(ds_data$eigen_centrality, freq=FALSE, col="lightyellow")
  
  egn_cen_prob = ds_data%>%
    dplyr::summarize(prob_1 = sum(eigen_centrality == 1)/n(), 
                     prob_num = sum(eigen_centrality > 1)/n())
  sum(egn_cen_prob)
  egn_cen_prob
  
  hist(log(ds_data$eigen_centrality[ds_data$eigen_centrality > 1]))
  
  mean_egn_cen = mean(log(ds_data$eigen_centrality[ds_data$eigen_centrality > 1]))
  sd_egn_cen = sd(log(ds_data$eigen_centrality[ds_data$eigen_centrality > 1]))
  
 
#  Data simulation for Eigen Centrality
egn_cen_sim = matrix(, nrow = n, ncol = 6)

for(i in 1:n){
  node_len = cluster_sim[i]
  egn_cen_per_cluster = vector(length=node_len)
  num_1 = round(cluster_sim[i]*egn_cen_prob[1],0)
  
  j = 1
  while(j <= num_1){
    egn_cen_per_cluster[j] = 1
    j=j+1
  }
  
  while(j <= cluster_sim[i]){
    val = rlnorm(1, meanlog=mean_egn_cen+.7, sdlog=sd_egn_cen-0.08)
    if(val>1 && val<11.5){
      egn_cen_per_cluster[j] = val
      j=j+1 
    } 
  }
  
  for(k in 1:6){
    egn_cen_sim[i,k] = summary(egn_cen_per_cluster)[k]
  }
}



par(mfrow=c(2,1))
hist(as.numeric(mat_egn_cen[,2:7]), freq=FALSE, col="lightyellow", 
     main="Density of Observed Eigen Centrality data", xlab="Observed Eigen Centrality")
hist(as.numeric(egn_cen_sim), freq=FALSE,  main="Density of Simulated Eigen Centrality data", 
     col="lightgreen", xlab="Simulated Eigen Centrality")


#--------------
# Central Eigen Data Simulation
#--------------  
set.seed(81)
summary(ds_data$centr_eigen)
par(mfrow=c(1,1))
hist(ds_data$centr_eigen, freq=FALSE, col="lightyellow")
hist(ds_data$centr_eigen[ds_data$centr_eigen>0], freq=FALSE, col="lightyellow")

cen_egn_prob = ds_data%>%
  dplyr::summarize(prob_NA = sum(is.na(centr_eigen))/n(), 
                   prob_num = sum(centr_eigen > 0, na.rm = T)/n(),
                   prob_0 = sum(centr_eigen == 0, na.rm = T)/n())
sum(cen_egn_prob)
cen_egn_prob

mean_cen_egn = mean(ds_data$centr_eigen, na.rm=TRUE)
sd_cen_egn = sd(ds_data$centr_eigen, na.rm=TRUE)


#Data simulation for Central Eigen
cen_egn_sim = matrix(, nrow = n, ncol = 7)

for(i in 1:n){
  node_len = cluster_sim[i]
  cen_egn_per_cluster = vector(length=node_len)
  num_na = round(cluster_sim[i]*cen_egn_prob[1],0)
  num_0 = round(cluster_sim[i]*cen_egn_prob[3],0)
  
  j = 1
  while(j <= num_na){
    cen_egn_per_cluster[j] = NA
    j=j+1
  }
  
  while(j <= sum(num_na,num_0)){
    cen_egn_per_cluster[j] = 0
    j=j+1
  }
  
  while(j <= cluster_sim[i]){
    val = rnorm(1, mean=mean_cen_egn-0.02, sd=sd_cen_egn-0.06)
    if(val>0 && val<1){
      cen_egn_per_cluster[j] = val
      j=j+1 
    } 
  }
  
  cen_egn_sim[i,1] = summary(cen_egn_per_cluster)[7]/node_len
  for(k in 2:7){
    cen_egn_sim[i,k] = summary(cen_egn_per_cluster)[k-1]
  }
}

par(mfrow=c(2,1))
hist(as.numeric(mat_cen_egn[,2:8]), freq=FALSE, col="lightyellow", xlim=c(0,1),
     main="Density of Observed Central Eigen data", xlab="Observed Central Eigen")
hist(as.numeric(cen_egn_sim), freq=FALSE, main="Density of Simulated Central Eigen data", xlim=c(0,1), 
     col="lightgreen", xlab="Simulated Central Eigen")

#########################################################

#-------------- 
# Bind simulated matrices
#-------------- 
X_sim = cbind(cluster_sim, nodes_sim, deg_sim, aa_len_sim, pre_infu_sim, 
              dose2_sim, davg_sim, dia_sim, assort_sim, transtv_sim, 
              edge_den_sim, cen_deg_sim, clo_sim, egn_cen_sim, cen_egn_sim)

colnames(X_sim) = c("count_cluster",     "min_node_count",    "q1_node_count" ,    "med_node_count",    "mean_node_count"  ,
                    "q3_node_count",     "max_node_count" ,   "min_deg"  ,         "q1_deg"    ,        "med_deg"   ,       
                    "mean_deg"  ,        "q3_deg"   ,         "max_deg"  ,         "min_AA_len"  ,      "q1_AA_len"  ,      
                    "med_AA_len" ,       "mean_AA_len" ,      "q3_AA_len"  ,       "max_AA_len"  ,      "min_cnt_pre_infu" ,
                    "q1_cnt_pre_infu",   "med_cnt_pre_infu" , "mean_cnt_pre_infu", "q3_cnt_pre_infu",   "max_cnt_pre_infu" ,
                    "min_cnt_dose2" ,    "q1_cnt_dose2" ,     "med_cnt_dose2" ,    "mean_cnt_dose2" ,   "q3_cnt_dose2"   ,  
                    "max_cnt_dose2",     "min_deg_avg" ,      "q1_deg_avg"  ,      "med_deg_avg" ,      "mean_deg_avg"  ,   
                    "q3_deg_avg" ,       "max_deg_avg",       "min_dia_len" ,      "q1_dia_len" ,       "med_dia_len"   ,   
                    "mean_dia_len" ,     "q3_dia_len" ,       "max_dia_len" ,      "prob_NA_assort" ,   "min_assort"   ,    
                    "q1_assort" ,        "med_assort",        "mean_assort" ,      "q3_assort" ,        "max_assort"  ,     
                    "prob_NA_trans",     "min_trans" ,        "q1_trans" ,         "med_trans"  ,       "mean_trans" ,      
                    "q3_trans" ,         "max_trans" ,        "min_edg_den" ,      "q1_edg_den"  ,      "med_edg_den"  ,    
                    "mean_edg_den" ,     "q3_edg_den"  ,      "max_edg_den",       "min_cen_deg"  ,     "q1_cen_deg"    ,   
                    "med_cen_deg" ,      "mean_cen_deg" ,     "q3_cen_deg" ,       "max_cen_deg"  ,     "prob_NA_cen_clo"  ,
                    "min_cen_clo",       "q1_cen_clo" ,       "med_cen_clo" ,      "mean_cen_clo"  ,    "q3_cen_clo"  ,     
                    "max_cen_clo",       "min_egn_cen" ,      "q1_egn_cen" ,       "med_egn_cen"  ,     "mean_egn_cen" ,    
                    "q3_egn_cen" ,       "max_egn_cen" ,      "prob_NA_cen_egn" ,  "min_cen_egn"  ,     "q1_cen_egn" ,      
                    "med_cen_egn" ,      "mean_cen_egn",      "q3_cen_egn" ,       "max_cen_egn")

#-----------
# Standardize the simulated data
#-----------  
X_sim_scaled = apply(X_sim, 2, function(y_lim) (y_lim - mean(y_lim)) / sd(y_lim) ^ as.logical(sd(y_lim)))

#View(X_sim_scaled)


#--------------
# Response Variable Data Simulation
#--------------  
# y is Bernoulli rvs
# y = 1 with probability pr(i)
# y = 0 with probability 1-pr(i)
# P(Y=1|X=i) = pr(i) = exp(z)/(1+exp(z))
# P(Y=0|X=i) = 1-pr(i) = 1/(1+exp(z)) 
# z = b0 + (b1 * max_cnt_pre_infu) + (b2 * max_egn_cen) + (b3 * max_cen_egn) + (b4 * max_dia_len) 
# b0: coefficient of intercept 
# b1: coefficient of max_cnt_pre_infu 
# b2: coefficient of max_egn_cen
# b3: coefficient of max_cen_egn
# b4: coefficient of max_dia_len
#---------  

# Using LASSO coefficients
beta_coef = coef(la)[abs(coef(la)[,1])>0,]  # fetch lasso regression coefficients

# Using artificial lasso coefficients
beta_artificial = c(-1, rep(0.3,4)) # 0.21

#y_sim_initial = y_sim_from_lasso_coeffs(la, X_sim_scaled, artifical=NULL)
y_sim_initial = y_sim_from_lasso_coeffs(la, X_sim_scaled, artifical=beta_artificial)
df_sim = data.frame(X_sim_scaled, y_sim_initial)


#trainUp = UpSampling(df_sim, "y_sim_initial")
trainDown = DownSampling(df_sim, "y_sim_initial")


# The true model here is the the features from the earlier lasso (25,43,82,89) model
 
################################################################################  
# Iterations
  itr = 10
  y_sim_dwn = trainDown$y_sim_initial
  dim(X_sim_scaled_dwn)
  y_sim_dwn[y_sim_dwn == 0] = -1
  X_sim_scaled_dwn = subset(trainDown, select = -y_sim_initial)
  
  pB = 10
  SS = 0.4 # threshold
  
  #-----------------
  # Group Lasso with cross validation and permutation tuning
  #-----------------
  
  #true_features_idx = grp_lasso_grp_idx # Using group_lasso_cv as reference
  true_features_idx = c(5,8,14,15) # Using lasso_cv as the true causal variables
  #true_features_idx = c(1,5,6,8,14,15) # Using lasso_cv as reference
  len_true_features = length(true_features_idx)
  
  grp_plasso_sensitivity = vector(length=itr)
  grp_plasso_fdr = vector(length=itr)
  grp_plasso_f1 = vector(length=itr)
  grp_plasso_power = vector(length=len_true_features)
  grp_plasso_df = data.frame()
  
  grp_lasso_cv_sensitivity = vector(length=itr)
  grp_lasso_cv_fdr = vector(length=itr)
  grp_lasso_cv_f1 = vector(length=itr)
  grp_lasso_cv_power = vector(length=len_true_features)
  grp_lasso_cv_df = data.frame()
  
  set.seed(25)
  for(i in 1:itr)
  {
    # For group lasso with cross validation
    grp_lasso_cv_features = grp_lasso_cv(X_sim_scaled_dwn, y_sim_dwn, 0)
    print(paste("GRP LASSO CV - Iteration#", i))
    #print("Feature Indexes: ")
    #print(paste(grp_lasso_cv_features, collapse=", "))
    grp_lasso_cv_grp_idx = select_grps(grp_lasso_cv_features)
    print("Group Indexes: ")
    print(paste(grp_lasso_cv_grp_idx, collapse=", "))
    len_feature_idx = length(grp_lasso_cv_grp_idx)
    
    mat_pred_features = matrix(cbind(grp_lasso_cv_grp_idx,rep(i,len_feature_idx)), nrow = len_feature_idx)
    grp_lasso_cv_df = rbind(grp_lasso_cv_df, mat_pred_features)
    
    grp_lasso_cv_sensitivity[i] = calc_sensitivity(true_features_idx, grp_lasso_cv_grp_idx)
    grp_lasso_cv_fdr[i] = calc_fdr(true_features_idx, grp_lasso_cv_grp_idx)
    grp_lasso_cv_power = grp_lasso_cv_power + calc_power(true_features_idx, grp_lasso_cv_grp_idx)
    
    pre_grp_lasso_cv = 1-grp_lasso_cv_fdr[i]
    rec_grp_lasso_cv = grp_lasso_cv_sensitivity[i]
    if(pre_grp_lasso_cv + rec_grp_lasso_cv > 0){
      grp_lasso_cv_f1[i] = 2*pre_grp_lasso_cv*rec_grp_lasso_cv/(pre_grp_lasso_cv+rec_grp_lasso_cv)
    }else{
      grp_lasso_cv_f1[i] = 0
    }
  }
  
  set.seed(25)
  for(i in 1:itr)
  {
    #i=1
    # For group lasso with permutation tuning
    grp_plasso_features = grp_plassob(X_sim_scaled_dwn, y_sim_dwn, pB, SS)
    print(paste("GRP PLASSO - Iteration#", i))
    #print("Feature Indexes: ")
    #print(paste(grp_plasso_features, collapse=","))
    grp_plasso_grp_idx = select_grps(grp_plasso_features)
    print("Group Indexes: ")
    print(paste(grp_plasso_grp_idx, collapse=","))
    len_feature_idx = length(grp_plasso_grp_idx)
    
    mat_pred_features = matrix(cbind(grp_plasso_grp_idx,rep(i,len_feature_idx)), nrow = len_feature_idx)
    grp_plasso_df = rbind(grp_plasso_df, mat_pred_features)
    
    grp_plasso_sensitivity[i] = calc_sensitivity(true_features_idx, grp_plasso_grp_idx)
    grp_plasso_fdr[i] = calc_fdr(true_features_idx, grp_plasso_grp_idx)
    grp_plasso_power = grp_plasso_power + calc_power(true_features_idx, grp_plasso_grp_idx)
    
    pre_grp_plasso = 1-grp_plasso_fdr[i]
    rec_grp_plasso = grp_plasso_sensitivity[i]
    if((pre_grp_plasso + rec_grp_plasso)> 0){
      grp_plasso_f1[i] = 2*pre_grp_plasso*rec_grp_plasso/(pre_grp_plasso+rec_grp_plasso)
    }else{
      grp_plasso_f1[i] = 0
    }
  }
  
  grp_lasso_cv_df%>%
    dplyr::group_by(V1)%>%
    dplyr::summarize(n())
  
  grp_plasso_df%>%
    dplyr::group_by(V1)%>%
    dplyr::summarize(n())
  
  grp_lasso_cv_stability = calc_stability(grp_lasso_cv_df, itr)
  colNames = c("Model", 
               "True_Group_Indexes", 
               "Sensitivity", 
               "FDR", 
               "F-1", 
               "Power", 
               "Stability")
  result_grp_lasso_cv = data.frame("GROUP_LASSO_CV", 
                                   paste("Grp-",true_features_idx, collapse=","),
                                   mean(grp_lasso_cv_sensitivity),
                                   mean(mean(grp_lasso_cv_fdr)),
                                   mean(grp_lasso_cv_f1),
                                   paste((grp_lasso_cv_power)/itr, collapse = ","),
                                   grp_lasso_cv_stability)
  colnames(result_grp_lasso_cv)=colNames
  
  
  grp_plasso_stability = calc_stability(grp_plasso_df, itr)
  colNames = c("Model", 
               "True_Group_Indexes", 
               "Sensitivity", 
               "FDR", 
               "F-1", 
               "Power", 
               "Stability")
  result_grp_plasso = data.frame("GROUP_PLASSO", 
                                 paste("Grp-",true_features_idx, collapse=","),
                                 mean(grp_plasso_sensitivity),
                                 mean(mean(grp_plasso_fdr)),
                                 mean(grp_plasso_f1),
                                 paste((grp_plasso_power)/itr, collapse = ","),
                                 grp_plasso_stability)
  colnames(result_grp_plasso)=colNames
  
  result_grp_lasso_cv
  result_grp_plasso
  
  
 
  
  #-----------------
  # LASSO
  #-----------------
  true_features_idx = plasso_feature_idx
  len_true_features = length(true_features_idx)
  
  plasso_sensitivity = vector(length=itr)
  plasso_fdr = vector(length=itr)
  plasso_f1 = vector(length=itr)
  plasso_power = vector(length=len_true_features)
  plasso_df = data.frame() 
  
  lasso_cv_sensitivity = vector(length=itr)
  lasso_cv_fdr = vector(length=itr)
  lasso_cv_f1 = vector(length=itr)
  lasso_cv_power = vector(length=len_true_features)
  lasso_cv_df = data.frame() 
  
  # Lasso with cross validation
  set.seed(25)
  for(i in 1:itr)
  {
    true_features_idx = plasso_feature_idx
    la_cv_sim = cv.glmnet(as.matrix(X_sim_scaled_dwn), y_sim_dwn, 
                          family="binomial", nfolds=5, 
                          alpha = 1)
    la_sim = glmnet(as.matrix(X_sim_scaled_dwn), y_sim_dwn, family="binomial", 
                lambda = la_cv_sim$lambda.min, 
                alpha = 1)
    
    lasso_cv_feature_idx = NULL
    # extract non-zero coefficients 
    for(j in 2:length(coef(la_sim))){
      if(abs(coef(la_sim)[j,1])>0){
        lasso_cv_feature_idx = c(lasso_cv_feature_idx,j-1)
      }
    }
    if(length(lasso_cv_feature_idx)==0){
      lasso_cv_feature_idx = 0
    }
    
    print(paste("LASSO - Iteration#", i))
    print("Feature Indexes: ")
    print(lasso_cv_feature_idx)
    len_feature_idx = length(lasso_cv_feature_idx)
    
    mat_pred_features = matrix(cbind(lasso_cv_feature_idx,
                                     rep(i,len_feature_idx)), nrow = len_feature_idx)
    lasso_cv_df = rbind(lasso_cv_df, mat_pred_features)
    
    lasso_cv_sensitivity[i] = calc_sensitivity(true_features_idx, lasso_cv_feature_idx)
    #if(length(lasso_cv_feature_idx) == 1 && lasso_cv_feature_idx == 0){
      #lasso_cv_fdr[i] = 1
    #}else{
      lasso_cv_fdr[i] = calc_fdr(true_features_idx, lasso_cv_feature_idx)
    #}
    
    lasso_cv_power = lasso_cv_power + calc_power(true_features_idx, lasso_cv_feature_idx)
    # Precision = 1-FDR
    # Recall = Sensitivity
    # F1 = 2*precision*recall/(precision+recall) = 2*(1-FDR)*Sensitivity/{(1-FDR)+Sensitivity}
    pre_lasso_cv = 1-lasso_cv_fdr[i]
    rec_lasso_cv = lasso_cv_sensitivity[i]
    if(pre_lasso_cv + rec_lasso_cv > 0){
      lasso_cv_f1[i] = 2*pre_lasso_cv*rec_lasso_cv/(pre_lasso_cv+rec_lasso_cv)
    }else if(rec_lasso_cv == FALSE){
      lasso_cv_f1[i] = 0
    }else{
      lasso_cv_f1[i] = 0
    }
  }
  
  lasso_cv_stability = calc_stability(lasso_cv_df, itr)
  colNames = c("Model", 
               "True_Feature_Indexes", 
               "Sensitivity", 
               "FDR", 
               "F-1", 
               "Power", 
               "Stability")
  result_lasso_cv = data.frame("LASSO_CV", 
                               paste(plasso_feature_idx,collapse=","), 
                               mean(lasso_cv_sensitivity),
                               mean(lasso_cv_fdr),
                               mean(lasso_cv_f1),
                               paste((lasso_cv_power)/itr, collapse = ","),
                               lasso_cv_stability)
  colnames(result_lasso_cv)=colNames
  
  
  # Lasso with permutation tuning
  set.seed(25)
  for(i in 1:itr)
  {
    true_features_idx = plasso_feature_idx
    plasso_features_idx = plassob(X_sim_scaled_dwn, y_sim_dwn, pB, SS)
    print(paste("PLASSO - Iteration#", i))
    print("Feature Indexes: ")
    print(plasso_features_idx)
    len_feature_idx = length(plasso_features_idx)
    
    mat_pred_features = matrix(cbind(plasso_features_idx,rep(i,len_feature_idx)), 
                               nrow = len_feature_idx)
    plasso_df = rbind(plasso_df, mat_pred_features)
    
    
    plasso_sensitivity[i] = calc_sensitivity(true_features_idx, plasso_features_idx)
    plasso_fdr[i] = calc_fdr(true_features_idx, plasso_features_idx)
    plasso_power = plasso_power + calc_power(true_features_idx, plasso_features_idx)
    # Precision = 1-FDR
    # Recall = Sensitivity
    # F1 = 2*precision*recall/(precision+recall) = 2*(1-FDR)*Sensitivity/{(1-FDR)+Sensitivity}
    pre_plasso = 1-plasso_fdr[i]
    rec_plasso = plasso_sensitivity[i]
    if(pre_plasso + rec_plasso > 0){
      plasso_f1[i] = 2*pre_plasso*rec_plasso/(pre_plasso+rec_plasso)
    }else{
      plasso_f1[i] = 0
    }
  }
  plasso_stability = calc_stability(plasso_df, itr)
  colNames = c("Model", 
               "True_Feature_Indexes", 
               "Sensitivity", 
               "FDR", 
               "F-1", 
               "Power", 
               "Stability")
  result_plasso = data.frame("PLASSO", 
                             paste(plasso_feature_idx,collapse=","), 
                             mean(plasso_sensitivity),
                             mean(plasso_fdr),
                             mean(plasso_f1),
                             paste((plasso_power)/itr, collapse = ","),
                             plasso_stability)
  colnames(result_plasso)=colNames
  
 
  # Display the frequency of the selected features across the 10 iterations
  lasso_cv_feature_freq = lasso_cv_df%>%
    dplyr::group_by(V1)%>%
    dplyr::summarize(freq=n())

  plasso_feature_freq = plasso_df%>%
    dplyr::group_by(V1)%>%
    dplyr::summarize(freq=n())
  
  lasso_cv_feature_freq
  plasso_feature_freq
  
  
 # Display Lasso_CV and Plasso results from simulation study
  result_lasso_cv
  result_plasso
  
 
  #-----------------
  #Exclusive Lasso with Cross Validation
  #-----------------
  exclsv_lasso_sensitivity = vector(length=itr)
  exclsv_lasso_fdr = vector(length=itr)
  exclsv_lasso_f1 = vector(length=itr)
  exclsv_lasso_power = vector(length=len_true_features)
  exclsv_lasso_df = data.frame()
  
  for(i in 1:itr)
  {
    exclsv_lasso_features_idx = exclsv_lasso(X_sim_scaled_dwn, y_sim_dwn, nlambda = 50)
    print(paste("EXCLUSIVE LASSO - Iteration#", i))
    print("Feature Indexes: ")
    print(exclsv_lasso_features_idx)
    len_feature_idx = length(exclsv_lasso_features_idx)
    
    mat_pred_features = matrix(cbind(exclsv_lasso_features_idx,rep(i,len_feature_idx)), 
                               nrow = len_feature_idx)
    exclsv_lasso_df = rbind(exclsv_lasso_df, mat_pred_features)
    
    exclsv_lasso_sensitivity[i] = calc_sensitivity(true_features_idx, exclsv_lasso_features_idx)
    exclsv_lasso_fdr[i] = calc_fdr(true_features_idx, exclsv_lasso_features_idx)
    exclsv_lasso_power = exclsv_lasso_power + calc_power(true_features_idx, 
                                                         exclsv_lasso_features_idx)
    
    exclsv_pre_lasso = 1-exclsv_lasso_fdr[i]
    exclsv_rec_lasso = exclsv_lasso_sensitivity[i]
    if(exclsv_pre_lasso + exclsv_rec_lasso > 0){
      exclsv_lasso_f1[i] = 2*exclsv_pre_lasso*exclsv_rec_lasso/
        (exclsv_pre_lasso+exclsv_rec_lasso)
    }else{
      exclsv_lasso_f1[i] = 0
    }
    
  }
  
  exclsv_lasso_stability = calc_stability(exclsv_lasso_df, itr)
  colNames = c("Model", 
               "True_Feature_Indexes",
               "Sensitivity", 
               "FDR", 
               "F-1", 
               "Power", 
               "Stability")
  result_exclusv_lasso = data.frame("EXCLUSIVE_LASSO", 
                                    paste(true_features_idx, collapse=","), 
                                    mean(exclsv_lasso_sensitivity),               
                                    mean(mean(exclsv_lasso_fdr)),                                 
                                    mean(exclsv_lasso_f1),
                                    paste((exclsv_lasso_power)/itr, collapse = ","),
                                    exclsv_lasso_stability)
  colnames(result_exclusv_lasso)=colNames
  result_exclusv_lasso
 
  
  
  col_cor = c(25,43,82,89)
  X_corr = X_scaled[,col_cor]
  pairs(X_corr)
  
  
  cor_matrix = data.frame(cor(X_corr))
  ggcorrplot(cor_matrix)+ggtitle("Heat Map for the true causal variables")
  
  
}

#############     END OF MAIN     ####################

#######################    FUNCTIONS     ########################

### Generate more samples and then down size

# DownSampling the majority class
DownSampling = function(df,y){
  #df=df_sim
  #y="y_sim_initial"
  tbl = sort(table(df[y]))
  
  y_majority_val = names(tbl)[2]
  y_minority_count = as.numeric(tbl)[1]
  y_majority_count = as.numeric(tbl)[2]
  
  df_majority = df[df[y]==0,]
  df_majority_downsized = sample_n(df_majority, y_minority_count, replace=FALSE)
  
  return(rbind(df[df[y]==1,],df_majority_downsized))
}

# UpSampling the minority class
UpSampling = function(df,y){
  tbl = sort(table(df[y]))
  y_minority_val = names(tbl)[1]
  y_minority_count = as.numeric(tbl)[1]
  y_majority_count = as.numeric(tbl)[2]
  
  num_rows_toadd = y_majority_count-y_minority_count
  df_minority = df[df[y]==1,]
  df_minority_added = sample_n(df_minority, num_rows_toadd, replace=TRUE)
  
  return(rbind(df,df_minority_added))
}

# Simulate Y from LASSO coefficients
y_sim_from_lasso_coeffs = function(la, X_sim_scaled, artifical){
  n = dim(X_sim_scaled)[1]
  
  #artifical = c(-1, 0.15, 0.15, 0.15, 0.15)
  
  if(is.null(artifical)){
    beta_coef = coef(la)[abs(coef(la)[,1])>0,]  # fetch lasso regression coefficients
  } else{
    beta_coef = artifical  # fetch artificial coefficients
  }
  
  
  #View(X_sim_scaled)
  
  df = cbind(X_sim_scaled[,25],X_sim_scaled[,43],X_sim_scaled[,82],X_sim_scaled[,89])
  z = beta_coef[1] + beta_coef[2]*df[,1] + beta_coef[3]*df[,2] + 
    beta_coef[4]*df[,3] + beta_coef[5]*df[,4]
  
  prob = exp(z)/(1+exp(z))
  
  set.seed(81)
  y_sim = rbern(n, prob)
  return(y_sim)
}

# Calculate Sensitivity
calc_sensitivity = function(true_features_idx, predicted_features_idx){
  sentvt = length(intersect(true_features_idx, predicted_features_idx))/length(true_features_idx)
  return(sentvt)
}

# Calculate False Discovery Rate FDR
calc_fdr = function(true_features_idx, predicted_features_idx){
  fdr = length(setdiff(predicted_features_idx, true_features_idx))/length(predicted_features_idx)
  return(fdr)
}

# Calculate Power for each causal variable
calc_power = function(true_features_idx, predicted_features_idx){
  #predicted_features_idx = grp_lasso_cv_features_idx
  
  true_features_idx = sort(true_features_idx)
  
  true_len = length(true_features_idx)
  output = rep(0, true_len)
  
  for(i in 1:true_len){
    val = length(intersect(true_features_idx[i], predicted_features_idx))
    if(val>0){
      output[i]=1
    }
  }
  return(output)
}

# Calculate Stability
calc_stability = function(df, itr){
  #df = plasso_df
  colnames(df) = c("feature_idx","iter")
  Jaccard_index = c()
  for(i in 1:(itr-1)){
    a = df[df$iter == i,1]
    #print(a)
    
    for(j in (i+1):itr){
      b = df[df$iter == j,1]
      #print(b)
      Jaccard_index = c(Jaccard_index,length(intersect(a,b))/length(union(a,b)))
    }
  }
  Jaccard_index
  return(mean(Jaccard_index))
}

# plasso for a binary phenotype:  
plassob = function(X, y, pB, SS){
  #X = X_sim_scaled_dwn
  #y = y_sim_dwn
  n = nrow(X)
  p = ncol(X)
  Spi = NULL
  
  
  # Calculate selection stability by running the logistic regression for B permutations
  for(k in 1:pB)  # pB: # of permutations
  {
    X_ko1 = X[c(sample(nrow(X))), ]  # permutation copy
    
    # fitting logistic regression with penalized maximum likelihood
    a = glmnet(cbind(X, X_ko1), y, family = "binomial", alpha = 1, maxit = 1e+07, thresh = 1e-05)  
    lasso = as.matrix(a$beta)
    ii = 1
    while (max(abs(lasso[(p+1):(2*p),ii])) == 0 & ii < dim(lasso)[2]){
      ii = ii+1 
    }
    selected_lasso = which(abs(lasso[,ii-1]) > 0)
    Spi = c(Spi, selected_lasso)
  }
  freq = tabulate(Spi)/pB
  out = which(freq > SS | freq == SS)
  
  if(length(out)>0){
    return(out)
  }
  else{
    return(0)
  }
  
} 

## GROUP LASSO using permutation tuning for a binary phenotype:  
grp_plassob = function(X, y, pB, SS){
  #X = X_sim_scaled_dwn
  #y = y_sim_dwn
  n = nrow(X)
  p = ncol(X)
  Spi = NULL
  # group index for X variables
  v.group = c(1,	rep(2,6),	rep(3,6),	rep(4,6),	rep(5,6),	rep(6,6),	rep(7,6),	rep(8,6),	
              rep(9,7),	rep(10,7),	rep(11,6),	rep(12,6),	rep(13,7),	rep(14,6),	rep(15,7),	
              16,	rep(17,6),	rep(18,6),	rep(19,6),	rep(20,6),	rep(21,6),	rep(22,6),	
              rep(23,6),	rep(24,7),	rep(25,7),	rep(26,6),	rep(27,6),	rep(28,7),	rep(29,6),	
              rep(30,7))

  
  # Calculate selection stability by running the logistic regression for B permutations
  for(k in 1:pB)  # pB: # of permutations
  {
    #set.seed(25)
    X_ko1 = X[c(sample(nrow(X))), ]  # permutation copy
    
    mat_X = as.matrix(cbind(X, X_ko1))
    a = gglasso(x=mat_X, y, loss="logit", intercept=T, group = v.group, eps=10^(-4))  
    
    # Note that a default of 100 different lambda values from max to min are applied
    
    grp_lasso = as.matrix(a$beta)
    #View(grp_lasso)
    
    ii = 1
    while (max(abs(grp_lasso[(p+1):(2*p),ii])) == 0 & ii < dim(grp_lasso)[2]){
      ii = ii+1 
    }
    selected_grp_lasso = which(abs(grp_lasso[,ii-1]) > 0)
    Spi = c(Spi, selected_grp_lasso)
  }
  freq = tabulate(Spi)/pB
  out = which(freq > SS | freq == SS)
  
  if(length(out)>0){
    return(out)
  }
  else{
    return(0)
  }
} 

# To find the groups for the selected features
select_grps = function(features_idx){
  select_grp = NULL
  v.group = c(1,	rep(2,6),	rep(3,6),	rep(4,6),	rep(5,6),	rep(6,6),	rep(7,6),	rep(8,6),	
              rep(9,7),	rep(10,7),	rep(11,6),	rep(12,6),	rep(13,7),	rep(14,6),	rep(15,7),	
              16,	rep(17,6),	rep(18,6),	rep(19,6),	rep(20,6),	rep(21,6),	rep(22,6),	
              rep(23,6),	rep(24,7),	rep(25,7),	rep(26,6),	rep(27,6),	rep(28,7),	rep(29,6),	
              rep(30,7))
  
  if(length(features_idx)>0){
    for(i in 1:length(features_idx)){
      idx = features_idx[i]
      select_grp = c(select_grp, v.group[idx])
    }
    return(unique(select_grp))
  }
  else{
    return(0)
  }
}

## GROUP LASSO using cross validation:  
grp_lasso_cv = function(X, y, lambda_val){
  #X = X_sim_scaled_dwn
  #y = y_sim_dwn
  #View(X_scaled)
 # X = X_scaled
  #y = y_obsvd
  Spi = NULL
  
  # group index for X variables
  v.group = c(1, rep(2,6),	rep(3,6),	rep(4,6),	rep(5,6),	rep(6,6),	rep(7,6),	
              rep(8,6),	rep(9,7),	rep(10,7),	rep(11,6),	rep(12,6),	
              rep(13,7),	rep(14,6),	rep(15,7))
  
  # Cross Validation
  grp_cv = cv.gglasso(as.matrix(X), y, loss="logit", group = v.group, intercept=T, 
                      eps=10^(-5), nfolds=5)
  
  
  #paste(round(grp_cv$lambda.min,6), round(grp_cv$lambda.1se,6))
  
  lambda_use = grp_cv$lambda.min-lambda_val
  grp_lasso = gglasso(as.matrix(X), y, loss="logit", intercept=F, 
                      group = v.group, lambda=lambda_use, 
                      eps=10^(-4)) 
  #round(coef(grp_lasso)[abs(coef(grp_lasso)[,1]) > 0,],6) # extract non-zero coefficients
  
  
  for(i in 2:length(coef(grp_lasso))){
    if(abs(coef(grp_lasso)[i]) > 0){
      Spi = c(Spi,i-1)
    }
  }
  
  if(length(Spi)>0){
    return(Spi)
  }
  else{
    return(0)
  }
} 

## Exclusive LASSO for a binary phenotype:  
exclsv_lasso = function(X, y, nlambda_val){
  #X = X_sim_scaled_dwn
  #y = y_sim_dwn
  nlambda_val = 50
  
  y[y == -1] = 0
  n = nrow(X)
  p = ncol(X)
  Spi = NULL
  # group index for X variables
  v.group = c(1, rep(2,6), rep(3,6),	rep(4,6),	rep(5,6),	rep(6,6),	rep(7,6),	rep(8,6),	
              rep(9,7),	rep(10,7), rep(11,6), rep(12,6), rep(13,7), rep(14,6), rep(15,7))
  
  
  mat_X = as.matrix(X)
  
  # Cross Validation
  exclusv_cv = cv.exclusive_lasso(mat_X, y, groups = v.group, intercept = F, 
                                  family="binomial",nfolds=5, standardize = FALSE, 
                                  algorithm = "pg", nlambda=nlambda_val)
  
  
  a = exclusive_lasso(mat_X, y,lambda = exclusv_cv$lambda.1se,
                      groups = v.group, family="binomial", intercept = FALSE,
                      standardize = FALSE, algorithm = "pg")
  
  for(i in 1:length(a$coef[,1])){
    if(abs(a$coef[i,1]) > 0){
      Spi = c(Spi,i)
    }
  }
  
  if(length(Spi)>0){
    return(Spi)
  }
  else{
    return(0)
  }
} 


# Q_Q plot
qq_function = function(ds, n, main_text, ylab_text){
  ds = sort(ds)
  marginal_quantiles = vector(length=n)
  
  for(j in 1:n){
    prob_val = (j-0.5)/n
    marginal_quantiles[j] = qnorm(prob_val, mean=0, sd=1)
  }
  
  plot(marginal_quantiles, ds, main = main_text, xlab = "Quantiles", 
       ylab = ylab_text, pch = 21, bg = "green", col = "blue")
  
  rq_coeff_num = 0
  rq_coeff_denom_1 = 0
  rq_coeff_denom_2 = 0
  
  for(i in 1:n)
  {
    rq_coeff_num = rq_coeff_num + (ds[i]-mean(ds)) * marginal_quantiles[i]
    rq_coeff_denom_1 = rq_coeff_denom_1 + ((ds[i]-mean(ds))^2)  
    rq_coeff_denom_2 = rq_coeff_denom_2 + (marginal_quantiles[i]^2)
  }
  
  return(rq_coeff_num/(sqrt(rq_coeff_denom_1)*sqrt(rq_coeff_denom_2)))
    
}

####################################################################

summary(ds_os)
OS_median = 20.3


ds = read.csv("C:\\Users\\sengu\\Dropbox\\PC\\Desktop\\SFSU_MasterThesis\\GeneratedData\\OS_assortativity.csv", 
              header=TRUE, sep=",")
View(ds)
summary(ds)


#----  Histogram plots  ----#
ds_treated = ds[which(ds$OS_mon >= OS_median),]
hist(ds_treated$assortativity, xlab = "Higher OS_mon patients", 
     ylab = "Assortativity", main = "Assortativity of patients with higher OS_mon",
     col = c("pink"))
summary(ds_treated)


ds_untreated = ds[which(ds$OS_mon < OS_median),]
hist(ds_untreated$assortativity, xlab = "Lower OS_mon patients", 
     ylab = "Assortativity", main = "Assortativity of patients with lower OS_mon", 
     col = c("yellow"))
summary(ds_untreated)

n = dim(ds)[1]
n_treated = dim(ds_treated)[1]
n_untreated = dim(ds_untreated)[1]
n
n_treated
n_untreated


#####  Checking for Normalancy
qq_function(sort(ds_treated$assortativity), n_treated,"Q-Q plot of Assortativity for higher OS_mon patients", "Assortativity")
qq_function(sort(ds_untreated$assortativity), n_untreated,"Q-Q plot of Assortativity for lower OS_mon patients", "Assortativity")



ds_cluster = read.csv("C:\\Users\\sengu\\Dropbox\\PC\\Desktop\\SFSU_MasterThesis\\GeneratedData\\Cluster_info.csv", 
                      header=TRUE, sep=",")
#View(ds_cluster)
#summary(ds_cluster)





# For assortativity NA

ds_assort_NA = read.csv("C:\\Users\\sengu\\Dropbox\\PC\\Desktop\\SFSU_MasterThesis\\GeneratedData\\Prob_NA_per_Patient.csv", 
                        header=TRUE, sep=",")
#View(ds_assort_NA)
summary(ds_assort_NA)

ds_assort_NA_treated = ds_assort_NA[which(ds_assort_NA$Patient_Survival_Time == "Increased"),]

ds_assort_NA_untreated = ds_assort_NA[which(ds_assort_NA$Patient_Survival_Time != "Increased"),]

summary(ds_assort_NA_treated)
summary(ds_assort_NA_untreated)

hist(ds_assort_NA_treated$prob_NA_per_patient, xlab = "Higher OS_mon patients", 
     ylab = "Assortativity", main = "Patients with higher OS_mon and Assortativity = NA",
     col = c("pink"))

hist(ds_assort_NA_untreated$prob_NA_per_patient, xlab = "Lower OS_mon patients", 
     ylab = "Assortativity", main = "Patients with lower OS_mon and Assortativity = NA",
     col = c("yellow"), ylim=c(0,25))

####################################################################################

# For positive assortativity

ds_assort_pos = read.csv("C:\\Users\\sengu\\Dropbox\\PC\\Desktop\\SFSU_MasterThesis\\GeneratedData\\Positive_Assort_Prob_per_Patient.csv", 
                         header=TRUE, sep=",")
#View(ds_assort_pos)
summary(ds_assort_pos)

ds_assort_pos_treated = ds_assort_pos[which(ds_assort_pos$Patient_Survival_Time == "Increased"),]

ds_assort_pos_untreated = ds_assort_pos[which(ds_assort_pos$Patient_Survival_Time != "Increased"),]

summary(ds_assort_pos_treated)
summary(ds_assort_pos_untreated)

hist(ds_assort_pos_treated$prob_pos_per_patient, xlab = "Higher OS_mon patients", 
     ylab = "Assortativity", main = "Patients with higher OS_mon and positive Assortativity",
     col = c("green"))

hist(ds_assort_pos_untreated$prob_pos_per_patient, xlab = "Lower OS_mon patients", 
     ylab = "Assortativity", main = "Patients with lower OS_mon and positive Assortativity",
     col = c("blue"), ylim=c(0,20), xlim=c(0.60,0.9))

#########################################################################################





#--------- 
#-----------  Data Simulation and Hypothesis Testing --------------#
#---------

# y is Bernoulli rvs
# y = 1 with probability pr(i)
# y = 0 with probability 1-pr(i)
# P(Y|X=i) = pr(i) = 1/(1+exp(-(z)))
# 1-pr(i) = exp(-(z))/(1+exp(-(z))) 
# z = b0 + (b1 * max_cnt_pre_infu) + (b2 * max_deg_avg) + (b3 * max_egn_cen) + (b4 * max_cen_egn)
# b0: coefficient of intercept = -1.061324
# b1: coefficient of max_cnt_pre_infu = 0.205565
# b2: coefficient of max_deg_avg = 0.000069
# b3: coefficient of max_egn_cen = 0.148353
# b4: coefficient of max_cen_egn = 0.100524
#---------  
set.seed(81)
n = 1000

# For LASSO model coefficients
beta_coef = coef(la)[abs(coef(la)[,1])>0,]  # fetch lasso regression coefficients


# Data Simulation 

df = normal_distr(n) # Assuming normal distribution of the predictors
#df = expo_distr(n) # Assuming exponential distribution of the predictors
#df = gamma_distr(n) # Assuming gamma distribution of the predictors

z = beta_coef[1] + beta_coef[2]*df[,1] + beta_coef[3]*df[,2] + beta_coef[4]*df[,3] + 
  beta_coef[5]*df[,4]

prob = 1/(1+exp(-z))

set.seed(81)
y = rbinom(n, 1, prob)

model_1 = glm(y~., data=data.frame(df,y), family="binomial")
summary(model_1)




#######################################################################

treated_patients = ds_os%>%
  dplyr::filter(OS_mon>=20.3)%>%
  dplyr::select(Patient_id)

cluster_level_treated_patients = all_data%>%
  dplyr::filter(Patient_id %in% treated_patients$Patient_id)%>%
  dplyr::group_by(Patient_id)%>%
  dplyr::summarise(max_clusters = max(membership),
    max_node_count = max(node_count))

cluster_level_untreated_patients = all_data%>%
  dplyr::filter(!(Patient_id %in% treated_patients$Patient_id))%>%
  dplyr::group_by(Patient_id)%>%
  dplyr::summarise(max_clusters = max(membership),
                   max_node_count = max(node_count))

ds1 = data.frame(cluster_level_treated_patients, 
                                  rep("Longer Survival",dim(cluster_level_treated_patients)[1]))
ds2 = data.frame(cluster_level_untreated_patients, 
           rep("Shorter Survival",dim(cluster_level_untreated_patients)[1]))

colnames(ds1)[4]="class"
colnames(ds2)[4]="class"

aggreg_cluster_level = rbind(ds1,ds2)



ggplot(aggreg_cluster_level, aes(x = class, y = max_clusters, 
                                 group=class,fill = class)) + 
  geom_boxplot()+scale_fill_manual(values=c("#00C8FF","#FF80CA"))+ 
  ggtitle("# of Clusters in Longer and Shorter Survival Groups") +
  ylab("# of Clusters")+
  theme(axis.text=element_text(size=rel(1.2),face="bold"),plot.title = element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)), #change legend title font size
        legend.text = element_text(size=rel(1.5)))


ggplot(aggreg_cluster_level, aes(x = class, y = max_node_count, fill = class, group=class)) + 
  geom_boxplot()+
  scale_fill_manual(values=c("#00FFC6","#FF5D00"))+ 
  ggtitle("# of Nodes in Longer and Shorter Survival Groups") +
  ylab("# of Nodes")+
  theme(axis.text=element_text(size=rel(1.2),face="bold"),plot.title = element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)), #change legend title font size
        legend.text = element_text(size=rel(1.5)))

+
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white")
  

# Simulate number of samples
n = 200


# Node_Count influence on other attributes
val_nd_cnt = ds_data%>%
  dplyr::filter(node_count == 2)%>%
  dplyr::group_by(node_count, assortativity, deg, transitivity, 
                  AA_length, diam_length, deg_avg, edge_density, centr_degree,
                  centr_clo, eigen_centrality, centr_eigen)%>%
  dplyr::select(node_count, assortativity, deg, transitivity, 
                AA_length, diam_length, deg_avg, edge_density, centr_degree,
                centr_clo, eigen_centrality, centr_eigen)
summary(val_nd_cnt)


hist(ds_data$assortativity[ds_data$Patient_id==10025011361], freq=FALSE, col="lightyellow")
# Assortativity probability at patient/sample level
Patient_level_assort = ds_data%>%
  dplyr::summarise(prob_NA = sum(is.na(assortativity))/n(), 
                   prob_neg_1 = sum(assortativity == -1, na.rm = T)/n(), 
                   prob_num = sum(assortativity > -1, na.rm = T)/n())

sum(ds_data$assortativity>0.5, na.rm=TRUE)/sum(ds_data$assortativity>= -1 , na.rm=TRUE)
summary(Patient_level_assort)
View(Patient_level_assort)

###################################################




#----------------------------------------------------------#  
unique(ds_data$Patient_id[ds_data$OS_mon < 20.3])
unique(ds_data$Patient_id[ds_data$OS_mon >= 20.3])
par(mfrow=c(2,1))
hist(ds_data$node_count[ds_data$Patient_id==10562011565], col="lightyellow")
hist(ds_data$node_count[ds_data$Patient_id==2000206980], col="lightyellow")

ind_sample_mem_nodes = ds_data%>%
  dplyr::filter(Patient_id==10025011361)%>%
  dplyr::select(node_count)

length(which(ind_sample_mem_nodes >= 2))  
#----------------------------------------------------------#

