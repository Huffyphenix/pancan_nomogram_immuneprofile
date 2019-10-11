####################################
# TCGA pancan survival analysis
# include TIL, expression data, and clinical data

source(file.path("/home/huff/project/github/pancan_nomogram_immuneprofile/3.survival_function.R"))
# path --------------------------------------------------------------------

basic_path <- "/home/huff/project"
data_path <- file.path(basic_path,"TCGA_nomogram/data")
res_path <- file.path(basic_path,"TCGA_nomogram/result/muli_uni_survival")


# load data ---------------------------------------------------------------

combined_clinical_TIL_exp_data <-
  readr::read_rds(file.path(data_path,"TCGA_combined_clinical_TIL_exp_data.rds.gz"))

# calculation -------------------------------------------------------------
# OS, univariable
combined_clinical_TIL_exp_data %>%
  # head(1) %>%
  dplyr::mutate(combine_data = purrr::pmap(list(PFS,OS_stage,Infiltration,exp_filter),fn_data_process,type="OS")) %>%
  dplyr::select(cancer_types, combine_data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(uni_OS = purrr::map2(data,features, fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> OS_uni_res
OS_uni_res %>%
  readr::write_tsv(file.path(res_path,"OS_uni_res.tsv"))


# OS, multi-variable
OS_uni_res %>%
  dplyr::filter(coxp<=0.05 | kmp <= 0.05) %>%
  dplyr::select(cancer_types,features) %>% 
  tidyr::nest(-cancer_types,.key="sig_features") -> OS_uni_res.sig

OS_uni_res %>%
  dplyr::filter(coxp<=0.05 | kmp <= 0.05) %>%
  readr::write_tsv(file.path(res_path,"OS_uni_res.sig.tsv"))

combined_clinical_TIL_exp_data %>%
  # head(1) %>%
  dplyr::mutate(combine_data = purrr::pmap(list(PFS,OS_stage,Infiltration,exp_filter),fn_data_process,type="OS")) %>%
  dplyr::select(cancer_types, combine_data) %>%
  tidyr::unnest() %>%
  tidyr::unnest() %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(combine_data.multi = purrr::map(data,.f=function(.x){
    .x %>%
      tidyr::spread(key="features",value="group")
  })) %>%
  dplyr::select(-data) %>%
  dplyr::inner_join(OS_uni_res.sig,by="cancer_types") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(combine_data.multi,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-combine_data.multi,-sig_features) %>%
  tidyr::unnest() -> OS_multi_res
OS_multi_res %>%
  readr::write_tsv(file.path(res_path,"OS_multi_res.tsv"))

OS_multi_res %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<10) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-10),-10,hr_l),hr_h=ifelse(hr_h>10,10,hr_h)) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(filename = paste("2.OS.multi-variable.cox",cancer_types,sep="_")) %>%
  dplyr::mutate(title = paste("Multi-variable cox model (OS)",cancer_types,sep=";")) %>%
  dplyr::mutate(res = purrr::pmap(list(data,filename,title),fn_cox_plot_2,dir=file.path(res_path,"plot"),w=6,h=6))

# PFS, univariable
combined_clinical_TIL_exp_data %>%
  # head(1) %>%
  dplyr::mutate(combine_data = purrr::pmap(list(PFS,OS_stage,Infiltration,exp_filter),fn_data_process,type="PFS")) %>%
  dplyr::select(cancer_types, combine_data) %>%
  tidyr::unnest() %>%
  dplyr::mutate(uni_OS = purrr::map2(data,features, fn_survival_test)) %>%
  dplyr::select(-data) %>%
  tidyr::unnest() -> PFS_uni_res
PFS_uni_res %>%
  readr::write_tsv(file.path(res_path,"PFS_uni_res.tsv"))

# OS, multi-variable
PFS_uni_res %>%
  dplyr::filter(coxp<=0.05 | kmp <= 0.05) %>%
  dplyr::select(cancer_types,features) %>% 
  tidyr::nest(-cancer_types,.key="sig_features") -> PFS_uni_res.sig

PFS_uni_res %>%
  dplyr::filter(coxp<=0.05 | kmp <= 0.05) %>%
  readr::write_tsv(file.path(res_path,"PFS_uni_res.sig.tsv"))

combined_clinical_TIL_exp_data %>%
  # head(1) %>%
  dplyr::mutate(combine_data = purrr::pmap(list(PFS,OS_stage,Infiltration,exp_filter),fn_data_process,type="PFS")) %>%
  dplyr::select(cancer_types, combine_data) %>%
  tidyr::unnest() %>%
  tidyr::unnest() %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(combine_data.multi = purrr::map(data,.f=function(.x){
    .x %>%
      tidyr::spread(key="features",value="group")
  })) %>%
  dplyr::select(-data) %>%
  dplyr::inner_join(PFS_uni_res.sig,by="cancer_types") %>%
  dplyr::mutate(surv_res.multi = purrr::map2(combine_data.multi,sig_features,fn_survival_test.multiCox)) %>%
  dplyr::select(-combine_data.multi,-sig_features) %>%
  tidyr::unnest() -> PFS_multi_res

PFS_multi_res%>%
  readr::write_tsv(file.path(res_path,"PFS_multi_res.tsv"))

PFS_multi_res %>%
  dplyr::mutate(cox_sig = ifelse(coxp<0.1,"1yes","2no")) %>%
  dplyr::mutate(hr=log2(hr)+1,hr_l=log2(hr_l)+1,hr_h=log2(hr_h)+1) %>%
  dplyr::filter(abs(hr)<10) %>%
  dplyr::mutate(hr_l=ifelse(hr_l<(-10),-10,hr_l),hr_h=ifelse(hr_h>10,10,hr_h)) %>%
  tidyr::nest(-cancer_types) %>%
  dplyr::mutate(filename = paste("2.PFS.multi-variable.cox",cancer_types,sep="_")) %>%
  dplyr::mutate(title = paste("Multi-variable cox model (PFS)",cancer_types,sep=";")) %>%
  dplyr::mutate(res = purrr::pmap(list(data,filename,title),fn_cox_plot_2,dir=file.path(res_path,"plot"),w=6,h=6))
