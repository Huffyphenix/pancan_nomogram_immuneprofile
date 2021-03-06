library(magrittr)
library(survival)
library(ggplot2)
# function to do survival analysis ----------------------------------------
my_theme <-   theme(
  panel.background = element_rect(fill = "white",colour = "black"),
  panel.grid.major=element_line(colour=NA),
  axis.text.y = element_text(size = 10,colour = "black"),
  axis.text.x = element_text(size = 10,colour = "black"),
  # legend.position = "none",
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.background = element_blank(),
  legend.key = element_rect(fill = "white", colour = "black"),
  plot.title = element_text(size = 20),
  axis.text = element_text(colour = "black"),
  strip.background = element_rect(fill = "white",colour = "black"),
  strip.text = element_text(size = 10),
  text = element_text(color = "black")
)
# data process
fn_data_process <- function(PFS,OS_stage,Infiltration,exp_filter,type){
  Infiltration %>%
    # tibble::rowid_to_column(var = "rowid") %>%
    tidyr::gather(-barcode,key="features", value="value") %>%
    # dplyr::filter(substr(barcode,14,15)=="01") %>%
    # dplyr::mutate(barcode = substr(barcode,1,12)) %>%
    # dplyr::group_by(barcode,features) %>%
    # dplyr::filter(rowid == min(rowid)) %>%
    # dplyr::ungroup() %>%
    dplyr::select(barcode,features,value) %>%
    dplyr::group_by(features) %>%
    dplyr::mutate(group = ifelse(value > quantile(value,0.5,na.rm=T),"2_high","1_low")) %>%
    dplyr::ungroup() -> TIL_gather
  exp_filter %>%
    tidyr::gather(-symbol, key="barcode", value="value") %>%
    # dplyr::filter(substr(barcode,14,15)=="01") %>%
    # dplyr::mutate(barcode = substr(barcode,1,12)) %>%
    # dplyr::group_by(barcode,symbol) %>%
    # dplyr::mutate(value = max(value)) %>%
    # unique() %>%
    # dplyr::select(-entrez_id) %>%
    dplyr::rename("features"="symbol") %>%
    dplyr::select(barcode,features,value) %>%
    dplyr::group_by(features) %>%
    dplyr::mutate(group = ifelse(value > quantile(value,0.5,na.rm=T),"2_high","1_low")) %>%
    dplyr::ungroup()  -> exp_gather
  OS_stage %>%
    dplyr::select(barcode,Age,Stage) %>%
    dplyr::mutate(Age=as.numeric(Age)) %>%
    tidyr::gather(-barcode,key="features",value="value") %>%
    tidyr::nest(-features) %>%
    dplyr::mutate(group = purrr::map2(features,data,.f=function(.x,.y){
      # print(.x)
      if (.x=="Age") {
        .y %>%
          dplyr::mutate(value = as.numeric(value)) %>%
          dplyr::mutate(group = ifelse(value > quantile(value,0.5,na.rm=T),"2_high","1_low")) %>%
          dplyr::select(-value)
      } else{
        .y %>%
          dplyr::mutate(group = value) %>%
          dplyr::select(-value)
      }
    })) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() -> .stage
  TIL_gather %>%
    rbind(exp_gather) %>%
    dplyr::select(-value) %>%
    rbind(.stage) %>%
    dplyr::inner_join(OS_stage %>% dplyr::select(barcode,OS,Status),by="barcode") %>%
    dplyr::inner_join(PFS, by="barcode") -> .combine_data
  if (type == "OS") {
    .combine_data %>%
      dplyr::rename("time"="OS","status"="Status") %>%
      dplyr::mutate(features = gsub("-",".",features)) %>%
      dplyr::mutate(time = as.numeric(time), status = as.numeric(status)) %>%
      dplyr::filter(!is.na(group)) %>%
      tidyr::nest(-features) 
  }else{
    .combine_data %>%
      dplyr::rename("time"="PFS.time","status"="PFS")%>%
      dplyr::mutate(features = gsub("-",".",features)) %>% 
      dplyr::filter(!is.na(group)) %>%
      tidyr::nest(-features)
  }
}

fn_filter_TIL_marker <- function(sig_features,TIL_markers){
  sig_features %>%
    dplyr::filter(features %in% unique(TIL_markers$Cell_type)) %>%
    .$features -> sig_cells
  TIL_markers %>%
    dplyr::filter(Cell_type %in% sig_cells) %>%
    tidyr::nest(-Cell_type) %>%
    dplyr::mutate(narrow_markers = purrr::map(data,.f=function(.x){
      .x %>%
        dplyr::filter(markers %in% sig_features$features)
    })) %>%
    dplyr::select(Cell_type,narrow_markers) %>%
    tidyr::unnest() -> markers_overlapped
  sig_features %>%
    dplyr::filter(!features %in% markers_overlapped$markers)
}
# 3.2.function to do cox survival analysis  ----
# 3.2.1.univariable cox analysis ---------
fn_survival_test <- function(data,feature,cancer_types){
  print(paste(feature,cancer_types,sep="-"))
  if(length(unique(data$group))>=2 && nrow(data)>=10){
    .cox <- survival::coxph(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)
    summary(.cox) -> .z
    
    # KM pvalue
    kmp <- 1 - pchisq(survival::survdiff(survival::Surv(time, status) ~ group, data = data, na.action = na.exclude)$chisq,df = length(levels(as.factor(data$group))) - 1)
    
    # mearge results
    tibble::tibble(
      group = rownames(.z$conf.int),
      n = .z$n,
      coef = .z$coefficients[,1],
      hr = .z$conf.int[1],
      hr_l = .z$conf.int[3],
      hr_h = .z$conf.int[4],
      coxp = .z$waldtest[3],
      kmp = kmp) %>%
      dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
  }else{
    tibble::tibble()
  }
  
  # multi-variable cox analysis
}

# 3.2.2.multi-variable cox analysis ---------
fn_survival_test.multiCox <- function(cancer_types,data,uni_sig_feature){
  print(cancer_types)
  covariates <- uni_sig_feature$features %>% unique()
  if(length(covariates)>=1){
    multi_formulas <- as.formula(paste('Surv(time, status)~', paste(covariates,collapse = " + ")))
    model <- coxph(multi_formulas, data = data, na.action = na.exclude)
    model.s <- summary(model)
    
    # Extract data
    tibble::tibble(
      Features = rownames(model.s$coefficients),
      n = model.s$n,
      coef = model.s$coefficients[,1],
      hr = model.s$conf.int[,1],
      hr_l = model.s$conf.int[,3],
      hr_h = model.s$conf.int[,4],
      coxp = model.s$coefficients[,5]) %>%
      dplyr::mutate(status = ifelse(hr > 1, "High_risk", "Low_risk"))
  }else{
    tibble::tibble()
  }
}

# 3.2.3.drawing cox plot ------
fn_cox_plot_1 <- function(data,filename,title,facet, dir,w=4,h=4){
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    # dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = Features, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6),
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2",1,2,3,4,5,6)) +
    # facet_grid(as.formula(facet),scales = "free", space = "free") +
    facet_wrap(as.formula(facet),scales = "free") +
    coord_flip() +
    # ggthemes::theme_gdocs() +
    my_theme +
    labs(y = "Hazard Ratio (High vs. low GSVA score)", x = "Features",title = title)+
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black"),
      title = element_text(size = 12)
    )  -> p;p
  ggsave(file.path(dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
}

fn_cox_plot_2 <- function(data,filename,title, dir,w=4,h=4){
  data %>%
    dplyr::arrange(hr) %>%
    .$Features -> features_rank
  data <- within(data,Features<- factor(Features,levels=features_rank))
  with(data,levels(Features))
  data %>% 
    # dplyr::mutate(functionWithImmune=functionWithImmune) %>%
    # dplyr::mutate(cancer_types = reorder(cancer_types,hr,sort)) %>%
    ggplot(aes(y = hr, x = Features, ymin=hr_l,ymax=hr_h)) +
    geom_pointrange(aes(color=cox_sig),size=0.5) +
    scale_color_manual(values=c("red","black")) +
    geom_hline(aes(yintercept = 1), linetype =2) +
    scale_size(name = "p-value") +
    scale_y_continuous(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6),
                       labels = c("1/64","1/32","1/16","1/8","1/4","1/2",1,2,3,4,5,6)) +
    # facet_grid(as.formula(facet),scales = "free", space = "free") +
    # facet_wrap(as.formula(facet),scales = "free") +
    coord_flip() +
    # ggthemes::theme_gdocs() +
    my_theme +
    theme(
      legend.position = "none",
      axis.line.y = element_line(color="black"),
      axis.text = element_text(color = "black",size=8),
      axis.title = element_text(color = "black",size=10),
      text = element_text(color = "black"),
      plot.title = element_text(size = 12)
    ) +
    labs(y = "Hazard Ratio (High vs. low value)", x = "Features",title = title) -> p;p
  ggsave(file.path(dir,paste(filename,"png",sep=".")),device = "png",width = w,height = h)
  ggsave(file.path(dir,paste(filename,"pdf",sep=".")),device = "pdf",width = w,height = h)
}
