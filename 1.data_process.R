####################################
# data process for TCGA data
# include TIL, expression data, and clinical data


# path --------------------------------------------------------------------

basic_path <- "/home/huff/project"
TIL_path <- file.path(basic_path,"/data/TCGA/immune_infiltration/miao_TCAP_prediction")
TCGA_path <- file.path("/home/huff/data/TCGA/TCGA_data") # pancan33_expr.rds.gz
data_path <- file.path(basic_path,"TCGA_nomogram/data")


# load data ---------------------------------------------------------------

# gene list exlude TIL prediction markers
inhibit_cells <- c("Exhausted","iTreg","Neutrophil","Monocyte","nTreg","Tr1")
non_inhibit_TIL_markers <- readr::read_tsv(file.path(data_path,"markers_used_to_predict_TIL_miao_TCAP.txt")) %>%
  tidyr::gather(key="Cell_type",value="markers")%>%
  dplyr::filter(!is.na(markers)) %>%
  dplyr::filter(! Cell_type %in% inhibit_cells)

checkpoints_as_TILMarker <- readr::read_tsv(file.path(data_path,"ICPs_all_info_class.tsv")) %>%
  dplyr::filter(! symbol %in% non_inhibit_TIL_markers$markers)

checkpoints <- readr::read_tsv(file.path(data_path,"ICPs_all_info_class.tsv"))

# expression data
genelist_exp <- readr::read_rds(file.path(TCGA_path,"pancan33_expr.rds.gz")) %>%
  dplyr::mutate(exp_filter = purrr::map2(expr,cancer_types,.f=function(.x,.y){
    print(.y)
    .x %>%
      dplyr::filter(symbol %in% checkpoints$symbol)
  })) %>%
  dplyr::select(-expr) 

# TIL data
# TIL <- readr::read_rds(file.path(TIL_path,"pancan33_immune_infiltration_by_TCAP.rds.gz")) 
TIL.new <- readr::read_rds(file.path(data_path,"pancan33_cancer_TIL-191108-liull.rds.gz")) %>%
  dplyr::mutate(Infiltration = purrr::map(cancer_TIL, .f=function(.x){
    .x %>%
      dplyr::rename("barcode"="Sample_ID")
  })) %>%
  dplyr::select(-cancer_TIL)

# clinical data
clinical_2018cell <- readr::read_rds(file.path("/home/huff/project/data/TCGA-survival-time/cell.2018.survival","TCGA_pancan_cancer_cell_survival_time.rds.gz"))
clinical_TCGA <- readr::read_rds(file.path("/home/huff/project/TCGA_survival/data","Pancan.Merge.clinical.rds.gz"))

stage_class <- tibble::tibble(Stage=c("stage i","stage ia","stage ib","stage ic","stage ii","stage iia",
                       "stage iib","stage iic","stage iii","stage iiia","stage iiib",
                       "stage iiic","stage iv","stage iva","stage ivb","stage ivc","stage x" ),
               Stage1 = c("stage i","stage i","stage i","stage i","stage ii","stage ii",
                          "stage ii","stage ii","stage iii","stage iii","stage iii",
                          "stage iii","stage iv","stage iv","stage iv","stage iv","stage x"))
clinical_2018cell %>%
  dplyr::mutate(PFS = purrr::map(data,.f=function(.x){
    .x %>%
      dplyr::select(bcr_patient_barcode, PFS, PFS.time) %>%
      dplyr::rename("barcode"="bcr_patient_barcode")
  })) %>%
  dplyr::select(-data) %>%
  dplyr::rename("cancer_types"="type") %>%
  dplyr::inner_join(clinical_TCGA,by="cancer_types") %>%
  dplyr::mutate(OS_stage = purrr::map(clinical_data,.f=function(.x){
    .x %>%
      dplyr::inner_join(stage_class,by="Stage") %>%
      dplyr::select(-Stage) %>%
      dplyr::rename("Stage"="Stage1")
  })) %>%
  dplyr::select(-clinical_data)-> clinical

# combine data
clinical %>%
  dplyr::inner_join(TIL.new, by = "cancer_types") %>%
  dplyr::inner_join(genelist_exp, by="cancer_types") -> combined_clinical_TIL_exp_data

fn_data_process_main <- function(cancer_types,PFS,OS_stage,Infiltration,exp_filter){
  print(cancer_types)
  Infiltration %>%
    tibble::rowid_to_column(var = "rowid") %>%
    tidyr::gather(-barcode,-rowid, key="features", value="value") %>%
    dplyr::filter(substr(barcode,14,15)=="01") %>%
    dplyr::mutate(barcode_1 = substr(barcode,1,12)) %>%
    dplyr::group_by(barcode_1,features) %>%
    dplyr::filter(rowid == min(rowid)) %>%
    .$barcode %>%
    unique() -> barcode_uniq
  
  exp_filter %>%
    tidyr::gather(-symbol, -entrez_id, key="barcode", value="value") %>%
    dplyr::filter(barcode %in% barcode_uniq) %>%
    dplyr::mutate(barcode = substr(barcode,1,12)) %>%
    dplyr::select(-entrez_id) %>%
    dplyr::select(barcode,symbol,value) %>%
    tidyr::spread(key="barcode", value="value")  %>%
    dplyr::mutate(x="x") %>%
    tidyr::nest(-x,.key="exp_filter") -> exp
  Infiltration %>%
    tidyr::gather(-barcode,key="features", value="value") %>%
    dplyr::filter(barcode %in% barcode_uniq) %>%
    dplyr::mutate(barcode = substr(barcode,1,12)) %>%
    dplyr::select(barcode,features,value) %>%
    tidyr::spread(key="features", value="value") %>%
    dplyr::mutate(x="x") %>%
    tidyr::nest(-x,.key="Infiltration")-> TIL
  PFS %>%
    dplyr::mutate(x="x") %>%
    tidyr::nest(-x,.key="PFS")-> PFS
  OS_stage %>%
    dplyr::mutate(x="x") %>%
    tidyr::nest(-x,.key="OS_stage")-> OS_stage
  TIL %>%
    dplyr::inner_join(exp,by="x") %>%
    dplyr::inner_join(OS_stage,by="x") %>%
    dplyr::inner_join(PFS,by="x") %>%
    dplyr::select(-x)
}

combined_clinical_TIL_exp_data %>%
  dplyr::mutate(combine_data = purrr::pmap(list(cancer_types,PFS,OS_stage,Infiltration,exp_filter),fn_data_process_main)) %>%
  dplyr::select(cancer_types,combine_data) %>%
  tidyr::unnest() -> combined_clinical_TIL_exp_data.new
  
combined_clinical_TIL_exp_data.new %>%
  readr::write_rds(file.path(data_path,"TCGA_combined_clinical_TIL_exp_data.rds.gz"))
