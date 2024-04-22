arg = commandArgs(T);
print(arg)
input_trait = arg[1];
enigma_trait = arg[2];
bfile_reference = arg[3]; # "/Volumes/faculty/njahansh/nerds/ravi/genetics/LAVA-main/support_data/eur/g1000_eur"
input_trait_name_conv = unlist(strsplit(arg[4]," ")); #[snp,beta,se,A1,A2,p]
out.fname = arg[5];
output_path=arg[6];

input_trait_name=  tail(strsplit(strsplit(tail(strsplit(input_trait, "/")[[1]], 1), "\\.")[[1]][1], "____")[[1]], 1)
file_enigma_name=  unlist(strsplit(strsplit(enigma_trait, "ENIGMA3_mixed_se_")[[1]][2], "_"))
Enigma_trait <- paste(file_enigma_name[3], file_enigma_name[4], sep = "_")

# input_trait ="//smb-ifs.ini.usc.edu/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/METAL_Out/METAL_output_files/ukbb_white_abcd_white____Total_MeanThickness1.tbl.gz" 
# enigma_trait ="//smb-ifs.ini.usc.edu/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/LAVA_Module/Enigma//ENIGMA3_mixed_se_wSA_Mean_inferiortemporal_surfavg_20190429.txt.gz" 
# bfile_reference ="//smb-ifs.ini.usc.edu/faculty/njahansh/nerds/ravi/genetics/LAVA-main/support_data/eur/g1000_eur" 
# input_trait_name_conv =unlist(strsplit("MarkerName Effect StdErr Allele1 Allele2 PvalueARE"," ")) 
# out.fname ="Total_MeanThickness1__inferiortemporal" 
# output_path= "//smb-ifs.ini.usc.edu/faculty/njahansh/nerds/ankush/GiNi_post_GWAS_processing/Exta_temp_files/CAUSE_Out/Total_MeanThickness1/wSA"

output_file_name <- paste0(output_path,"/", out.fname);

print(input_trait);
print(enigma_trait);
print(bfile_reference);
print(input_trait_name_conv[[6]]);
print(out.fname);
print(output_path);
print(output_file_name);


library(cause)
library(tidyverse)
library(genetics.binaRies)
library(ieugwasr)



################################################################################Get Gene name
library(httr)
library(jsonlite)
library(xml2)

################################################################################


################################################################################Cause Plot Updated
# install.packages("scales")
# install.packages("gridExtra")


# if (!requireNamespace("gridExtra", quietly = TRUE)) {
#   # If not installed, install the package
#   install.packages("gridExtra")
# }
# 
# if (!requireNamespace("scales", quietly = TRUE)) {
#   # If not installed, install the package
#   install.packages("scales")
# }

library(gridExtra)

library(scales)
library(ggrepel)
library(ggpubr)


################################################################################



set.seed(100);
r2_thresh = 0.01;
pval_thresh = 1e-3;

input_trait_file <- read_tsv(input_trait);
Enigma_trait_file <- read_tsv(enigma_trait);



# Format Data for CAUSE
input_trait_Enigma <- gwas_merge(input_trait_file, Enigma_trait_file, snp_name_cols = c(input_trait_name_conv[[1]], "SNP"),
                                 beta_hat_cols = c(input_trait_name_conv[[2]], "BETA"),
                                 se_cols = c(input_trait_name_conv[[3]], "se"),
                                 A1_cols = c(input_trait_name_conv[[4]], "A1"),
                                 A2_cols = c(input_trait_name_conv[[5]], "A2"),
                                 pval_cols = c(input_trait_name_conv[[6]], "P"));

Enigma_input_trait <- gwas_merge(Enigma_trait_file, input_trait_file, snp_name_cols = c("SNP", input_trait_name_conv[[1]]),
                                 beta_hat_cols = c("BETA", input_trait_name_conv[[2]]),
                                 se_cols = c("se", input_trait_name_conv[[3]]),
                                 A1_cols = c("A1", input_trait_name_conv[[4]]),
                                 A2_cols = c("A2", input_trait_name_conv[[5]]),
                                 pval_cols = c("P", input_trait_name_conv[[6]]));


# (COPC, INT, snp_name_cols = c("SNP", "SNP"), 
#                        beta_hat_cols = c("b", "b"), 
#                        se_cols = c("se", "se"), 
#                        A1_cols = c("A1", "A1"), 
#                        A2_cols = c("A2", "A2"), 
#                        pval_cols = c("p", "p")
# )

# Calculate nuisance parameters
varlist_Enigma_input_trait <- with(Enigma_input_trait, sample(snp, size=1000000, replace=FALSE));
params_Enigma_input_trait <- est_cause_params(Enigma_input_trait, varlist_Enigma_input_trait);

varlist_input_trait_Enigma <- with(input_trait_Enigma, sample(snp, size=1000000, replace=FALSE));
params_input_trait_Enigma <- est_cause_params(input_trait_Enigma, varlist_input_trait_Enigma);

# LD Pruning
input_trait_Enigma_clump <- input_trait_Enigma %>%
  rename(rsid = snp,
         pval = p1,
  ) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile_reference);

top_vars_input_trait_Enigma <- input_trait_Enigma_clump$rsid;


Enigma_input_trait_clump <- Enigma_input_trait %>%
  rename(rsid = snp,
         pval = p1,
  ) %>%
  ieugwasr::ld_clump(dat = .,
                     clump_r2 = r2_thresh,
                     clump_p = pval_thresh,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile_reference);


top_vars_Enigma_input_trait <- Enigma_input_trait_clump$rsid;


################################################################################Get Gene name

server <- "https://rest.ensembl.org"
ext <- "/vep/human/id"

process_gene_name <- function(res_snp_df) {
  # Extract IDs from the input dataframe
  ids <- res_snp_df$snp
  
  # Convert the IDs to a JSON array
  json_data <- toJSON(list(ids = ids))
  
  # Make the POST request with the JSON data
  r <- POST(paste(server, ext, sep = ""), 
            content_type("application/json"), 
            accept("application/json"), 
            body = json_data)
  
  # Check for HTTP status code
  stop_for_status(r)
  
  # Parse JSON response
  parent_df <- fromJSON(toJSON(content(r)))
  
  # Denormalize the dataframe
  denormalized_df <- parent_df %>%
    mutate(row_number = row_number()) %>%
    unnest(cols = transcript_consequences, names_sep = "_") %>%
    select(-row_number)
  
  # Select subset of columns
  denormalized_df_sub <- denormalized_df %>%
    select(transcript_consequences_gene_symbol, transcript_consequences_gene_symbol_source, input, id)
  
  # Replace <Null> values with empty strings and extract first element from lists
  denormalized_df_sub_chr <- denormalized_df_sub %>%
    mutate(transcript_consequences_gene_symbol = map_chr(transcript_consequences_gene_symbol, ~ ifelse(is.null(.x), "", .x[[1]])),
           transcript_consequences_gene_symbol_source = map_chr(transcript_consequences_gene_symbol_source, ~ ifelse(is.null(.x), "", .x[[1]])),
           input = sapply(input, `[`, 1),
           id = sapply(id, `[`, 1))
  
  return(denormalized_df_sub_chr)
}


############################# ONLY ELPD PLOT + Gene Name
#'@title Plot posteriors for one CAUSE fit
#'@param fit object of class cause_post `
#'@param type Either "posteriors" or "data". See details.
#'@param data If type=="data" then an object of type cause_data_fit must be supplied
#'@param pval_thresh p-value threshold used when type=="data".
#'@details If type == "posteriors", the function plots the marginal posterior
#'distribution of each parameter. If type=="data", the function plots data colored
#'by the probability that Z = 1.
#'@export
plot.cause_post <- function(fit, type=c("posteriors", "data"), data=NULL, pval_thresh=5e-8){
  type = match.arg(type)
  if(type == "posteriors"){
    marge_post <- map_df(seq_along(fit$marge_post), function(i){
      p <- fit$marge_post[[i]]
      p$param <- fit$params[i]
      return(p)})
    marge_post <- select(marge_post, mid, width, param, post, prior) %>%
      gather("dist", "pdf", -mid, -param, -width)
    marge_post$param <- factor(marge_post$param, levels =c("gamma", "eta", "q"))
    plt <- ggplot(marge_post) + geom_line(aes(x=mid, y=pdf/width, linetype=dist)) +
      xlab("Parameter Value") + ylab("Density") +
      theme_bw() + theme(legend.title = element_blank()) +
      facet_wrap(~param, scale="free")
    return(plt)
  }
  if(type=="data"){
    stopifnot(inherits(data, "cause_data_fit"))
    if(all(c("gamma", "eta") %in% fit$params)){
      #Causal model
      medians <- data.frame(param = c("gamma", "eta"), med = summary(fit)$quants[1, c(1, 2)])
      var <- "prob_Z1_causal"
      title <- "Causal Model"
    }else if("eta" %in% fit$params){
      #Sharing model
      medians <- data.frame(param = c( "eta"), med = summary(fit)$quants[1, 2])
      var <- "prob_Z1_sharing"
      title <- "Sharing Model"
    }
    
    plt <- data %>% mutate(pval1 = 2*pnorm(-abs(beta_hat_1/seb1))) %>%
      filter(pval1 < pval_thresh) %>%
      ggplot(.) +
      geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
      geom_errorbar(aes(ymin = beta_hat_2 -1.96*seb2, ymax = beta_hat_2 + 1.96*seb2, x = beta_hat_1 ), color="grey") +
      geom_errorbarh(aes(y = beta_hat_2, xmin = beta_hat_1 - 1.96*seb1, xmax = beta_hat_1 + 1.96*seb1), color="grey") +
      geom_abline(aes(slope = med, linetype=param, intercept=0), data=medians) +
      geom_point(aes(x=beta_hat_1, y=beta_hat_2, col=get(var), size = -log10(pval1))) +
      scale_color_continuous(limits=c(0, 1), name = "P(Z = 1)", low="grey", high=muted("red")) +
      ggtitle(title) +
      theme_bw()
    return(plt)
  }
}

#'@title Plot CAUSE
#'@description Plot posteriors for sharing and causal models, display summary tables
#'@param res object of class cause
#'@param intern If TRUE, function returns a list of grobs. Otherwise it plots.
#'@param type Either "posteriors" or "data". See details.
#'@param pval_thresh p-value threshold used if type=="data".
#'@details If type == "posteriors", the function will plot the posterior distributions of the parameters
#'and display summary tables giving medians and credible intervals. If type == "data" the function
#'will plot the data thresholded on the trait 1 p-value if using pval_thresh.
#'@export
plot.cause <- function(res,gene_txt, intern=FALSE, type=c("posteriors", "data"), pval_thresh = 5e-8){
  type <- match.arg(type)
  if(type == "posteriors"){
    plts <- lapply(c("sharing", "causal"), function(i){
      plt <- plot(res[[i]])
      return(plt)})
    elpd <- res$elpd %>% mutate(p = pnorm(z),
                                delta_elpd = signif(delta_elpd, digits=2),
                                se_delta_elpd = signif(se_delta_elpd, digits=2),
                                z = signif(z, digits=2),
                                p = signif(p, digits=2))
    elpd <- tableGrob(elpd)
    tab <- tableGrob(summary(res)$tab)
    plts[[3]] <- tab
    plts[[4]] <- elpd
    if(intern) return(plts)
    h <- arrangeGrob(grobs = plts,
                     layout_matrix = rbind(c(4, 4, 4, 3, 3, 3),
                                           c(NA, NA, 1, 1, 1, 1),
                                           c(2, 2, 2, 2,2,2)))
    plot(h)
  }
  if(type=="data"){
    
    
    gene_txt2 <- distinct(gene_txt[, c("input","transcript_consequences_gene_symbol")])
    
    gene_txt2 <- gene_txt2 %>%
      group_by(`input`) %>%
      slice(1) %>%
      ungroup()
    
    func_fixer <- function(row) {
      if (is.na(row$transcript_consequences_gene_symbol) || row$transcript_consequences_gene_symbol == "-") {
        return(row$`input`)
      }
      return(row$transcript_consequences_gene_symbol)
    }
    
    gene_txt2$F_SYMBOL <- sapply(split(gene_txt2, seq(nrow(gene_txt2))), func_fixer)
    
    
    gene_txt2 <- gene_txt2 %>%
      rename(snp =`input` )
    
    
    res$data$SYMBOL <- left_join(res$data, gene_txt2, by ="snp")$F_SYMBOL
    # print(res)
    medians2 <- data.frame(param = c("gamma"), med = summary(res[["causal"]])$quants[1, c(1)])
    
    
    # plts <- lapply(c("sharing", "causal"), function(i){
    #   plt <- plot(res[[i]], type="data", data=res$data, pval_thresh=pval_thresh)
    #   return(plt)})
    max_delta_elpd <- max(abs(res$data$delta_elpd))
    plts <-  res$data %>%
      mutate(pval1 = 2*pnorm(-abs(beta_hat_1/seb1))) %>%
      filter(pval1 < pval_thresh) %>%
      
      ggplot(data=.,aes(x=beta_hat_1,y=beta_hat_2, label = SYMBOL)) +
      # geom_point()+
      geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
      geom_errorbar(aes(ymin = beta_hat_2 -1.96*seb2, ymax = beta_hat_2 + 1.96*seb2, x = beta_hat_1 ), color="grey") +
      geom_errorbarh(aes(y = beta_hat_2, xmin = beta_hat_1 - 1.96*seb1, xmax = beta_hat_1 + 1.96*seb1), color="grey") +
      geom_abline(aes(slope = med, intercept=0), data=medians2, linetype="dashed")+
      geom_point(aes(x=beta_hat_1, y=beta_hat_2, col=delta_elpd, size = -log10(pval1)))+
      geom_label_repel(color = "#E18727") +
      scale_color_gradient2(name = "Contribution\nto test statisitc", mid = "grey", limits=c(-1, 1)*max_delta_elpd) +
      ggtitle("ELPD Contribution") +
      theme_bw()
    if(intern) return(plts)
    # h <- arrangeGrob(grobs = plts, nrow=1)
    plot(plts)
  }
}



###########################################################Updated_Plot_end
################################################################################Get Gene name_end


# Fit CAUSE
res_input_trait_Enigma <- cause(X=input_trait_Enigma, variants = top_vars_input_trait_Enigma, param_ests = params_input_trait_Enigma)
summary(res_input_trait_Enigma, ci_size=0.95)

res_input_trait_Enigma_snps <- res_input_trait_Enigma$data %>%  filter(delta_elpd < 0.0) %>% filter(snp %in% top_vars_input_trait_Enigma[which(top_vars_input_trait_Enigma %in% (input_trait_Enigma %>% filter(p1 < 5e-8) %>% select(snp))$snp)])






res_Enigma_input_trait <- cause(X=Enigma_input_trait, variants = top_vars_Enigma_input_trait, param_ests = params_Enigma_input_trait)
summary(res_Enigma_input_trait, ci_size=0.95)

res_Enigma_input_trait_snps <- res_Enigma_input_trait$data %>% filter(snp %in% top_vars_Enigma_input_trait[which(top_vars_Enigma_input_trait %in% (Enigma_input_trait %>% filter(p1 < 5e-8) %>% select(snp))$snp)])


res_input_trait_Enigma_gene <-process_gene_name(res_input_trait_Enigma_snps)
res_Enigma_input_trait_gene <-process_gene_name(res_Enigma_input_trait_snps)


ggsave(

  paste(output_file_name,"_input_trait_Enigma",".png",sep=""),
  plot = plot(res_input_trait_Enigma,res_input_trait_Enigma_gene, type="data")+ labs(y = paste0("Association with ", input_trait_name," (b)"), x = paste0("Association with ", Enigma_trait," (b)"), title = ""),
  dpi = 600,
  bg = "transparent"
)

ggsave(

  paste(output_file_name,"_Enigma_input_trait",".png",sep=""),
  plot = plot(res_Enigma_input_trait,res_Enigma_input_trait_gene, type="data")+ labs(y = paste0("Association with ", Enigma_trait," (b)"), x = paste0("Association with ", input_trait_name," (b)"), title = ""),
  dpi = 600,
  bg = "transparent"
)


