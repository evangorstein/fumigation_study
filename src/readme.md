# Analyses
<!-- See Notion for more details. -->

All the following analyses are for bacteria.

1. PERMANOVA analysis: See `permanova.Rmd` for both bacteria and fungi.
    - Conclusion: There are significant differences by fumigation but not by treatment
2. Identify taxa that are different per group: See `bact_heatmap.Rmd` and `fung_heatmap.Rmd` for the heatmap visualization. Also, see `nmds.Rmd` for NMDS analysis for both bacteria and fungi.
3. Compare alpha-diversity (Shannon index and inverse Simpson index) by fumigation and treatment with statistical testing and linear regression: See `shannon_diversity_analysis.Rmd` and `invsimpson_diversity_analysis.Rmd`. 
   - Conclusion: For linear regression analysis, the model is appropriate for alpha-diversity and fumigation, meaning, a change in fumigation status affects alpha-diversity. The model is not appropriate for treatment.
   - For statistical testing, the results are there is significant difference in alpha-diversity for different fumigation status, but not for treatment.
4. Microbiome network analysis: See `network_analysis_bact.Rmd` and `network_analysis_fung.Rmd`.
   - Contains networks constructed for each fumigation status using glasso method, subsequent sub-network analysis for taxa unique to each fumigation status, and taxa in all fumigation status.
   - Identified taxa that changed the most between fumigation status and constructed CARlasso for those taxa.
  
<!-- The same analyses are applied to fungi.
1. PERMANOVA: See `permanova.Rmd`
    - Conclusion: There are significant differences by fumigation but not by treatment
2. Identify taxa that are different per group: See `fung_visual.Rmd` for the heatmaps and `nmds.Rmd` for NMDS analysis.
3. Compare alpha-diversity (Shannon index and inverse Simpson index) by fumigation and treatment with a linear model: See `regression_shannon.Rmd` and `regression_simpson.Rmd`.
   - Conclusion: Linear model is appropriate for alpha-diversity and fumigation, meaning, a change in fumigation status affects alpha-diversity. The linear model is not appropriate for treatment.
4. Microbiome network analysis: See `combined_graph_fungi.Rmd`
   - Contains networks constructed for each fumigation status using 3 methods.
   - Identified taxa that changed the most between fumigation status and constructed CARlasso for those taxa. -->
   
