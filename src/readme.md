# Analyses
See Notion for more details.

All the following analyses are for bacteria.
1. PERMANOVA: See `permanova.Rmd`
    - Conclusion: There are significant differences by fumigation but not by treatment
2. Identify taxa that are different per group: See `bac_visual.Rmd` for the heatmaps and `nmds.Rmd` for NMDS analysis.
3. Compare alpha-diversity by fumigation and treatment with a linear model: See `regression.Rmd`
   - Conclusion: Linear model is appropriate for alpha-diversity and fumigation, meaning, a change in fumigation status affects alpha-diversity. The linear model is not appropriate for treatment.
4. Microbiome network analysis: See `combined_graph_bact.Rmd`
   - Contains networks constructed for each fumigation status using 3 methods.
   - Identified taxa that changed the most between fumigation status and constructed CARlasso for those taxa.
  
The same analyses are applied to fungi.
1. PERMANOVA: See `permanova.Rmd`
    - Conclusion: There are significant differences by fumigation but not by treatment
2. Identify taxa that are different per group: See `fung_visual.Rmd` for the heatmaps and `nmds.Rmd` for NMDS analysis.
3. Compare alpha-diversity by fumigation and treatment with a linear model: See `regression.Rmd`
   - Conclusion: Linear model is appropriate for alpha-diversity and fumigation, meaning, a change in fumigation status affects alpha-diversity. The linear model is not appropriate for treatment.
4. Microbiome network analysis: See `combined_graph_fungi.Rmd`
   - Contains networks constructed for each fumigation status using 3 methods.
   - Identified taxa that changed the most between fumigation status and constructed CARlasso for those taxa.
   
