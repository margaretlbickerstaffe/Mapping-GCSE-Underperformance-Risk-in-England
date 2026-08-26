DataWrangling Files: These notebooks were used to produce the final msoaAgg{period_label}.csv files used in the R scripts for the cross-sectional analysis. Each folder contains the dataset outputs as well as the Python scripts used to create said outputs. 
SecondarySchools.ipynb
Filters the most up-to-date list of schools within the UK to non-selective, non-PRU secondary schools, sourced from the Department for Education's Key Stage 4 Performance Data by school.
Outputs finalSecondarySchools.csv
YouthProvisionList.ipynb (updating work done by Irtania Ismira Dewi in her open-access code)
Reads in data from the Charity Commission, Companies House, and ONSPD datasets, which are linked in the dissertation data source section but too large to upload to GitHub.
Filters based on broad and strict inclusion criteria.
Produces four charity lists for each year in the case study.
Outputs:
df_final_chc_broad{period_label}.csv
df_final_chc_strict{period_label}.csv
df_coh_final_broad{period_label}.csv
df_coh_final_strict{period_label}.csv
WalkableCharities.ipynb
Reads in finalSecondarySchools.csv, each of the df_final_coh/chc_strict/broad{period_label}.csv files, and the 2016/2017 Key Stage 4 performance data.
Filters the 2016/2017 schools against finalSecondarySchools.csv to establish a list of non-selective, non-PRU secondary schools that operated continuously from 2016/2017 to the 2024/2025 school year, assuming admissions type and PRU status remained stable across the study period.
Adds coordinates to both schools and youth provisions in each list using the ONSPD file.
For each school, a list of youth provisions within a one-mile radius is developed, then cross-referenced with walking distance via the OpenStreetMap API to produce a count of walkable charities.
For each school and study year, the count of walkable charities across the four lists is Borda ranked to produce a final ranking of schools most saturated by youth provision across England.
Outputs: finalDf{period_label}.csv
OutsideVariables.ipynb
Reads in Key Stage 4 performance data for each year, as well as IDACI, SAMHI, and race data linked in the dissertation data source section.
Groups each covariate by MSOA.
Merges by secondary school as well.
Outputs: msoaAgg{period_label}.csv and schoolData{period_label}.csv for secondary-school-level analysis.
RScripts Files: This folder contains five R scripts used to complete the repeated cross-sectional Bayesian risk analysis. These scripts use data from the msoaAgg{period_label}.csv files, as well as various shapefiles linked in the dissertation data source section.
CaseStudyAnalysis Files: This folder contains the outputs from the repeated cross-sectional Bayesian risk analysis. These CSV files were used to perform case study analysis in Python. The Python notebook used for this analysis, RiskProfilingAnalysisByMSOA.ipynb, is also included in this folder.
