# Case-control-causal-review
 
[R](https://www.r-project.org/) Code to accompany the paper:

Mésidor, M., Xu, M., Diop, A., Fantodji, C., Parent, M. É., & Keil, A. (2025). Use of causal inference methods in case–control studies: a methodology review. American Journal of Epidemiology, kwaf182. [link](https://academic.oup.com/aje/advance-article/doi/10.1093/aje/kwaf182/8239042)


### Contents
- **code/ccc_interceptcorrection.R**: intercept-adjusted logistic regression
- **code/ccc_pscores.R**: propensity score based methods
- **code/ccc_weighted.R**: case-control weighted maximum likelihood
- **code/ccc_tmle.R**: case-control weighted targete maximum likelihood
- **slides/SER-CSEB Webinar_PartI/22.pptx**: Slides from a 3-Dec-2025 seminar on the paper, co-sponsored by SER and CSEB


### Requirements

R version: 4.0+

R packages:

- dplyr
- survival
- modi