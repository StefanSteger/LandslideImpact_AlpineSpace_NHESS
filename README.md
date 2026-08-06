# Impact-based Early Warning of Mass Movements

<p align="center">
    <a href="https://doi.org/10.5194/nhess-26-3637-2026">
        <img alt="Paper DOI" src="https://img.shields.io/badge/Paper-10.5194%2Fnhess--26--3637--2026-sienna?style=flat-square"></a>
    <a href="https://doi.org/10.6084/m9.figshare.30271795">
        <img alt="Figshare DOI" src="https://img.shields.io/badge/Supplement-10.6084%2Fm9.figshare.30271795-556472?style=flat-square&logo=figshare&logoColor=white"></a>
    <a href="https://doi.org/10.5281/zenodo.21645145">
        <img alt="ZENODO DOI" src="https://img.shields.io/badge/PID-doi:10.5281/zenodo.21645145-1682D4?style=flat-square&logo=zenodo&logoColor=white"></a>
    <a href="https://style.tidyverse.org">
        <img alt="Code style: tidyverse" src="https://img.shields.io/badge/codestyle-tidyverse-1a162d?style=flat-square&logo=r&logoColor=white"></a>
    </a>
</p>

This repository supplements the manuscript by
Stefan Steger<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-0886-5191)</sup>,
Raphael Spiekermann<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-4772-9750)</sup>,
Mateo Moreno<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-9530-3076)</sup>,
Sebastian Lehner<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-7562-8172)</sup>,
Katharina Enigl<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-3263-0918)</sup>,
Alice Crespi<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-4186-8474)</sup>
and
Matthias Schlögl<sup>[![](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0002-4357-523X)</sup>
(2026):
**Dynamic spatial modelling of mass movement impacts for large areas: A data-driven framework for impact-based early warning**.
*Natural Hazards and Earth System Sciences*, 26(8), 3637–3665. [doi:10.5194/nhess-26-3637-2026](https://doi.org/10.5194/nhess-26-3637-2026).


## Overview

The R code in this repository uses the provided data to:

- Fit and analyze models separately for **slide-types**, **flow-types**, and **fall-types**.  
- Visualize data and results.  
- Include example snippets for specific analyses, such as:
  - **Basin-based cross validation**  
  - **Variable importance assessment**  


## Script Structure

1. **Step 1:** Load data and create basic visualizations  
2. **Step 2:** Fit the models  
3. **Step 3:** Calculate fitting performance and plot ROC curves *(Fig. 4a in publication)*  
4. **Step 4:** Cross validation: basin-based CV and visualization of final performance *(Fig. 4a–d)*  
5. **Step 5:** Variable importance plots, including calculation examples *(Fig. 5)*  
6. **Step 6:** Visualize partial effects from fitted models *(Fig. 6, 7, 8)*
