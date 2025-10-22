# Impact-based Early Warning of Mass Movements

<p align="center">
    <a href="https://doi.org/10.5194/egusphere-2025-4940">
        <img alt="Preprint DOI" src="https://img.shields.io/badge/Preprint-10.5194/egusphere--20255--4940-sienna"></a>
    <a href="https://doi.org/10.6084/m9.figshare.30271795.v1">
        <img alt="Figshare DOI" src="https://img.shields.io/badge/Supplement-m9.figshare.30271795.v1-a60845"></a>
    <a href="https://style.tidyverse.org">
        <img alt="Code style: tidyverse" src="https://img.shields.io/badge/codestyle-tidyverse-darkblue"></a>
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
(2025):
**Impact-based early warning of mass movements – A dynamic spatial modelling approach for the Alpine region**.
*EGUsphere* [preprint]. [doi:10.5194/egusphere-2025-4940](https://doi.org/10.5194/egusphere-2025-4940).


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
