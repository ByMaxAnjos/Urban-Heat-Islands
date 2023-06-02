# Zoom City Carbon Model::High spatio-temporal resolution of Urban-Heat-Islands

## Introduction

As part of **Zoom City Carbon Model (ZCCM)**, we present the **ZCCM::UHI**, a set of R functions which models Urban Heat Island at high-definition and provides hourly interpolated air temperature data sourced from local air temperature readings, including community-sourced data and Local Climate Zone. The **UHI Geographic Information platform** communicates the outcomes of the model in an interactive manner, enabling users, stakeholders, research community, and civil society to explore summary statistic, zoom maps of air temperature and dashboard that be accessed via this [link]

Please note that **ZCCM::UHI** is currently undergoing peer-review, and caution is advised when interpreting its outcomes. Our methodology is based on Anjos, M.; Meier, F. City Carbon Budget and hourly net CO2 fluxes at 0.01º resolution for informed climate action(in preparation).

### People

The development of the ZCCM::UHI was led by [Dr. Max Anjos](https://www.researchgate.net/profile/Max-Anjos/research) and joined by Dr.Fred Meier, and it is hosted at the [Chair of Climatology, Institute of Ecology, Technische Universität Berlin](https://www.klima.tu-berlin.de/index.php?show=home_start&lan=en).

### Funding

This project is was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES) – Finance Code 001, and by the Alexander Von Humboldt Foundation.

### Contact

Please feel free to contact us if you have any questions or suggestions by emailing [maxanjos\@campus.ul.pt](mailto:maxanjos@campus.ul.pt). If you are interested in contributing to the development of the model, we welcome you to join our team.

Happy coding!

## Input and setting requirements

To ensure the model runs correctly, it is necessary to load the following inputs:

1.  Air temperature data .csv (required) with a minimum of four columns labeled *date*, *Latitude*, *Longitude*, *airT*. 
2.  Other variables (optional) should have the same date column recommendation.

Note that the model converts the date-time into a R-formatted version, e.g., "2023-03-13 11:00:00" or "2023-03-13".

The following dataframe is presented as follows:

```{r setup, include=TRUE}
air_UCON %>% head(10)
```
<img width="691" alt="Screenshot 2023-05-03 at 08 57 51" src="https://user-images.githubusercontent.com/94705218/235909499-82427b94-5f35-4e58-b08b-0418d6fb4f44.png">

