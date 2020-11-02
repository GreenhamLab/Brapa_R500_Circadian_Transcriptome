# _Brassica rapa_ R500 Circadian Transcriptome Analysis 

## A detailed analysis of two separate transcriptomic data sets in _Brassica rapa_ using the novel R package 'DiPALM'.

This repository contains an R markdown script (Brapa\_CircadianTranscriptome\_Markdown.Rmd) detailing the bioinformatics, statistical analysis and figure plotting that is published in the following manuscript (available [here](https://elifesciences.org/articles/58993)):

```
Greenham K*, Sartor RC*, Zorich S, Lou P, Mockler TC, McClung CR (2020).
Expansion of the circadian transcriptome in Brassica rapa and genome wide diversification of paralog expression patterns.
eLife. 9:e58993 doi: 10.7554/eLife.58993
*These authors contributed equally. 
```
All other files in this repository are supplemental data or code files that are utilized in the R markdown script.

This analysis demonstrates the use of a novel R package to detect Differential Patterning Analysis via Linear Modeling (DiPALM). This package, along with a detailed vignette are available through the Comprehensive R Archive Network (CRAN) (https://cran.r-project.org/). 

Details of the two high-resolution time-course RNA-seq experiments used in this analysis are available in the R markdown file. The first data set is published along with this paper and represents transcriptomic abundance (RNA-seq) of _B. rapa_ over a 2-day time course after entrainment in two separate conditions (cycling light or cycling tempurature). The second data set is from a previously published study and represents transcriptomic abundance (RNA-seq) of _B. rapa_ over a 2-day time course in well watered (control) and mild drought (treatment) conditions. This data set is published in the following manuscript (available [here](https://elifesciences.org/articles/29655)):

```
Greenham K*, Guadagno CR*, Gehan MA, Mockler TC, Weinig C, Ewers, BE, McClung CR (2017). Temporal network analysis identifies early physiological and transcriptomic indicators of mild drought in Brassica rapa. 
eLife. 18(6). doi: 10.7554/eLife.29655 
*These authors contributed equally. 
```
### Getting Started

The file called 'Brapa\_CircadianTranscriptome\_Markdown.Rmd' is an R markdown file that contains detailed explanations and code to run this analysis. We recommend using the [Rstudio](https://rstudio.com/) integrated development environment to view and run this file.
More details can be found in the .Rmd file or in the more aesthetically pleasing .html file. All other files in this repository should be stored in a directory and pointed to by the 'rawDataPath' variable (line 32) in the R markdown file. This repository should contain everything needed to run through the entire analysis. 
