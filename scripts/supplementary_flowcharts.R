library(DiagrammeR)
library(DiagrammeRsvg)
library(tidyverse)
library(rsvg)

# Compositional turnover flowchart
turnover <- grViz(
  diagram = "digraph flowchart {
      # define node aesthetics
      node [fontname = Arial, style = filled, fontsize = 19]                
      tab1 [label = '@@1', shape = rectangle, style = filled, fillcolor = beige]
      tab2 [label = '@@2', shape = rectangle, style = filled, fillcolor = beige]
      tab3 [label = '@@3', shape = rectangle, style = filled, fillcolor = beige]
      tab4 [label = '@@4', shape = rectangle, style = filled, fillcolor = beige]
      tab5 [label = '@@5', shape = rectangle, style = filled, fillcolor = beige]
      tab6 [label = '@@6', shape = rectangle, style = filled, fillcolor = beige]
      tab7 [label = '@@7', shape = rectangle, style = filled, fillcolor = beige]
      tab8 [label = '@@8', shape = parallelogram, style = filled, fillcolor = lightgrey]
# set up node layout
      tab1 -> tab2;
      tab2 -> tab3;
      tab3 -> tab4;
      tab4 -> tab5;
      tab5 -> tab6;
      tab6 -> tab7;
      tab8 -> tab1
     }
      [1]: 'Bray-Curtis \\n (Time[1] vs Time[2])'
      [2]: 'Calculate residual Bray-Curtis \\n (accounts for duration)'
      [3]: 'Add constant to residual Bray-Curtis\\n (bounds values between [0,1])'
      [4]: 'Logit transformation'
      [5]: 'Diversity modelling'
      [6]: 'Make counterfactual \\n predictions'
      [7]: 'Inverse logit transformation \\n of predictions'
      [8]: 'Harmonised site-level pollen matrix'
      
      "
)
turnover

# Whole study flowchart
study <- grViz(
  diagram = "digraph flowchart {
      # define node aesthetics
      node [fontname = Arial, style = filled, fontsize = 26]                
      tab1 [label = '@@1', shape = parallelogram, style = filled, fillcolor = lightgrey]
      tab2 [label = '@@2', shape = rectangle, style = filled, fillcolor = beige]
      tab3 [label = '@@3', shape = rectangle, style = filled, fillcolor = beige]
      tab4 [label = '@@4', shape = rectangle, style = filled, fillcolor = beige]
      tab5 [label = '@@5', shape = parallelogram, style = filled, fillcolor = lightgrey]
      tab6 [label = '@@6', shape = rectangle style = filled, fillcolor = beige]
      tab7 [label = '@@7', shape = rectangle, style = filled, fillcolor = beige]
      tab8 [label = '@@8', shape = parallelogram, style = filled, fillcolor = lightgrey]
      tab10 [label = '@@10', shape = rectangle, style = filled, fillcolor = beige]
      tab11 [label = '@@11', shape = rectangle, style = filled, fillcolor = beige]
      tab12 [label = '@@12', shape = rectangle, style = filled, fillcolor = beige]
     

# set up node layout
      tab1 -> tab2;
      tab2 -> tab3;
      tab8 -> tab7;
      tab10 -> tab6;
      tab2 -> tab4;
      tab3 -> tab10;
      tab4 -> tab10;
      tab6 -> tab7;
      tab7 -> tab11;
      tab5 -> tab2
      tab11 -> tab12
      }
      [1]: 'LegacyPollen\\n 1.0 pollen \\ndataset'
      [2]: 'Filter high quality pollen\\n samples and records'
      [3]: 'Resample pollen samples'
      [4]: 'Build age-depth models'
      [5]: 'LegacyAge\\n 1.0 chronological \\ndataset'
      [6]: 'Calculate raw diversity metrics\\n (richness, evenness, compositional turnover)'
      [7]: 'Diversity modelling'
      [8]: 'Palaeo data \\n (temperature change, \\nprecipitation change, \\narch-dates)'
      [9]: 'LegacyAge\\n 1.0 datings \\ndataset'
      [10]: 'Join 1000 resampled global pollen datasets\\n with draw from age-depth models'
      [11]: 'Make counterfactual \\n predictions'
      [12]: 'Compare empirical data with full and \\ncounterfactual model predictions'

      

      "
)
study




