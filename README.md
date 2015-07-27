### OrthoMCL validation pipeline ###
========

This repository contains scripts for working with OrthoMCL clustering results. 


You can use the pipeline to identify the core, shared and lineage specific proteins in the pan-proteome. Core proteins are proteins present in all members of a clade while shared proteins are proteins that are present in one or more members. Lineage specific proteins do not have any orthologs in other genomes. We recommend runing OrthoMCL with varying stringency settings (E value 1e-5, percent identity 30%, 50% and 60%) and comparing the results of all runs to identify the consistent and reproducible clusters.

Please see Jones et al. 2015 for an example of how this pipeline can be used.
http://apsjournals.apsnet.org/doi/abs/10.1094/PDIS-08-14-0833-RE

Please create an issue or contact suryasaha at cornell dot edu if you have any questions. Thanks!

Here is an example workflow, warts and all
https://gist.github.com/suryasaha/949741ec6b541bf6f1ce
