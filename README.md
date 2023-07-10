# ResMAG

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


ResMAG is a state-of-the-art and user-friendly Snakemake workflow designed for the analysis of metagenomic data. It integrates multiple bioinformatics tools and algorithms to facilitate key steps in metagenome analysis, including bin refinement, metagenome-assembled genome (MAG) reconstruction, taxonomic classification of MAGs, and identification of antibiotic resistance genes.<br />

**Key Features:**<br />

**Binning Techniques**: Employ a collection of five state-of-the-art binning tools to partition metagenomic contigs into individual bins, allowing for comprehensive and accurate analysis.<br />

**MAG Reconstruction**: Utilize cutting-edge algorithms to reconstruct high-quality metagenome-assembled genomes (MAGs) from sequencing data.<br />

**Taxonomic Classification**: Apply advanced taxonomic classification methods to assign taxonomic labels to MAGs and identify the microbial community composition within the metagenomic samples.<br />

**Antibiotic Resistance Gene Identification**: Perform in-depth analysis to detect and characterize antibiotic resistance genes within the metagenomic data, providing valuable insights into antimicrobial resistance profiles.<br />

**Performance Refinement**: Continuously optimize the pipeline by incorporating the latest advancements in metagenomics research, ensuring the highest accuracy and efficiency in metagenomic data analysis.<br />


## Usage
### Preparations
To prepare the workflow
1. Clone it to your desired working folder via git or your preferred IDE
2. Edit the `config/config.yaml` file:
   - Specify a project name (`project-name`)
   - Download a kraken reference database for the host read filtering and specify the path to it (`kraken-db`)
3. Provide a sample information in the `config/pep/samples.csv` file with keeping the header and format as `.csv`:

```
sample_name,fq1,fq2
sample1,path/to/your/fastq/sample1_R1.fastq.gz,path/to/your/fastq/sample1_R2.fastq.gz
```
### Run the workflow
```snakemake --use-conda --cores all --rerun-incomplete```

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

## Output

## Contributing

**Bug report**

**Feature request**

## License

ResMAG is released under the [BSD-2 Clause](https://www.open-xchange.com/hubfs/2_Clause_BSD_License.pdf?hsLang=en). Please review the license file for more details.

## Contact Information

For any questions, or feedback, please contact the project maintainer at katharina.block@uk-essen.de or josefa.welling@uk-essen.de. We appreciate your input and support in using and improving ResMAG.

## Acknowledgements

We would like to express our gratitude towards Adrian Doerr, Alexander Thomas, Johannes Köster, Ann-Kathrin Brüggemann and the IKIM who have contributed to the development and testing of ResMAG. Their valuable insights and feedback have been helpful throughout the creation of the workflow.
IKIM, Johannes Köster, Alexander Thomas, ?

## References

**Tools**
[CoverM](https://github.com/wwood/CoverM)<br />
[DAS Tool](https://doi.org/10.1038/s41564-018-0171-1)<br />
[fastp](https://doi.org/10.1093/bioinformatics/bty560)<br />
[FastQC](https://github.com/s-andrews/FastQC)<br />
[Kraken 2](https://doi.org/10.1186/s13059-019-1891-0)<br />
[MEGAHIT](https://doi.org/10.1093/bioinformatics/btv033)<br />
[MetaBAT 2](https://doi.org/10.7717%2Fpeerj.7359)<br />
[MetaBinner](https://doi.org/10.1186/s13059-022-02832-6)<br />
[MetaCoAG](https://doi.org/10.1101/2021.09.10.459728)<br />
[minimap2](https://doi.org/10.1093/bioinformatics/bty191)<br />
[pandas](https://doi.org/10.5281/zenodo.3509134)<br />
[Rosella](https://github.com/rhysnewell/rosella)<br />
[samtools](https://doi.org/10.1093/gigascience/giab008)<br />
[VAMB](https://doi.org/10.1038/s41587-020-00777-4)<br />
[MultiQC](https://doi.org/10.1093%2Fbioinformatics%2Fbtw354)<br />



**Literature**

[Not here yet](https://www.lipsum.com/feed/html)

## Citation

A paper is on its way. If you use ResMAG in your work before the paper, then please consider citing this GitHub.
