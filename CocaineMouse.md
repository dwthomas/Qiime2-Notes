# 16S Metabarcoding with Qiime 2
In this example we'll go over how to use QIIME 2 to analyze metabarcoding data.
These data are from set of mouse fecal samples provided by [Jason Bubier from The Jackson Laboratory](https://www.jax.org/research-and-faculty/faculty/research-scientists/jason-bubier).
The samples were run targeting the V1-V3 region of the 16S gene using the 27F 5'-AGAGTTTGATCCTGGCTCAG-3' and 534R 5'-TTACCGCGGCTGCTGG-3' primers.

Our metadata is available to view as a [google sheet](https://docs.google.com/spreadsheets/d/1Tj_wj9h19JU6ShbxPu-2kQ3g0r2hf-YRQbm-R4-xnFk/edit?usp=sharing), we can use the [keemei](https://keemei.qiime2.org/) sheets add-on to validate the metadata.
This data was collected to measure change in the fecal bacteria in laboratory mice who self administered cocaine.  The mice are from two strains, CC004 mice are high-responding, while CC041 are low responding.  More details on the strains responses can be found [here](https://link-springer-com.unh.idm.oclc.org/article/10.1007/s00213-019-05429-3). There were male and female mice of both strains, and the mice were sampled before and after treatment.  
The treatments were the self administration of cocaine, or a saline sham.  

We will primarily use the [Qiime 2](https://qiime2.org/) bioinformatics platform.
Qiime 2 is free and open source and available from Linux and OSX.
We will use the Qiime2 command line interface, there is also the ["Artifact" python API](https://docs.qiime2.org/2019.4/interfaces/artifact-api/) which can be more powerful.
## Getting the Data
We start by activating the Qiime 2 environment.  The server is a shared resource and we may want to be able to use different version of programs, like blast or R or Python than Qiime 2 requires.  To enable this Qiime 2 is given its own working environment with the exact version of all the programs it requires.  Qiime 2 currently puts out a new version about every 3 months.  You should upgrade versions as they come available, even if you began with an earlier version.
~~~bash
#      version: qiime2-year.month
source activate qiime2-2021.4
~~~
Now lets grab a copy of the data!  Notice that the copy command will warn us that it is skipping the reads directory, that is OK!
~~~bash
cp -r /home/share/examples/cocaine_mouse .
cd cocaine_mouse
ls
~~~

Let's begin by confirming the primers are as expected, we can do that by looking at the raw reads.  We want to confirm that we are seeing the primers at the beginning of the majority of the reads.  

~~~bash
zcat reads/JBCDJ00OLL1STT0B00000191771C7M7FGT1904906_GATCAAGG_ACTCCATC_S3_R1_001.fastq.gz | less -S
zcat reads/JBCDJ00OLL1STT0B00000191771C7M7FGT1904906_GATCAAGG_ACTCCATC_S3_R2_001.fastq.gz | less -S
~~~




You should see the forward primers starting the forward reads, and the reverse primer starting the reverse reads.

Now we're ready to import the data into qiime.  First let's have a look at the anatomy of most qiime2 commands.  There is a consistent scheme to the command and the arguments (a big improvement from qiime1!).  Inputs always begin with --i-, outputs with --o-, metadata with --m- and other parameters with --p-.  We can use this and tab completion to easily guess and type out most commands, for instance if I forget the inputs a command takes, I can start typing --i-, then use tab-completion to display the options!  
~~~bash
# Anatomy of a qiime command
qiime plugin action\
   --i-inputs  foo\       ## input arguments start with --i
   --p-parameters bar\    ## paramaters start with --p
   --m-metadata mdat\     ## metadata options start with --m
   --o-outputs out        ## and output starts with --o
~~~
Qiime works on two types of files, Qiime Zipped Archives (.qza) and Qiime Zipped Visualizations (.qzv).  Both are simply renamed .zip archives that hold the appropriate qiime data in a structured format.  This includes a "provenance" for that object which tracks the history of commands that led to it.  The qza files contain data, while the qzv files contain visualizations displaying some data.  We'll look at a quality summary of our reads to help decide how to trim and truncate them.
~~~bash
qiime tools import\
  --type 'SampleData[PairedEndSequencesWithQuality]'\
  --input-path manifest.csv\
  --output-path demux.qza\
  --input-format PairedEndFastqManifestPhred33
   ## the correct extension is automatically added for the output by qiime.
~~~

## Quality Control
Now we want to look at the quality profile of our reads.  Our goal is to determine how much we should truncate the reads before the paired end reads are joined.  This will depend on the length of our amplicon, and the quality of the reads.
~~~bash
qiime demux summarize\
   --i-data demux.qza\
   --o-visualization demux
~~~
[demux.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2Fmkhsq5pexllos38%2Fdemux.qzv%3Fdl%3D1)

When looking we want to answer these questions:

How much total sequence do we need to preserve an sufficient overlap to merge the paired end reads?

How much poor quality sequence can we truncate before trying to merge?

In this case we know our amplicons are about 390 bp long, and we want to preserve approximately 50 bp combined overlap.  So our target is to retain ~450 bp of total sequence from the two reads.  450 bp/2 = 225 bp but looking at the demux.qzv, the forward reads seem to be higher quality than the reverse, so let's retain more of the forward and less of the reverse.

## Denoising
We're now ready to denoise our data. Through qiime we will be using the program DADA2, the goal is to take our imperfectly sequenced reads, and recover the "real" sequence composition of the sample that went into the sequencer.
DADA2 does this by learning the error rates for each transition between bases at each quality score.  It then assumes that all of the sequences are errors off the same original sequence.  Then using the error rates it calculates the likelihood of each sequence arising.  Sequences with a likelihood falling below a threshold are split off into their own groups and the algorithm is iteratively applied.  Because of the error model we should only run samples which were sequenced together through dada2 together, as different runs may have different error profiles.  We can merge multiple runs together after dada2.  During this process dada2 also merges paired end reads, and checks for chimeric sequences.
~~~bash
qiime dada2 denoise-paired\
  --i-demultiplexed-seqs demux.qza\
  --p-trim-left-f 20  --p-trim-left-r 16\
  --p-trunc-len-f 290 --p-trunc-len-r 270\
  --p-n-threads 72\
  --o-denoising-stats dns\
  --o-table table\
  --o-representative-sequences rep-seqs
~~~
Now lets visualize the results of Dada2.
~~~bash
## Metadata on denoising
qiime metadata tabulate\
   --m-input-file dns.qza\
   --o-visualization dns
## Unique sequences accross all samples
qiime feature-table tabulate-seqs\
   --i-data rep-seqs.qza\
   --o-visualization rep-seqs
## Table of per-sample sequence counts
qiime feature-table summarize\
   --i-table table.qza\
   --m-sample-metadata-file mdat.tsv\
   --o-visualization table
~~~
[dns.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2F826xzk6umkcbyv4%2Fdns.qzv%3Fdl%3D1)

[rep-seqs.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2Fblv9535pizui7dj%2Frep-seqs.qzv%3Fdl%3D1)

[table.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2Fowp855n17ex7cx9%2Ftable.qzv%3Fdl%3D1)

Looking at dns.qzv first we can see how many sequences passed muster for each sample at each step performed by dada2.  Here we are seeing great final sequence counts, and most of the sequences being filtered in the initial quality filtering stage.  Relatively few are failing to merge, which suggests we did a good job selecting our truncation lengths.

In the table.qzv we can see some stats on our samples.  We have millions of counts spread across thousands of unique sequences and tens of samples.  We'll come back to the table.qzv file when we want to select our rarefaction depth.

In the rep-seqs.qzv we can see the sequences and the distribution of sequence lengths.  Each sequence is a link to a web-blast against the ncbi nucleotide database.
The majority of the sequences we observe are in our expected length range.
Later on we can use this to blast specific sequences we are interested in against the whole nucleotide database.

## Taxonomic Assignment
To assign taxonomy we will use a naive bayes classifier trained by the qiime2 authors on our gene region.
If we were using a different primer pair we would want to use a different method, like vsearch.

~~~bash
nohup qiime feature-classifier classify-consensus-vsearch\
    --i-query rep-seqs.qza\
    --i-reference-reads /home/share/databases/SILVA_databases/silva-138-99-seqs.qza\
    --i-reference-taxonomy /home/share/databases/SILVA_databases/silva-138-99-tax.qza\
    --p-maxaccepts 5\
    --p-query-cov 0.4\
    --p-perc-identity 0.7\
    --p-threads 16\
    --o-classification taxonomy &
~~~

Let's visualize the taxonomy in a few different ways.
~~~bash
qiime metadata tabulate\
   --m-input-file taxonomy.qza\
   --o-visualization taxonomy.qzv

qiime taxa barplot --i-table table.qza\
   --i-taxonomy taxonomy.qza\
   --o-visualization taxa-barplot\
   --m-metadata-file mdat.tsv
~~~
[taxonomy.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2F3p43mv2bp3z16pd%2Ftaxonomy.qzv%3Fdl%3D1)

[taxa-barplot.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2F6050v2wrr8jfgd8%2Ftaxa-barplot.qzv%3Fdl%3D1)

## Diversity analysis
Our next step is to look at the diversity in the sequences of these samples.
Here we will use the differences between the sequences in the sample, and metrics to quantify those differences to tell us about the diversity, richness and evenness of the sequence variants found in the samples.
In doing so we will construct a de novo phylogenetic tree, which works much better if we first remove any spurious sequences that are not actually the target region of our 16S gene.
To do that we will use our taxonomic assignments to filter out sequences that remained Unassigned, are assigned only as Bacteria or are Eukaryotes.  We should look at what we are filtering out and try and find out what it is.

~~~bash
## exact match to filter out unassigned and Bacteria
## exact because bacteria is part of many other that we want to keep.
qiime taxa filter-table\
   --i-table table.qza\
   --i-taxonomy taxonomy.qza\
   --p-exclude "Unassigned"\
   --p-mode exact\
   --o-filtered-table 16S-table

qiime feature-table filter-samples\
    --i-table 16S-table.qza\
    --m-metadata-file to_keep.tsv\
    --p-no-exclude-ids\
    --o-filtered-table filtered-16S-table

## Filter the sequences to reflect the new table.
qiime feature-table filter-seqs\
   --i-table filtered-16S-table.qza\
   --i-data rep-seqs.qza\
   --o-filtered-data filtered-16S-rep-seqs

qiime feature-table tabulate-seqs\
   --i-data filtered-16S-rep-seqs.qza\
   --o-visualization bacteria-rep-seqs
~~~

Now that we have only our target region we can create the de novo phylogenetic tree.
We'll use the default qiime2 pipeline because it is quick and easy to run, while providing good results.
This pipeline first performs a multi sequence alignment with mafft, this alignment would be significantly worse if we had not removed the non target sequences.
It then masks highly variable parts of the sequence as they add noise to the tree.
It then uses FastTree to create an unrooted phylogenetic tree which is then midpoint rooted.

~~~bash
qiime phylogeny align-to-tree-mafft-fasttree\
   --i-sequences filtered-16S-rep-seqs.qza\
   --o-alignment aligned-rep-seqs.qza\
   --o-masked-alignment masked-aligned-rep-seqs.qza\
   --o-tree unrooted-tree.qza\
   --o-rooted-tree rooted-tree.qza\
   --p-n-threads 16
~~~

Now we can look at the [tree we created on iToL](https://itol.embl.de/tree/20922221311082651562009639).
And for reference here is [the tree if we had not filtered it](https://itol.embl.de/tree/209222213110293351562082149).
We can see that the filtering upped the contrast between different groups.

Now we are ready to run some diversity analysis!
We are going to start by running qiimes core phylogenetic pipeline, this will take into account the relationships between sequences, as represented by our phylogenetic tree.
It will calculate a few key metrics for us, faiths-pd a measure of phylogenetic diversity, evenness, a measure of evenness and several beta statistics, like weighted and unweighted unifracs.
To do these comparisons we need to make our samples comparable to each other.
The way this is generally done is to rarefy the samples to the same sampling depth.
We can use the bacteria-table.qzv we made earlier to inform this decision.
We want to balance setting as high as possible of a rarefaction depth to preserve as many reads as possible, while setting it low enough to preserve as many samples as possible.

~~~bash
qiime diversity core-metrics-phylogenetic\
   --i-phylogeny rooted-tree.qza\
   --i-table filtered-16S-table.qza\
   --p-sampling-depth 16950\
   --m-metadata-file mdat.tsv\
   --output-dir core-metrics-results
~~~

[weighted_unifrac_emperor.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2F7vka91y2qvexljw%2Fweighted_unifrac_emperor.qzv%3Fdl%3D1)

We can see that the strains separate well, which implies that we should be able to find some separating distances in our data.

Lets start looking for those differences by looking at differences in diversity as a whole.
For numeric metadata categories we can plot our favorite metrics with the value of that metadata.

~~~bash
qiime diversity alpha-correlation\
   --i-alpha-diversity core-metrics-results/faith_pd_vector.qza\
   --m-metadata-file mdat.tsv\
   --o-visualization core-metrics-results/faith-alpha-correlation
~~~
[faith-alpha-correlation.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2F363h8ohpq8uyixw%2Ffaith-alpha-correlation.qzv%3Fdl%3D1)

Then for the categorical metadata catagories we can plot some box and whisker plots.
~~~bash
qiime diversity alpha-group-significance\
   --i-alpha-diversity core-metrics-results/faith_pd_vector.qza\
   --m-metadata-file mdat.tsv\
   --o-visualization core-metrics-results/faith-group-significance
~~~

[faith-group-significance.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2Fqlwmijytz3ivx1k%2Ffaith-group-significance.qzv%3Fdl%3D1)

## Differential Abundance Analysis
Now lets combine the taxonomy with the diversity analysis to see if there are related groups of organisms that are differentially abundant groups within the samples.
We'll start by combining our table and tree into a hierarchy and set of balances.
Balances are the weighted log ratios of sets of features for samples.
And we will be looking for significant differences in the balances between groups of samples.

Songbird is currently running behind Qiime2, so for now we need to travel back to the distant past, of June 2020...

~~~bash
source activate qiime2-2020.6

qiime songbird multinomial\
  --i-table ./filtered-16S-table.qza\
  --p-formula "Strain"\
  --p-epochs 1000\
  --p-num-random-test-examples 5\
  --p-differential-prior 0.5\
  --p-summary-interval 1\
  --o-differentials differentials\
  --o-regression-stats regression-stats\
  --m-metadata-file mdat.tsv\
  --o-regression-biplot regression-biplot

qiime songbird multinomial\
  --i-table ./filtered-16S-table.qza\
  --p-formula "1"\
  --p-epochs 1000\
  --p-num-random-test-examples 5\
  --p-differential-prior 0.5\
  --p-summary-interval 1\
  --o-differentials null-differentials\
  --o-regression-stats null-regression-stats\
  --m-metadata-file mdat.tsv\
  --o-regression-biplot null-regression-biplot

qiime songbird summarize-paired\
  --i-baseline-stats null-regression-stats.qza\
  --i-regression-stats regression-stats.qza\
  --o-visualization regression-summary
~~~

[regression-summary.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2Fepzlw713vysapgj%2Fregression-summary.qzv%3Fdl%3D1)


~~~bash
qiime qurro differential-plot\
    --i-ranks differentials.qza\
    --i-table filtered-16S-table.qza\
    --m-sample-metadata-file mdat.tsv\
    --m-feature-metadata-file  taxonomy.qza\
    --o-visualization Strain-qurro
~~~

[Strain-qurro.qzv](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fdl.dropbox.com%2Fs%2Ffcydarlze8auwrw%2FStrain-qurro.qzv%3Fdl%3D1)
