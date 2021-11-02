#! /bin/bash

#################################################
#################################################
#####                                        ####
##### Amplicon analysis using QIIME2 & DADA2 ####
#####         Several libraries              ####
#####                                        ####
#################################################
#################################################


#activate qiime (later in script)
conda activate qiime2-2021.2
cd /home/mfbeuq/Amplicon/V4a_spat

#Metadata path (requirements: https://docs.qiime2.org/2020.2/tutorials/metadata/)
map=/home/mfbeuq/Amplicon/V4a_spat/scripts_maps/map_file_v4a_spat.tsv


#### create manifest file (2020.02 version) ####
for j in Lauf1_*
do
   cd $j/demux-all
    #header line
    echo "sample-id,forward-absolute-filepath,reverse-absolute-filepath" > MANIFEST.tsv
    #write info for forward reads into MANIFEST
    for i in *_R1_*.fastq.gz
    do
        #echo $i 
        var2=$(echo "$i" | (grep -o -a  "S[0-9]*" | head -1))
        #make corresponding R2-file name
        var3=$(echo "$i" | (sed "s/\(^S[0-9]*_.*_R\)[1-2]\(_.*fastq.gz$\)/\12\2/" ))
        echo $var2,"/home/mfbeuq/Amplicon/V4c/$j/demux-all/$i","/home/mfbeuq/Amplicon/V4c/$j/demux-all/$var3" >> MANIFEST.tsv
    done
    cat MANIFEST.tsv | sed 's/\,/\t/g' > MANIFEST ; rm MANIFEST.tsv ; cd ../
    cd ../
done


#import as whole artifact
for i in Lauf1_*
do
   cd $i
   qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /home/mfbeuq/Amplicon/V4a_spat/$i/demux-all/MANIFEST \
  --output-path reads-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2  
cd ../
done


#primer trimming
#sequences are
#799f: AACMGGATTAGATACCCKG
#1193r: ACGTCATCCCCACCTTCC
for i in Lauf1_*
do
   cd $i
   qiime cutadapt trim-paired \
    --i-demultiplexed-sequences reads-demux.qza \
    --p-front-r ACGTCATCCCCACCTTCC \
    --p-front-f AACMGGATTAGATACCCKG \
    --p-error-rate 0.15\
    --p-cores 12 \
    --o-trimmed-sequences trim-seqs.qza \
    --verbose
    
    qiime demux summarize \
    --i-data trim-seqs.qza \
    --o-visualization trim-seqs.qzv    
    
    cd ../
done

#visualise: check for decreasing quality
qiime tools view Lauf1_1/trim-seqs.qzv
qiime tools view Lauf1_2/trim-seqs.qzv


#denoising using dada2 - takes long to run!
#trun-len --> decide based on qualities from above
#change if required (truncate at q-score of > 20)
for i in Lauf1_*
do
   cd $i
   qiime dada2 denoise-paired \
    --i-demultiplexed-seqs trim-seqs.qza \
    --p-trunc-len-f 220 \
    --p-trunc-len-r 220 \
    --p-n-threads 0 \
    --o-table denoise.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats denoising-stats.qza \
    --verbose
    
#create visualisations
#first: overall denoising output ('OTU-table')
#incl. metadata

qiime feature-table summarize \
    --i-table denoise.qza \
    --m-sample-metadata-file $map \
    --o-visualization denoise.qzv \
    
    cd ../
done

#visualise
qiime tools view Lauf1_1/denoise.qzv
qiime tools view Lauf1_2/denoise.qzv


#third: denoising stats
for i in Lauf1_*
do
   cd $i
   qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv 
  cd ../
done

#visualise
qiime tools view Lauf1_1/denoising-stats.qzv
qiime tools view Lauf1_2/denoising-stats.qzv


#####################MERGING

mkdir QIIME_output ; cd QIIME_output

#if you have several runs, you should be able to combine them now, using the following commands

qiime feature-table merge \
  --i-tables ../Lauf1_1/denoise.qza \
  --i-tables ../Lauf1_2/denoise.qza \
  --p-overlap-method sum \
  --o-merged-table table.qza

qiime feature-table merge-seqs \
  --i-data ../Lauf1_1/rep-seqs.qza \
  --i-data ../Lauf1_2/rep-seqs.qza \
  --o-merged-data rep-seqs.qza

#FeatureTable and FeatureData summaries
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file $map
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
qiime tools view table.qzv
qiime tools view rep-seqs.qzv



### QIIME analysis

mkdir filtered_tables


#Filter out rare ASVs: features that appear less then 20 times and not in least 5 samples are thrown out
qiime feature-table filter-features \
    --i-table filtered_tables/table.qza \
    --p-min-frequency 20 \
    --p-min-samples 5 \
    --o-filtered-table filtered_tables/table-filtered.qza
#Filter samples with less than 10000 reads
qiime feature-table filter-samples \
    --i-table filtered_tables/table-filtered.qza \
    --p-min-frequency 10000 \
    --o-filtered-table filtered_tables/table-pruned-filtered.qza 
#Remove Singletons
qiime feature-table filter-features \
    --i-table filtered_tables/table-pruned-filtered.qza  \
    --p-min-samples 2 \
    --o-filtered-table filtered_tables/all-table-pruned-filtered.qza 
#Create a visualization
qiime feature-table summarize \
    --i-table filtered_tables/all-table-pruned-filtered.qza  \
    --o-visualization filtered_tables/all-table-pruned-filtered.qzv \
    --m-sample-metadata-file $map 
    
#create phylogenetic trees:
#sequence alignment, mafft ist die software
qiime alignment mafft \
    --i-sequences rep-seqs.qza \
    --o-alignment aligned-rep-seqs.qza \
    --verbose
qiime alignment mask \
    --i-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --verbose

#build unrooted tree
qiime phylogeny fasttree \
    --i-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --verbose
  
#rooted tree, using longest root
qiime phylogeny midpoint-root \
    --i-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza \
    --verbose

#Taxonomic assignment using a custom built QIIME naive Bayesian classifier (accessible at ...........................)
qiime feature-classifier classify-sklearn \
    --i-classifier /home/mfbeuq/Amplicon/bin/training-feature-classifiers/138/silva-138-ssu-nr99-799f-1193r-classifier.qza \
    --i-reads rep-seqs.qza \
    --o-classification taxonomy.qza \
    --p-n-jobs 2 \
    --p-confidence 0.75 \
    --verbose
qiime metadata tabulate \
    --m-input-file taxonomy.qza \
    --o-visualization taxonomy.qzv

    
    
#Filter the data table into three subsets, each belonging to one of the loosely associated (L), tightly associated root microbiota (T) and the bulk soil microbiota (B)
for comp in "L" "T" "B"; do
    qiime feature-table filter-samples \
        --i-table all-table-pruned-filtered.qza \
        --m-metadata-file $map \
        --p-where "[compartment]='"$comp"'" \
        --o-filtered-table filtered_tables/"$comp"-table-pruned-filtered.qza
    #Create a visualization for all three
    qiime feature-table summarize \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --o-visualization filtered_tables/"$comp"-table-pruned-filtered.qzv \
        --m-sample-metadata-file $map
     
    #Filter the each of the three data tables into four subsets, each belonging to one of the four individual trees
    for tree in "T1" "T2" "T3" "T4" ; do
        qiime feature-table filter-samples \
            --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
            --m-metadata-file $map \
            --p-where "[tree]='"$tree"'" \
            --o-filtered-table filtered_tables/tree/"$comp"-"$tree"-table.qza
done  


#collapse the original data table and all three sub datasets to different taxonomic levels
mkdir collapsed_taxa
for comp in "L" "T" "B" "all"; do
    #phylum 
    qiime taxa collapse \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --i-taxonomy taxonomy.qza \
        --p-level 2 \
        --o-collapsed-table collapsed_taxa/"$comp"-table-pruned-filtered-phylum.qza
    #class
    qiime taxa collapse \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --i-taxonomy taxonomy.qza \
        --p-level 3 \
        --o-collapsed-table collapsed_taxa/"$comp"-table-pruned-filtered-class.qza
    #order
    qiime taxa collapse \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --i-taxonomy taxonomy.qza \
        --p-level 4 \
        --o-collapsed-table collapsed_taxa/"$comp"-table-pruned-filtered-order.qza
    #family
    qiime taxa collapse \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --i-taxonomy taxonomy.qza \
        --p-level 5 \
        --o-collapsed-table collapsed_taxa/"$comp"-table-pruned-filtered-family.qza
    #genus
    qiime taxa collapse \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --i-taxonomy taxonomy.qza \
        --p-level 6 \
        --o-collapsed-table collapsed_taxa/"$comp"-table-pruned-filtered-genus.qza
    #species
    qiime taxa collapse \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --i-taxonomy taxonomy.qza \
        --p-level 7 \
        --o-collapsed-table collapsed_taxa/"$comp"-table-pruned-filtered-species.qza
done



#create core-metrics including Shannon's diversity index for each table
mkdir core-metrics
for comp in "L" "T" "B" "all"; do
    qiime diversity core-metrics-phylogenetic \
        --i-phylogeny rooted-tree.qza \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --p-sampling-depth none \
        --m-metadata-file $map \
        --p-n-jobs-or-threads auto \
        --output-dir core-metrics/"$comp"-core-metrics-results
done

#extract the Shannon's diversity calculation into a readable .tsv file
for i in core-metrics/*metrics-results; do
    cd $i
    qiime diversity alpha-group-significance \
        --i-alpha-diversity shannon_vector.qza \
        --m-metadata-file $map \
        --o-visualization shannon_vector.qzv
    unzip -p shannon_vector.qzv */data/metadata.tsv >shannon_metadata.tsv
    cd ../../
done


#For beta-diversity: Robust Aitchison PCA Beta Diversity with DEICODE
for comp in "L" "T" "B"; do
    mkdir DEICODE/"$comp"
    qiime deicode auto-rpca \
        --i-table filtered_tables/"$comp"-table-pruned-filtered.qza \
        --p-min-feature-count 10 \
        --p-min-sample-count 500 \
        --o-biplot DEICODE/"$comp"/RA_ordination.qza \
        --o-distance-matrix DEICODE/"$comp"/RA_distance.qza
    qiime emperor biplot \
        --i-biplot DEICODE/"$comp"/RA_ordination.qza \
        --m-sample-metadata-file $map \
        --m-feature-metadata-file taxonomy.qza \
        --o-visualization DEICODE/"$comp"/RA_biplot.qzv \
        --p-number-of-features 5
    for var1 in tree root_qdt; do
        qiime diversity beta-group-significance \
            --i-distance-matrix DEICODE/"$comp"/RA_distance.qza \
            --m-metadata-file $map \
            --m-metadata-column "$var1" \
            --p-method permdisp \
            --o-visualization DEICODE/"$comp"/"$var1"_permdisp_significance.qzv
done; done



###Export a biom file for further analysis in R
mkdir ../phyloseq

#Export OTU table
qiime tools export \
    --input-path filtered_tables/all-table-pruned-filtered.qza \
    --output-path ../phyloseq/
 

# OTU tables exports as feature-table.biom so convert to .tsv
biom convert \
  -i ../phyloseq/feature-table.biom \
  -o ../phyloseq/feature-table.tsv \
  --to-tsv
  
#Change #OTUID to OTUID
var="OTUID"
sed -i "2s/^#OTU ID/$var/" ../phyloseq/feature-table.tsv

#Export taxonomy table
qiime tools export \
  --input-path taxonomy.qza \
  --output-path ../phyloseq

#Change "feature ID" to "OTUID"
var="OTUID \t taxonomy \t confidence"
sed -i "1s/.*/$var/" ../phyloseq/taxonomy.tsv
 
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path ../phyloseq
