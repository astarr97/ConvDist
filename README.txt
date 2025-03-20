IMPORTANT: You should not have a species appear more than once as the primary related or outgroup species in a particular set of species you want to compare.  For example,
for the horse + rhino clade and cetaceans, you could reasonably have Tragulus_javanicus be the primary outgroup.  However, when you reach step(4) you will need to repull variants
and Tragulus_javanicus will appear twice in the command for halSNPs.  halSNPs removes duplicate species, so this well then cause problems in remHead_callConvDiv.py and you would
need to modify remHead_callConvDiv.py to add the second Tragulus_javanicus (I had to do this repeatedly).  Instead, it is much easier to just make Bos_taurus the primary outgroup for cetaceans 
and keep Tragulus_javanicus as the primary outgroup for  the horse + rhino clade.  This will result in nearly identical output.

IMPORTANT: You may want the PhyloP species to not be any of the focal species. That works fine, you just need to modify the string of species at step 4 to include a burn in of 
the PhyloP species repeated three times (e.g. "Mus_musculus,Mus_musculus,Mus_musculus;...) where ... represents the remaining species triplets.  It is vital that the PhyloP species not
appear as one of the primary related or outgroup species for any comparison for the same reason as above (no duplicate species for halSNPs). You will then need to delete all the 
Mus_musculus containing segments of code from the resulting .sh files BUT be careful to retain the gzip lines for the FiltConv and FiltPoly files that back it up.

IMPORTANT: You will need to modify remHead_callConvDiv.py to not remove the thing after the first "." from the scaffold name for Nycticebus_pygmaeus due to it having a weird scaffold naming scheme.

ComputingPhyloP and ComputingVariants just contains repeats of code in the main ConvDist folder to help with clarity.

All code is written in python, it also requires cactus and phast.  

This folder contains the pipeline to do various tests for convergent evolution, as well as decelerated evolution/positive selection for specific clades.
The goal of the pipeline is to:
    1. Compute PhyloP scores and species support with an arbitrary set of masked genomes.
    2. Pull the genetic variants associated with specific clades to test for accelerated evolution.
    2.5 Post-process the PhyloP scores
    3. Postprocess this information and intersect it all together as well as perform various other operations to prepare to do the tests.
    4. Intersect and identify convergent/divergent sites
    5. Actually do the tests
    
The first two steps take the vast majority of the compute power and must be run in parallel across many different cores to complete in a reasonable amount of time.
Notably, you can use precomputed PhyloP scores without masking any genomes without any issues.
The last two steps take longer, but can be run on one or a few cores with one part of step 3 requiring access to more compute.

The general setup is that for some convergently evolved trait, you need to have 3 clades (or single species depending on the trait).
We will use marine/aquatic adaptation as an example.
The clade with the trait is referred to as the focal clade.  The most related clade that does not have the trait is the related clade.  
The outgroup that does not have the trait and is equidistant from the focal and related clades is referred to as the outgroup.
For example, cetaceans (whales and dolphins) could be the focal clade, bovidae (cows, sheep, goats, antelope) the related clade, and pigs the outgroup.
Alternatively, sea otters could be the focal clade, weasels the related clade, and honey badgers the outgroup.
In this case, you might be unsure whether you want to include river otters along with sea otters.  You have the option to include it in the computationally intensive steps, but ignore after (see below)

In general, we require that all bases in the clade be the same to be considered fixed and used for the NearestDist/ConvVsDiv tests.
If they are at all variable in the clade, we call it polymorphic.

STEP 1
STEP 1 is done for all the focal species at once and separately for all the related species at once, referenced to the same genome
To compute the PhyloP scores:
IMPORTANT: In general, you must mask entire clades.  For example, for seals and walruses (pinnipeds), sea and river otters are in the related clade.
When computing PhyloP scores for the rel species, you must mask sea and river otters, in addition to the badgers, weasels, skunks, etc. as they are all in the related clade.
You must do this even if you are not including them in the variant pulling for seals and walruses!

Some notes before we begin.  Due to an unresolved bug in the phast code, phylop will occasionally error on single bases in large files, leading to no output for the entire file being computed on.
To combat this issue, we split the genome into 10 megabase chunks that are computed in parallel.  Each 10 megabase chunk is further split into 1 megabase chunks that get turned into a deduplicated, masked maf file.
Each of those 1 megabase chunks is further split into 1,000 base pair chunks.  These are what we compute PhyloP on.
If one of those computations errors, there is a system to catch that error, split the 1,000 base pair file into single maf blocks, and compute on those.
This way, we only lose 1 maf block to this bug.

The shell script that illustrates how to do this computation is available in Master_Shell_Scripts/compute_phylop_Aquatic.sh
It relies on ScriptsThatMakeScripts/make_scripts_mask_phylop.py to make the scripts
This then relies on check_error_write.py to search for errors and do the necessary recomputation.  This will only handle errors that start with "ERROR" (no quotes), so it is worth running check_errors.py
to see if there are any other errors that need to be handled.

Details:
    1. Define the species you want to mask, the species you want the PhyloP scores to be referenced to (should ideally be one with a contiguous assembly and good annotation), and the output folder name.
    2. Get the chromosome names and sizes for the reference genome with samtools (key input to the script that makes scripts)
    3. Run the script that makes scripts.  This will create folders for each run (run1-runN), an All folder where it will be output, 
    an Errored folder where the 1,000 base pair mafs that errored will be moved (to check error handling is working), and a gtf folder for the various gtfs needed
    
Results will be moved to the All folder as they come in.  You should check that there are as many bed/tsv files in the All directory as there are run folders in the parent directory.
You should also run check_error.py to see if there are any errors that do not start with ERROR as we will not catch those and they may need to be corrected
Overall, this will probably take a little over 24 hours to finish assuming there are enough cores available to do all the computation in parallel

Finally, it is often desirable to run the same tests on the related species as these did not convergently evolve a trait and so can act as a good control.  
To do this, you can generally just rerun the same shell, but with the related species masked instead of the focal species.  An example of this is in Master_Shell_Scripts/compute_phylop_Aquatic_Flipped.sh

STEP 2
STEP 2 is done for each clade of interest separately
NOTE: How much memory you need to allocate per job (which can be changed in the write_beg python function in make_scripts_get_variants_dist.py) depends on the number of species and how closely related they are
For the cetaceans + bovidae, you need 24-32GB of memory per job
For otter, 16 (or less, not tested) is fine
To pull the variants for each set of focal, related, and outgroup species/clades:
Pulling the genetic variants and then processing them requires two python scripts and a Config file, the example file is Master_Shell_Scripts/compute_variants_Aquatic.sh
Similar to the PhyloP file, this requires a bed file of chromosome sizes (how to get it is included in the example script)
The config file (example in this directory, Config_Aquatic.txt) has four columns Focal_species  Related_species	Outgroup_species	Contigs_file
The first three are comma delimited lists of the focal speces, related species, and outgroup species.  There can be two parts to this list, separated by a semi-colon.
The first part is the set of species we want to require all have the same base for it to be considered fixed (or to differ if we want to call it polymorphic)
The part after the semi-colon contains species that will have the variants pulled and put in the tsv file, but ignored after that.
For example, you can include the river otter after the semi-colon if you want to pull the variants but not do any filtering on it.
You can check if the config file you made is valid by using the script validate_config.py

IMPORTANT: The ordering, specifically the first species in the list, is important.  These are considered the "primary" focal, related, or outgroup species.  
The primary focal species is the genome to which the output will be referenced.  The primary focal, related, and outgroup species will be saved by filter_and_tag_variants.py
and will also be used as a first pass filter.  For example, if a site has an N or NA in the primary related species, that site will not be included.
However, if another non-primary species has an N or NA but the primary related species (and potentially other related species depending on the parameters used) has an actual base, that site will be kept.
In general, the primary focal, related, and outgroup species should have as high of genome quality as possible, assuming there is any choice.

Similar to the PhyloP code, this code splits things into blocks of either 20 million bases, or 750 contigs whichever comes first.  
It then pulls the bases per 500,000 base chunk, cats everything together, and then uses filter_and_tag_variants.py to process the variants for downstream steps.  

filter_and_tag_variants.py computes a few useful statistics and outputs them (see the script for more details).  The column layout in the resulting bedfile is:
"refSequence", "refPosition1", "refPosition2", spec_focus, spec_rel, spec_out, "Simple_Derived", "Focus_FixedOrPoly", "Focus_PropBase", "Related_FixedOrPoly", "Related_PropBase", "Outgroup_FixedOrPoly", "Outgroup_PropBase", "FocusEqualRelated", "FocusEqualOutgroup", "RelatedEqualOutgroup", "FocusIsSingelton", "RelatedIsSingelton", "OutgroupIsSingelton"
The first 3 columns are the necessary bed columns (0 indexed chromosomal position)
The next 3 are the primary focal, related, and outgroup species
After that, Simple_Derived refers to whether the change is primary focal derived or primary related derived, or ambiguous, using only the 3 primary species to determine it
The FixedOrPoly and PropBase columns refer to whether a position is fixed or polymorphic in the clade of interest (always fixed if there is only one species) and the proportion of species with a valid A, C, G, or T base
FocusEqualRelated and the next two columns indicate whether there are any focal species that have a base that equals any found in the related species
Finally, IsSingelton indicates whether, if a site is polymorphic, multiple species have the minor allele (N) or if it is just one species (S)
If a site is tri-allelic and the two minor alleles are singeltons, it then outputs (MS)
This column is a potentially useful filter to guard against assembly errors that are genome-specific.  E.g. the protein-coding region of LYST in Orcinus_orca has very clear errors that can massively skew things.
For the fixed-within-clade vs polymorphic-within-clade approach, this is generally expected to make things slightly conservative, though it may make it anti-conservative for some genes.

After this finished, the bed and tsv file are moved to the All folder.

STEP 2.5
Processing the PhyloP scores:
The next step is handled by Master_Shell_Scripts/postprocess_phylop.sh.  This makes sure each line is unique, sorts, and backs up the PhyloP file after converting to bigwig.
It also converts the species support to a bed format, sorts, and backs up.

PRE STEP 3
You may need to generate a bed files of CDS for input to STEP 3.  
An example of how to do this is found in the "Mouse_Annotation_Example" folder
In general, this works by taking a gtf file, converting to bed, merging overlapping CDS for each gene, and then sorting.  
It may also require liftover to get things in the proper coordinates.
Outside of human and mouse (I use gencode for these), you can use the annotations made by TOGA (which again will need to be processed) https://genome.senckenberg.de/download/TOGA/
The contig/scaffold names used by TOGA are often not used in the 447-way hal file and you will need to use NCBI to match them

STEP 3
STEP 3 is done for each clade of interest separately
Processing the variants:
There are two shell scripts you might use to process the variants.  The first, process_non_PhyloP_Ref.sh, should be used if the focal species is not the one for which PhyloP scores are referenced.
The second, process_PhyloP_Ref.sh, should be used if the reference genome for the variants and PhyloP scores is the same. 
The main difference between the two scripts is that the non_PhyloP_Ref version does reciprocal liftover and then one more liftover to get to the right reference genome.

IMPORTANT: A note on checking you did everything correctly.  If something was not masked that should have been or vice versa, it will likely show up if you plot the 100,000-200,000 highest Conv/Poly PhyloP scores.
They should have a nearly identical distribution, notably with lots of values a little greater than or near 10.  
However, if something got messed up, then one of them (most likely Conv) will not have these very high PhyloP scores!
This is how I caught an error in the masking for the marine study.

IMPORTANT: During processing, we remove the ".1" (by removing anything after the first ".") from the contig name.  We do this to be compatible with TOGA annotations

IMORTANT: You need to set the parameters at the top of the file, as well as a different phylop_path around 50 lines from the end of the script as there is different masking done for the related clade/species.

When processing the non-reference, the strategy is to first cat everything together, make unique, sort, compress and backup.
Then we split into "Conv" (i.e. fixed within clade) and "Poly" (i.e. polymorphic within clade).  This is done with filter_for_conv.py and filter_for_poly.py.
If there is only one species, the Poly file will be empty.
The default parameters for Conv are to require the same base in all focal species, have greater than 25% of the focal species with a valid base in the alignment, and the base not to equal any of the rel or out species.
For Poly, we require at least one focal species to have a different base and that one of the alleles be equal to at least one of the rel or out species.

After this, we use reciprocal liftover to filter out sites that do not uniquely lift over.
We then lift back to the reference species, and intersect with the PhyloP, gene annotation, and species support then reformat and move to the output directory (backup will occur after doing stuff in the output directory)
We do all of this for both the focal clade/species and related clade/species so that we can check that any results from the focal clade/species do not also appear in the rel clade/species
Finally, we pull out singeltons and multi-singeltons (see step 2), lift over to the reference species genome, compress both lifted and unlifted, and backup.
Singeltons are (rarely) due to assembly errors that do not affect alignment (i.e. alignment is still correct) but lead to stretches of some genomes with a lot of singeltons.
As a result, when doing the fixed-within-clade vs poly-within-clade comparison, there is the option to remove singeltons
If there are no poly-within-clade sites, then these files will be empty

For the PhyloP reference species, we skip all the liftover and just do the intersection part.

STEP 4
STEP 4 is done with all clades/species together after STEP 3 is completed for all clades/species of interest
This is done pairwise for all species so with 4 clades/species there is 4!/(2!*2!) = 6 comparisons, with 5 clades/species it is 5!/(3!*2!) = 12 comparisons
(You can see how this might get out of hand if you have say 10 clades/species which would be 45 comparisons)

To do this, we use the python script 

The goal of STEP 4 is to identify convergent/divergent sites and determine whether the related/outgroup for both clades/species of interest in the comparison has the same base
Unfortunately, I have not figured out a good way to determine this (due to many sites being the reverse complement in one clade compared to another clade)
As a result, we just pull the variants again this time just taking the primary focal, related, and outgroup species for the two clades/species of interest in the comparison.

StatisticalTests contains the tests we can run on the output of this pipeline.  The first pair, run_fwc_pwc.py and run_fwc_pwc.sh compare the PhyloP scores of sites that are a unique
nucleotide in the focal set of species and are invariant in that clade to those that are not invariant in that clade to better understand changes in selection on the branch of interest
relative to all other branches that come after it.  You can then combine p-values across clades to test for convergence.  See the shell script and python script for implementation details.  

The second pair, run_nearest_site.py and run_nearest_site.sh directly test for convergent evolution by comparing the PhyloP score of convergent/divergent sites to the nearest site 
in which a change occurred only in one lineage.  It then combines p-values across all possible pairwise comparisons and combines the combined p-value for convergent and divergent.

Finally, make_scripts.py and get_fdr.py can be used to do permutations for run_nearest_site.py in order to check that the FDR is well-calibrated (in all cases so far it has been).
