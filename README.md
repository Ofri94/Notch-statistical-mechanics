# Notch-statistical-mechanics
Predicting the transcriptional response and sensitivity of Notch related Genetic elements with different architectures

Scripts can be devided to three groups:
Site composition modulations, Cryptic sites and distance.

nth_csl_mth_sps, nth_csl_sites, nth_csl_mth_sps_plotter are the first group.
nth_csl_mth_lesser_sps is second group.
Three_sites_1_Fr, Three_sites_2_Fr, two_csls_Fr are the third group.

The first group calculates Transcriptional response and sensitivity for any combination of CSL and SPS sites.
The second group calculates the same as the first, with the addition of 1 cryptic SPS site.
The third can calculate variations of only 2 or 3 sites (according to the name), with the ability to increase distance between sites.

n - number of CSL sites
m - number of SPS sites
C - cooperativity factor for SPS sites
phi - the factor by which cryptic sites are lessened in their affinity
Fr - a measure between 0-1 indicating distances affects on the repressor. 1 is no distance (fully operated repressor) and 0 is "infinite" distance (repressor doesn't work).
