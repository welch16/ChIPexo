#!/bin/sh

seed=12345

scripts/sample_exo_for_saturation.sh $seed &&
scripts/construct_exo_bins_for_saturation.sh $seed &&
scripts/call_exo_peaks_for_saturation.sh $seed &&
scripts/call_exo_binding_sites_for_saturation.sh $seed
