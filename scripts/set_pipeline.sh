#!/bin/sh

seed=12345

scripts/sample_set_for_saturation.sh $seed &&
scripts/construct_set_bins_for_saturation.sh $seed &&
scripts/call_set_peaks_for_saturation.sh $seed &&
scripts/call_set_binding_sites_for_saturation.sh $seed
