#!/bin/sh

seed=12345

scripts/sample_pet_for_saturation.sh $seed &&
scripts/construct_pet_bins_for_saturation.sh $seed &&
scripts/call_pet_peaks_for_saturation.sh $seed &&
scripts/call_pet_binding_sites_for_saturation.sh $seed
