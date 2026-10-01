#!/bin/bash
# Prepare the native-space retinotopic maps for neuropythy's register_retinotopy.
#
# For each subject, hemisphere and group (empirical, predicted), this script
#   1. transforms the polar angle map to the neuropythy convention
#      (LH: 0-180 referring to UVM -> RHM -> LVM) and saves it as ..._polarAngle_neuropythy;
#   2. converts polar angle, eccentricity and pRF size to .mgz;
#   3. converts the weight maps to .mgz (the empirical variance explained for the empirical
#      group, the mean and ones dummy weights for the predicted group).
#
# The predicted maps are read under the toolbox file names, whose tag depends on the map:
# polar angle and eccentricity come from the visual coordinate model (visualCoord-model) and
# pRF size from its own model (pRFsize-model). The outputs drop the tag, so
# hcp_neuropythy_execute.py finds <sub>.predicted_<map>.<hemi>.native.func.mgz.
#
# Usage:
#   hcp_neuropythy_transforms.sh -s <subjects dir> -r <validation repo> [-i "<sub> <sub> ..."] [-n]
#     -i   subjects to process, space separated (default: every subject directory)
#     -n   dry run, print the commands without running them
ml freesurfer/7.3.2

subjects_dir=""
path_to_validation_repo=""
subject_ids=""
dry_run=0
while getopts s:r:i:n flag
do
    case "${flag}" in
        s) subjects_dir=${OPTARG};;
        r) path_to_validation_repo=${OPTARG};;
        i) subject_ids=${OPTARG};;
        n) dry_run=1;;
        ?)
           echo "script usage: $(basename "$0") -s <path to subs> -r <path to validation repo> [-i \"<sub> ...\"] [-n]" >&2
           exit 1;;
    esac
done

if [ -z "$subjects_dir" ] || [ -z "$path_to_validation_repo" ]; then
    echo "script usage: $(basename "$0") -s <path to subs> -r <path to validation repo> [-i \"<sub> ...\"] [-n]" >&2
    exit 1
fi

failures=0
fail() {
    echo "  FAILED: $*" >&2
    failures=$((failures + 1))
}

run() {
    if [ "$dry_run" -eq 1 ]; then
        echo "  $*"
    else
        "$@" || { fail "$1 (see error above)"; return 1; }
    fi
}

# Model tag of the toolbox file names for each predicted map
tag_for_map() {
    case "$1" in
        polarAngle|eccentricity) echo "visualCoord-model";;
        pRFsize) echo "pRFsize-model";;
        *) echo "unknown map: $1" >&2; return 1;;
    esac
}

# list of subjects
if [ -n "$subject_ids" ]; then
    read -r -a file_names <<< "$subject_ids"
else
    file_names=()
    for d in "$subjects_dir"/* ; do
        sub=$(basename "$d")
        if [ -d "$d" ] && [ "$sub" != "fsaverage" ] && [[ "$sub" != .* ]] && [[ "$sub" != processed_* ]] && [ "$sub" != "logs" ]; then
            file_names+=("$sub")
        fi
    done
fi

# define hemispheres
hemisphere=("lh" "rh")

# define groups
groups=("empirical" "predicted")

for sub in "${file_names[@]}"; do
    path_sub="$subjects_dir/$sub/deepRetinotopy"
    if [ ! -d "$path_sub" ]; then
        fail "$sub: $path_sub does not exist"
        continue
    fi

    for hem in "${hemisphere[@]}"; do
        for group in "${groups[@]}"; do
            echo "Preparing $group maps of the $hem hemisphere in $sub..."

            # Input maps of this group. The predicted maps carry the toolbox model tag.
            declare -A input_map
            for param in polarAngle eccentricity pRFsize; do
                if [ "$group" = "predicted" ]; then
                    input_map[$param]="$path_sub/$sub.predicted_${param}_$(tag_for_map $param).$hem.native.func.gii"
                else
                    input_map[$param]="$path_sub/$sub.empirical_${param}.$hem.native.func.gii"
                fi
            done
            if [ "$group" = "predicted" ]; then
                weights=("$path_sub/$sub.mean_variance_explained.$hem.native.func.gii"
                         "$path_sub/$sub.ones_variance_explained.$hem.native.func.gii")
            else
                weights=("$path_sub/$sub.empirical_variance_explained.$hem.native.func.gii")
            fi

            missing=0
            for f in "${input_map[@]}" "${weights[@]}"; do
                if [ ! -f "$f" ]; then
                    echo "  missing input: $f" >&2
                    missing=1
                fi
            done
            if [ "$missing" -eq 1 ]; then
                fail "$sub $hem $group: missing inputs, skipped"
                continue
            fi

            # 1. Transform the polar angle map to the neuropythy convention
            path_to_save="$path_sub/$sub.${group}_polarAngle_neuropythy.$hem.native.func.gii"
            run python "$path_to_validation_repo/functions/preprocess.py" transform_polarangle_to_benson14 \
                --path_to_use "${input_map[polarAngle]}" --path_to_save "$path_to_save" --hemisphere "$hem" || continue

            # 2. Convert the maps to .mgz, without the model tag
            run mri_convert "$path_to_save" "$path_sub/$sub.${group}_polarAngle_neuropythy.$hem.native.func.mgz"
            run mri_convert "${input_map[eccentricity]}" "$path_sub/$sub.${group}_eccentricity.$hem.native.func.mgz"
            run mri_convert "${input_map[pRFsize]}" "$path_sub/$sub.${group}_pRFsize.$hem.native.func.mgz"

            # 3. Convert the weight maps to .mgz
            if [ "$group" = "predicted" ]; then
                for dummy_weight in mean ones; do
                    run mri_convert "$path_sub/$sub.${dummy_weight}_variance_explained.$hem.native.func.gii" \
                        "$path_sub/$sub.predicted_${dummy_weight}_weight.$hem.native.func.mgz"
                done
            else
                run mri_convert "${weights[0]}" "$path_sub/$sub.empirical_weight.$hem.native.func.mgz"
            fi
        done
    done
done

echo "Done: ${#file_names[@]} subjects, $failures failures."
echo "Verify with: ls $subjects_dir/<sub>/deepRetinotopy/*.native.func.mgz"
exit $((failures > 0))
