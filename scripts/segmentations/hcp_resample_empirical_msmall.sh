#!/usr/bin/env bash
# Resample the HCP empirical maps from 32k fs_LR to each subject's native surface through the HCP
# MSMAll registration.
#
# The HCP 7T retinotopy maps are MSMAll aligned, and the manual visual area labels (hcp-annot-vc)
# were transferred to native space through that registration. The toolbox route, FreeSurfer's
# sphere.reg with the fsaverage-deformed fs_LR sphere, places the same maps a few millimetres away
# in some subjects. This script uses the subject's MSMAll native sphere instead, fetched by
# scripts/hcp_msmall_data.sh. Polar angle is resampled as cosine and sine and recombined, so the
# 0/360 wrap of the left hemisphere is never interpolated.
#
# Outputs, per subject and hemisphere, in <output dir>/<sub>/deepRetinotopy/:
#   <sub>.empirical_{polarAngle,eccentricity,pRFsize,variance_explained}.<hemi>.native.func.gii
#
# Usage:
#   hcp_resample_empirical_msmall.sh -f <freesurfer dir with the 32k maps> -m <hcp_msmall dir> \
#       -t <dir with S1200_7T_Retinotopy181.[LR].sphere.32k_fs_LR.surf.gii> -r <validation repo> \
#       -o <output subjects dir> [-i "<sub> <sub> ..."] [-n]
#     -i   subjects to process (default: every subject directory in the output dir)
#     -n   dry run, print the commands without running them
ml connectomeworkbench/1.5.0

freesurfer_dir=""
msmall_dir=""
template_dir=""
path_to_validation_repo=""
output_dir=""
subject_ids=""
dry_run=0
while getopts f:m:t:r:o:i:n flag
do
    case "${flag}" in
        f) freesurfer_dir=$(realpath "${OPTARG}");;
        m) msmall_dir=$(realpath "${OPTARG}");;
        t) template_dir=$(realpath "${OPTARG}");;
        r) path_to_validation_repo=$(realpath "${OPTARG}");;
        o) output_dir=$(realpath -m "${OPTARG}");;
        i) subject_ids=${OPTARG};;
        n) dry_run=1;;
        ?)
           echo "script usage: $(basename "$0") -f <freesurfer dir> -m <hcp_msmall dir> -t <template dir> -r <validation repo> -o <output dir> [-i \"<sub> ...\"] [-n]" >&2
           exit 1;;
    esac
done

if [ -z "$freesurfer_dir" ] || [ -z "$msmall_dir" ] || [ -z "$template_dir" ] || [ -z "$path_to_validation_repo" ] || [ -z "$output_dir" ]; then
    echo "script usage: $(basename "$0") -f <freesurfer dir> -m <hcp_msmall dir> -t <template dir> -r <validation repo> -o <output dir> [-i \"<sub> ...\"] [-n]" >&2
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

# list of subjects
if [ -n "$subject_ids" ]; then
    read -r -a file_names <<< "$subject_ids"
else
    file_names=()
    for d in "$output_dir"/* ; do
        [ -d "$d" ] && file_names+=("$(basename "$d")")
    done
fi

# define hemispheres
hemisphere=("lh" "rh")

for sub in "${file_names[@]}"; do
    path_sub="$output_dir/$sub/deepRetinotopy"
    for hem in "${hemisphere[@]}"; do
        if [ "$hem" == "lh" ]; then hemi=L; else hemi=R; fi
        echo "Resampling empirical maps of the $hem hemisphere in $sub through the MSMAll registration..."

        sphere_32k="$template_dir/S1200_7T_Retinotopy181.$hemi.sphere.32k_fs_LR.surf.gii"
        sphere_native="$msmall_dir/$sub/$sub.$hemi.sphere.MSMAll.native.surf.gii"
        area_32k="$msmall_dir/$sub/$sub.$hemi.midthickness_MSMAll.32k_fs_LR.surf.gii"
        area_native="$msmall_dir/$sub/$sub.$hemi.midthickness.native.surf.gii"
        missing=0
        for f in "$sphere_32k" "$sphere_native" "$area_32k" "$area_native" \
                 "$freesurfer_dir/$sub/surf/$sub.fs_empirical_polarAngle_$hem.func.gii"; do
            if [ ! -f "$f" ]; then
                echo "  missing input: $f" >&2
                missing=1
            fi
        done
        if [ "$missing" -eq 1 ]; then
            fail "$sub $hem: missing inputs, skipped"
            continue
        fi
        [ "$dry_run" -eq 1 ] || mkdir -p "$path_sub"

        # 1. Polar angle, resampled as cosine and sine components and recombined afterwards
        prefix="$path_sub/$sub.empirical_polarAngle_${hem}.32k"
        run python "$path_to_validation_repo/functions/preprocess.py" polarangle_to_components \
            --path_to_use "$freesurfer_dir/$sub/surf/$sub.fs_empirical_polarAngle_$hem.func.gii" --path_prefix "$prefix" || continue
        for component in cos sin; do
            run wb_command -metric-resample "${prefix}_${component}.func.gii" \
                "$sphere_32k" \
                "$sphere_native" ADAP_BARY_AREA "$path_sub/$sub.empirical_polarAngle_${component}.$hem.native.func.gii" \
                -area-surfs "$area_32k" "$area_native"
        done
        run python "$path_to_validation_repo/functions/preprocess.py" components_to_polarangle \
            --cos_path "$path_sub/$sub.empirical_polarAngle_cos.$hem.native.func.gii" \
            --sin_path "$path_sub/$sub.empirical_polarAngle_sin.$hem.native.func.gii" \
            --path_to_save "$path_sub/$sub.empirical_polarAngle.$hem.native.func.gii"
        [ "$dry_run" -eq 1 ] || rm -f "${prefix}"_{cos,sin}.func.gii "$path_sub/$sub.empirical_polarAngle_"{cos,sin}".$hem.native.func.gii"

        # 2. Eccentricity, pRF size and variance explained, resampled directly
        for map in eccentricity pRFsize variance_explained; do
            run wb_command -metric-resample "$freesurfer_dir/$sub/surf/$sub.fs_empirical_${map}_$hem.func.gii" \
                "$sphere_32k" \
                "$sphere_native" ADAP_BARY_AREA "$path_sub/$sub.empirical_${map}.$hem.native.func.gii" \
                -area-surfs "$area_32k" "$area_native"
        done
    done
done

echo "Done: ${#file_names[@]} subjects, $failures failures."
echo "Verify with: ls $output_dir/<sub>/deepRetinotopy/<sub>.empirical_*.native.func.gii"
exit $((failures > 0))
