#!/usr/bin/env bash
# Download the HCP MSMAll registration files needed to resample the 32k fs_LR retinotopy maps to
# each subject's native surface the same way the HCP pipelines and neuropythy do.
#
# The manual visual area labels (hcp-annot-vc) were transferred to native space through the MSMAll
# registration, whereas our native maps go through FreeSurfer's sphere.reg and the fsaverage-deformed
# fs_LR sphere. Resampling with the files fetched here removes that mismatch. Per subject and
# hemisphere (L, R) the script fetches from the S1200 Structural Preprocessed package:
#   MNINonLinear/Native/<sub>.<H>.sphere.MSMAll.native.surf.gii      the native mesh on the fs_LR sphere (MSMAll)
#   MNINonLinear/Native/<sub>.<H>.midthickness.native.surf.gii       native-side area surface
#   MNINonLinear/fsaverage_LR32k/<sub>.<H>.midthickness_MSMAll.32k_fs_LR.surf.gii   32k-side area surface
#
# The bucket s3://hcp-openaccess needs HCP AWS credentials (ConnectomeDB profile -> "Amazon S3
# access"). Either export AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY, or store them once with
#   aws configure --profile hcp
# and pass -p hcp.
#
# Usage:
#   scripts/hcp_msmall_data.sh [-d <output dir>] [-s "<sub> <sub> ..."] [-p <aws profile>] [-n]
#     -d   where to save, one folder per subject (default: datasets/hcp_training/hcp_msmall)
#     -s   subjects to fetch (default: every numeric subject folder in datasets/hcp_training/freesurfer)
#     -p   aws cli profile holding the HCP credentials
#     -n   dry run, print the commands without running them
#
# Run with -n first, then on one subject, then check that the MSMAll sphere has as many vertices as
# the FreeSurfer surface (the summary prints the command).

set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUCKET=s3://hcp-openaccess/HCP_1200
output_dir="$REPO/datasets/hcp_training/hcp_msmall"
subject_ids=""
profile=""
dry_run=0

while getopts d:s:p:n flag; do
    case "${flag}" in
        d) output_dir=$(realpath -m "${OPTARG}");;
        s) subject_ids=${OPTARG};;
        p) profile=${OPTARG};;
        n) dry_run=1;;
        ?) echo "usage: $(basename "$0") [-d output dir] [-s \"<sub> ...\"] [-p aws profile] [-n]" >&2; exit 1;;
    esac
done

if ! command -v aws > /dev/null; then
    echo "aws cli not found; install it or add it to the PATH" >&2
    exit 1
fi
aws_options=()
if [ -n "$profile" ]; then
    aws_options+=(--profile "$profile")
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
    read -r -a subjects <<< "$subject_ids"
else
    subjects=()
    for d in "$REPO"/datasets/hcp_training/freesurfer/[0-9]*; do
        [ -d "$d" ] && subjects+=("$(basename "$d")")
    done
fi

keys_for_subject() {
    local sub=$1 hemi
    for hemi in L R; do
        echo "MNINonLinear/Native/$sub.$hemi.sphere.MSMAll.native.surf.gii"
        echo "MNINonLinear/Native/$sub.$hemi.midthickness.native.surf.gii"
        echo "MNINonLinear/fsaverage_LR32k/$sub.$hemi.midthickness_MSMAll.32k_fs_LR.surf.gii"
    done
}

# Nothing to fetch when every file is already there, so no credentials are needed in that case
missing=0
for sub in "${subjects[@]}"; do
    for key in $(keys_for_subject "$sub"); do
        [ -s "$output_dir/$sub/$(basename "$key")" ] || missing=$((missing + 1))
    done
done
if [ "$missing" -eq 0 ]; then
    echo "All $((${#subjects[@]} * 6)) files of ${#subjects[@]} subjects are already in $output_dir, nothing to download."
    exit 0
fi

# Check the credentials once, so that a missing key fails here and not once per file
if [ "$dry_run" -eq 0 ]; then
    if ! aws "${aws_options[@]}" s3 ls "$BUCKET/${subjects[0]}/MNINonLinear/Native/" > /dev/null 2>&1; then
        echo "cannot list $BUCKET/${subjects[0]}/MNINonLinear/Native/: check the HCP AWS credentials (see the header of this script)" >&2
        exit 1
    fi
fi

downloaded=0
skipped=0
for sub in "${subjects[@]}"; do
    echo "Fetching MSMAll files for $sub..."
    subject_dir="$output_dir/$sub"
    [ "$dry_run" -eq 1 ] || mkdir -p "$subject_dir"
    for key in $(keys_for_subject "$sub"); do
        {
            target="$subject_dir/$(basename "$key")"
            if [ -s "$target" ]; then
                skipped=$((skipped + 1))
                continue
            fi
            # Download to a temporary name and rename on success, so an interrupted transfer is
            # not mistaken for a complete file on the next run
            if run aws "${aws_options[@]}" s3 cp --only-show-errors "$BUCKET/$sub/$key" "$target.part"; then
                if [ "$dry_run" -eq 0 ]; then
                    mv "$target.part" "$target"
                fi
                downloaded=$((downloaded + 1))
            else
                rm -f "$target.part"
            fi
        }
    done
done

echo "Done: ${#subjects[@]} subjects, $downloaded files downloaded, $skipped already present, $failures failures."
echo "Output: $output_dir/<sub>/"
echo "Verify that the MSMAll sphere and the FreeSurfer surface are the same mesh (vertex counts must match), e.g."
echo "  python -c \"import nibabel as nib; from nilearn import surface; s = '${subjects[0]}'; print(nib.load('$output_dir/' + s + '/' + s + '.L.sphere.MSMAll.native.surf.gii').darrays[0].data.shape[0], surface.load_surf_mesh('$REPO/datasets/hcp_training/freesurfer/' + s + '/surf/lh.white')[0].shape[0])\""
exit $((failures > 0))
