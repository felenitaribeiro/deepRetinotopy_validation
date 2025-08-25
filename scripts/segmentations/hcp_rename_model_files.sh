#!/bin/bash

# this script is used for renaming the files

subjects_dir=""
while getopts s: flag
do
    case "${flag}" in
        s) subjects_dir=${OPTARG};;
        ?)
           echo "script usage: $(basename "$0") [-s path to subs]" >&2
           exit 1;;
    esac
done
# in my case the subjects_dir was "/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer"

# list of all subjects
file_names=()
for d in $subjects_dir/* ; do
    if [ -d "$d" ]; then
    file_names+=("$(basename "$d")")
  fi
done



for sub in "${file_names[@]}"; do

    # path to subjects
    path_sub=""$subjects_dir"/"$sub"/deepRetinotopy"

    # Durchlaufe alle Dateien im Ordner
    for file in "$path_sub"/*; do
        # Nur reguläre Dateien bearbeiten
        if [[ -f "$file" && "$file" == *model* ]]; then
            filename=$(basename "$file")

            # Entferne nur den ersten Teil '_model' (ohne weiteren Unterstrich)
            newname=$(echo "$filename" | sed 's/_model//')

            # Neuer Pfad
            newpath="$path_sub/$newname"

            echo "Renaming: $filename → $newname"
            mv "$file" "$newpath"
        fi
    done
done