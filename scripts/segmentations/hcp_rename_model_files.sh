#!/bin/bash


# list of all subjects
file_names=()
for d in /BULK/LABDATA/openneuro/nyu4christian/HCP/subjects/* ; do
    if [ -d "$d" ]; then
    file_names+=("$(basename "$d")")
  fi
done
#file_names=("100610")

for sub in "${file_names[@]}"; do

    # Ordnerpfad
    path_sub="/BULK/LABDATA/openneuro/nyu4christian/HCP/freesurfer/$sub/deepRetinotopy"

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