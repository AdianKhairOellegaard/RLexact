#!/bin/bash

# The exact name of your template file
TEMPLATE="TFIMQFT_10x1x1.h"

# 1. Safety Check: Abort if the file isn't in the exact same folder
if [ ! -f "$TEMPLATE" ]; then
    echo "Error: Could not find '$TEMPLATE'."
    echo "Make sure the script is running in the exact same directory as your input file."
    exit 1
fi

# The N values for your replicas
N_VALUES=(12 14 16 18 20 22 24 26 28 30)

for N in "${N_VALUES[@]}"; do
    OUTPUT="TFIMQFT_${N}x1x1.h"
    echo "Creating $OUTPUT..."

    # Generate the sequences for the custom symmetries
    sym1="Custom symmetry $(seq -s " " 0 $((N-1)))"
    sym2="Custom symmetry $(seq -s " " 1 $((N-1))) 0"

    # Clear the output file before writing
    > "$OUTPUT"

    # 2. Read the template line by line
    while IFS= read -r line || [ -n "$line" ]; do
        
        # Strip hidden Windows carriage returns (\r)
        line="${line%$'\r'}"

        # Specific check first so it doesn't get caught by the generic one below
        if [[ "$line" == "Number of spins in unit cell "* ]]; then
            echo "$line" >> "$OUTPUT"
        elif [[ "$line" == "Number of spins "* ]]; then
            echo "Number of spins $N" >> "$OUTPUT"
        elif [[ "$line" == "Number of couplings "* ]]; then
            echo "Number of couplings $N" >> "$OUTPUT"
        elif [[ "$line" == "Custom symmetry 0"* ]]; then
            echo "$sym1" >> "$OUTPUT"
        elif [[ "$line" == "Custom symmetry 1"* ]]; then
            echo "$sym2" >> "$OUTPUT"
        elif [[ "$line" == "Coupling vector "* ]]; then
            continue
        else
            echo "$line" >> "$OUTPUT"
        fi
    done < "$TEMPLATE"

    # 3. Append the new coupling vectors for the current N
    for (( i=0; i<N; i++ )); do
        j=$(( (i + 1) % N ))
        echo "Coupling vector $i $j 0" >> "$OUTPUT"
    done

done

echo "All replicas generated successfully with intact input data!"