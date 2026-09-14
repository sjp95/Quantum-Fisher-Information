#!/bin/bash

# Check if file exists and is empty
# if [ -s "$1" ]; then
#     echo "File is not empty"
# else
#     echo "File is empty. Deleting..."
#     rm "$1"
# fi

#!/bin/bash


# for kpp in {1..30}
#     do
#        for kpp1 in {1..30}
#        do
#             file_path="SZiSZj_30_100_100_"$kpp"_"$kpp1".dat"

#             # Check if the file exists
#             if [ -f "$file_path" ]; then
#                 # Get the number of lines in the file
#                 num_lines=$(wc -l < "$file_path")

#                 # Check if the file has more than 1 line (ignoring the first line)
#                 if [ "$num_lines" -gt 1 ]; then
#                    echo "File is not empty, keeping it."
#                 else
#                     # File is empty (excluding the first line), so delete it
#                     echo "File is empty, deleting it."
#                     echo $kpp $kpp1
#                     rm "$file_path"
#                 fi
#             else
#                 echo "File does not exist."
#             fi
#         done
# done


#!/bin/bash

for kpp in {1..30}; do
    for kpp1 in {1..30}; do
        file_path="SZiSZj_30_100_100_${kpp}_${kpp1}.dat"

        # Check if the file exists
        if [ -f "$file_path" ]; then
            # Get the number of lines in the file
            num_lines=$(wc -l < "$file_path")

            # Check if the file has more than 1 line (ignoring the first line)
            if [ "$num_lines" -le 1 ]; then
                # File is empty (excluding the first line), so delete it
                echo "File is empty, deleting it: $file_path"
                rm "$file_path"
            fi
        else
            echo "File does not exist: $file_path"
        fi
    done
done
