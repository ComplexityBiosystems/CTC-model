#! /bin/bash

# execute all notebooks 
for notebook in notebooks/*.ipynb
do
  echo "# Executing $notebook"
  jupyter nbconvert --ExecutePreprocessor.timeout=3600 --to notebook --execute $notebook  
done
# clean up
echo "# Cleaning up nbconvert files"
rm -f notebooks/*nbconvert* 
echo "# Done!"
