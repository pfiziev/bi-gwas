#INPUT=random_GWAS6.json
INPUT=CEU_GWAS.json

echo $(date): python convert_to_bimax_format.py ../$INPUT
python convert_to_bimax_format.py ../$INPUT

echo "$(date): ./bimax ./$INPUT.bimax > $INPUT.bimax.out"
./bimax ./$INPUT.bimax > $INPUT.bimax.out

echo $(date): python check_results.py ../$INPUT $INPUT.bimax.out
python check_results.py ../$INPUT $INPUT.bimax.out
