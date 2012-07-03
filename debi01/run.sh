INPUT=random_GWAS.json
#INPUT=CEU_GWAS.json

echo python convert_to_debi_format.py ../$INPUT
python convert_to_debi_format.py ../$INPUT

echo src/debi ./$INPUT.debi results/ -o0.5 -pu -s50
src/debi ./$INPUT.debi results/ -o0.5 -pu -s50

echo python check_results.py ../$INPUT results/$INPUT.debi.biclusters
python check_results.py ../$INPUT results/$INPUT.debi.biclusters
