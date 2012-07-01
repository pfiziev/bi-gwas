INPUT=random_GWAS.json
python convert_to_debi_format.py ../$INPUT
src/debi ./$INPUT.debi results/ -o0.5 -pu
python check_results.py ../$INPUT results/$INPUT.debi.biclusters
