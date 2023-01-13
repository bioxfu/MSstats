docker run --rm -v $PWD:/data msstats Rscript msstats.R \
 --input=/data/example_data/MaxQuantResults/ \
 --output=/data/msstats_output \
 --anno=/data/example_data/annotation.csv