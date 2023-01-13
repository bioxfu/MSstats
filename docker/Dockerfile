FROM r-base:4.2.2

COPY sources.list /etc/apt/sources.list

RUN apt-get update

RUN apt-get install -y cmake libcurl4-openssl-dev libssl-dev libxml2-dev libnetcdf-dev

RUN Rscript -e 'install.packages("BiocManager")'

RUN Rscript -e 'BiocManager::install("MSstats")'

RUN Rscript -e 'BiocManager::install("MSstatsTMT")'

RUN Rscript -e 'BiocManager::install("MSstatsPTM")'

RUN Rscript -e 'BiocManager::install("MSstatsLiP")'

RUN Rscript -e 'BiocManager::install("MSstatsLOBD")'

RUN Rscript -e 'BiocManager::install("MSstatsSampleSize")'

RUN Rscript -e 'BiocManager::install("MSstatsQC")'

RUN Rscript -e 'install.packages("docopt")'

COPY scripts/msstats.R /app/msstats.R

WORKDIR /app

CMD ["Rscript", "msstats.R"]


