FROM ghcr.io/bimberlabinternal/cellmembrane:latest


RUN Rscript -e "devtools::install_github(repo = 'eisascience/HISTA', dependencies = T, upgrade = 'always')" \
&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds


                
# select port
EXPOSE 3838
