FROM rocker/r-base:4.1.2

RUN cat /etc/apt/sources.list
RUN rm /etc/apt/sources.list.d/debian-unstable.list

RUN apt-get update && apt-get install -y \
  strace \
  libxml2-dev libhdf5-dev \
  libcurl4-openssl-dev \
  libgeos-dev \
  libssl-dev \
  libfontconfig1-dev

RUN R -e 'install.packages(c("xml2"))'
RUN R -e 'install.packages(c("openssl"))'
RUN R -e 'install.packages(c("rjson"))'
RUN R -e 'install.packages(c("SeuratObject"))'
RUN R -e 'install.packages(c("Seurat"))'
RUN R -e 'install.packages(c("shiny"))'
RUN R -e 'install.packages(c("shinyjs"))'
RUN R -e 'install.packages(c("shinydashboard"))'
RUN R -e 'install.packages(c("tidyverse"))'
RUN R -e 'install.packages(c("devtools"))'
RUN R -e 'install.packages(c("DT"))'
RUN R -e 'install.packages(c("varhandle"))'
RUN R -e 'install.packages(c("plyr"))'
RUN R -e 'install.packages(c("shinyjs"))'
RUN R -e 'install.packages(c("rlist"))'
RUN R -e 'install.packages(c("shinythemes"))'
RUN R -e 'install.packages(c("viridis"))'
RUN R -e 'install.packages(c("logging"))'
RUN R -e 'install.packages(c("patchwork"))'
RUN R -e 'install.packages(c("svglite"))'
RUN R -e 'install.packages(c("systemfonts"))'
RUN R -e 'install.packages(c("Seurat"))'

RUN ln -sf /usr/share/zoneinfo/Pacific/Auckland /etc/localtime

RUN mkdir -p /app/data

WORKDIR /app

VOLUME /app/data

ADD . /app/

EXPOSE 4242

CMD R -e "shiny::runApp('./', host='0.0.0.0', port=4242)"
