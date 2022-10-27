FROM rocker/shiny-verse:latest

RUN apt-get update && apt-get install -y \
  libxml2-dev libhdf5-dev \
  libcurl4-openssl-dev \
  libgeos-dev \
  libssl-dev \
  libfontconfig1-dev

RUN apt-get update && apt-get install -y \
  less \
  htop \
  aptitude \
  strace

RUN R -e 'install.packages(c("xml2", "openssl", "rjson"))'
RUN R -e 'install.packages(c("SeuratObject", "Seurat"))'
RUN R -e 'install.packages(c("shinydashboard"))'
#RUN R -e 'install.packages(c("devtools"))'
RUN R -e 'install.packages(c("DT"))'
RUN R -e 'install.packages(c("varhandle"))'
RUN R -e 'install.packages(c("shinyjs"))'
RUN R -e 'install.packages(c("rlist"))'
RUN R -e 'install.packages(c("shinythemes"))'
RUN R -e 'install.packages(c("viridis"))'
RUN R -e 'install.packages(c("logging"))'
RUN R -e 'install.packages(c("svglite"))'
RUN R -e 'install.packages(c("systemfonts"))'

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libigraph0-dev \
  libigraph0v5

RUN ln -sf /usr/share/zoneinfo/Pacific/Auckland /etc/localtime

RUN mkdir -p /app/data
RUN mkdir -p /app/logs

WORKDIR /app

ADD . /app/

RUN cp /app/shiny-server.conf /etc/shiny-server/shiny-server.conf
RUN sudo chown shiny:shiny /app/logs

VOLUME /app/data
VOLUME /app/logs

EXPOSE 3838
