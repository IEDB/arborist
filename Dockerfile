FROM ubuntu@sha256:c9cf959fd83770dfdefd8fb42cfef0761432af36a764c077aed54bbc5bb25368

RUN apt-get update
RUN apt-get install -y curl unzip
RUN apt-get install -y sqlite3 libsqlite3-dev
RUN apt-get install -y mariadb-client
RUN apt-get install -y python3-pip
RUN apt-get install -y openjdk-11-jre-headless
RUN apt-get install -y ncbi-blast+ hmmer

WORKDIR /arborist
ADD requirements.txt /arborist
RUN pip3 install -r requirements.txt

ADD Makefile /arborist
RUN make deps
