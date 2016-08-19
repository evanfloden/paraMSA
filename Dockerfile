# Dockerfile for lncRNA-Annotation-nf
# Based on debian
############################################################

# Set the base image to Ubuntu
FROM debian:jessie

# File Author / Maintainer
MAINTAINER Evan Floden <evanfloden@gmail.com>

# Install compiler and perl stuff
RUN apt-get update && apt-get install --yes --no-install-recommends \
 wget \
 ed \
 less \
 locales \
 vim-tiny \
 git \
 cmake \
 build-essential \
 gcc-multilib \
 apt-utils \
 perl \
 python \
 ruby \
 expat \
 libexpat-dev \
 libarchive-zip-perl \
 libdbd-mysql \
 libdbd-mysql-perl \
 libdbd-pgsql \
 libgd-gd2-perl \
 libgd2-noxpm-dev \
 libpixman-1-0 \
 libpixman-1-dev \
 graphviz \
 libxml-parser-perl \
 libsoap-lite-perl \
 libxml-libxml-perl \
 libxml-dom-xpath-perl \
 libxml-libxml-simple-perl \
 libxml-dom-perl \
 cpanminus \
 && rm -rf /var/lib/apt/lists/*


# Install perl modules
RUN cpanm --force CPAN::Meta \
 XML::Parser \
 readline \ 
 Term::ReadKey \
 YAML \
 Digest::SHA \
 Module::Build \
 ExtUtils::MakeMaker \
 Test::More \
 Data::Stag \
 Config::Simple \
 Statistics::Lite \
 Statistics::Descriptive \
 Parallel::ForkManager \
 GD \
 GD::Graph \
 GD::Graph::smoothlines \
 Test::Most \
 Algorithm::Munkres \
 Array::Compare Clone \
 PostScript::TextBlock \
 SVG \
 SVG::Graph \
 Set::Scalar \
 Sort::Naturally \
 Graph \
 GraphViz \
 HTML::TableExtract \
 Convert::Binary::C \
 Math::Random \
 Error \
 Spreadsheet::ParseExcel \
 XML::Parser::PerlSAX \
 XML::SAX::Writer \
 XML::Twig XML::Writer \
 && rm -rf /root/.cpanm/work


# Install BioPerl last built and Guidance v2.01
RUN cpanm -v  \
 CJFIELDS/BioPerl-1.6.924.tar.gz \
 && apt-get update -y --fix-missing && apt-get install -y gengetopt \
 libpng-dev \
 time \
 && rm -rf /root/.cpanm/work \
 && wget -q -O- http://guidance.tau.ac.il/ver2/guidance.v2.01.tar.gz | tar xz && cd guidance.v2.01 && make

ENV PERL5LIB="/guidance.v2.01/www/Selecton:/guidance.v2.01/www/bioSequence_scripts_and_constants:/guidance.v2.01/www/Guidance"

#
# Clustal2 
# 
RUN wget -q http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz -O- \
  | tar xz \
  && mv clustalw-*/clustalw2 /usr/local/bin \
  && rm -rf clustalw-*
  
#
# Prank 
# 
RUN wget -O- -q http://wasabiapp.org/download/prank/prank.source.150803.tgz | tar xz \
  && cd prank-msa/src/ \
  && make \
  && cp prank /usr/local/bin \
  && rm -rf /prank-msa
  
#
# Algorithm::Cluster
#
RUN cpanm -n Algorithm::Cluster  

#
# FastTree 
# 
RUN wget -q http://www.microbesonline.org/fasttree/FastTree \
  && chmod +x FastTree \
  && cp ./FastTree /usr/local/bin/
  
#
# HH-Suite 
# 
RUN wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/releases/hhsuite-latest.tar.gz -qO- | tar xz \
  && cd hhsuite-* \
  && make && make install
    
ENV PATH=$PATH:/hhsuite-2.0.16/bin:/hhsuite-2.0.16/scripts     
ENV HHLIB="/hhsuite-2.0.16"
ENV USER=root
  
#
# Download T-Coffee pre-built package
#
RUN wget -q http://www.tcoffee.org/Packages/Stable/Version_11.00.8cbe486/linux/T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz && \
  tar xf T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz && \
  mv T-COFFEE_installer_Version_11.00.8cbe486_linux_x64 /opt/tcoffee && \
  rm -rf T-COFFEE_installer_Version_11.00.8cbe486_linux_x64.tar.gz
    
ENV PATH=$PATH:/opt/tcoffee

#  
RUN wget -q http://mafft.cbrc.jp/alignment/software/mafft-7.130-with-extensions-src.tgz -O- \
  | tar xz && \
  cd mafft-7.130-with-extensions/core && \
  make && make install && \
  cd $HOME && \
  rm -rf mafft-7.130-with-extensions
