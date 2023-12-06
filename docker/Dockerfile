FROM centos:centos7.9.2009 

ENV PATH=/bin:/usr/bin:/sbin:/usr/sbin:/usr/local/bin:/usr/local/sbin:/opt/fsl/bin:/opt/fix:/opt/miniconda/bin:/opt/miniconda/envs/default/bin \
    LD_LIBRARY_PATH=/lib64:/usr/lib64:/usr/local/lib64 \
    FSLTCLSH=/opt/fsl/bin/fsltclsh \
    FSLWISH=/opt/fsl/bin/fslwish \
    FSLDIR=/opt/fsl \
    MCR_INHIBIT_CTF_LOCK=1 \
    FSLOUTPUTTYPE=NIFTI_GZ

RUN echo "Installing OS packages..." \
    && yum --enablerepo=extras install -y epel-release \
    && yum install -y -q bzip2 ca-certificates python-devel curl \
     libcurl libcurl-devel wget unzip bc dc gcc gcc-c++ openssl openssl-devel \
     libxml2 libxml2-devel lapack-3.4.2-8.el7 R R-devel R-java libRmath-devel \
     libjpeg libjpeg-devel git libgit2 libgit2-devel openblas openblas-devel perl-IO-Compress-Zlib perl-TermReadKey \
     perl-Clone perl-Config-General perl-Font-TTF \
     perl-List-MoreUtils perl-Math-Round perl-Text-Format perl-Params-Validate \
     perl-Regexp-Common perl-Statistics-Basic perl-Readonly perl-ExtUtils-MakeMaker cairo-devel libXt-devel \
    && yum install -y centos-release-scl \
    && yum install -y devtoolset-7-gcc-c++ devtoolset-7-gcc-gfortran \
    && source /opt/rh/devtoolset-7/enable bash \
    && yum install -y gmp gmp-devel harfbuzz-devel fribidi-devel udunits2-devel proj proj-devel proj-epsg \
    && wget https://slackware.uk/slackware/slackware64-14.2/slackware64/l/libpng-1.6.23-x86_64-1.txz \
    && tar Jxvf libpng-1.6.23-x86_64-1.txz \
    && ln -s /usr/lib64/libpng16.so.16.23.0 /usr/local/lib/libpng16.so.16 \
    && ln -s /usr/lib64/libpng16.so.16.23.0 /usr/local/lib/libpng16.so \
    && ln -s /usr/lib64/libpng16.la /usr/local/lib/libpng16.la \
    && yum-builddep -y gdal \
    && wget http://download.osgeo.org/gdal/2.4.0/gdal-2.4.0.tar.gz \
    && tar zxvf gdal-2.4.0.tar.gz && cd gdal-2.4.0 && ./configure && make -j8 && make install \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/* \
    && mkdir -p ~/.R \
    && echo 'CXX17 = g++ -std=gnu++17 -fPIC' > ~/.R/Makevars \
    && echo 'CC=/opt/rh/devtoolset-7/root/bin/gcc' >> ~/.R/Makevars \
    && echo 'CXX=/opt/rh/devtoolset-7/root/bin/g++' >> ~/.R/Makevars

RUN echo "Installing R..." \
    && source /opt/rh/devtoolset-7/enable bash \
    && wget https://cdn.rstudio.com/r/centos-7/pkgs/R-3.6.0-1-1.x86_64.rpm \
    && yum -y install ./R-3.6.0-1-1.x86_64.rpm \
    && ln -sf /opt/R/3.6.0/bin/R /usr/bin/R \
    && ln -sf /opt/R/3.6.0/bin/Rscript /usr/bin/Rscript \
    && ln -sf /opt/R/3.6.0/lib/R /usr/lib/R \
    && rm R-3.6.0-1-1.x86_64.rpm \
    && yum clean packages \
    && rm -rf /var/cache/yum/* /tmp/* /var/tmp/* \
    && mkdir -p ~/.R \
    && echo 'CXX17 = g++ -std=gnu++17 -fPIC' > ~/.R/Makevars \
    && echo 'CC=/opt/rh/devtoolset-7/root/bin/gcc' >> ~/.R/Makevars \
    && echo 'CXX=/opt/rh/devtoolset-7/root/bin/g++' >> ~/.R/Makevars \
    && echo "options(Ncpus=8);install.packages(c('devtools', 'remote'), repos='https://cran.rstudio.com');" \
    "r = getOption(\"repos\");r[\"CRAN\"] = \"https://cran.rstudio.com\";options(repos = r);" \
    "require(devtools);install_version(\"kernlab\", version=\"0.9-24\");" \
    "install_version(\"ROCR\", version=\"1.0-7\");" \
    "install_version(\"class\", version=\"7.3-14\",upgrade=FALSE);" \
    "install_version(\"mvtnorm\", version=\"1.0.8\",upgrade=FALSE);" \
    "install_version(\"multcomp\", version=\"1.4-8\",upgrade=FALSE);" \
    "install_version(\"coin\",version=\"1.2.2\",upgrade=FALSE);" \
    "install_version(\"party\", version=\"1.0-25\",upgrade=FALSE);" \
    "install_version(\"e1071\", version=\"1.6-7\",upgrade=FALSE);" \
    "install_version(\"randomForest\", version=\"4.6-12\",upgrade=FALSE);" | R --no-save


RUN echo "Install MCR 2019b..." \
    && curl -sSL -o MCR_R2019b_glnxa64_installer.zip https://ssd.mathworks.com/supportfiles/downloads/R2019b/Release/9/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019b_Update_9_glnxa64.zip \
    && unzip -d MCR_R2019b MCR_R2019b_glnxa64_installer.zip \
    && rm -f MCR_R2019b*.zip \
    && mkdir /opt/MCR2019b \
    && cd MCR_R2019b \
    && echo "destinationFolder=/opt/MCR2019b" > installer_input.txt \
    && echo "agreeToLicense=yes" >> installer_input.txt \
    && echo "mode=silent" >> installer_input.txt \
    && ./install -inputFile installer_input.txt \
    && cd .. && rm -rf MCR_R2019b


RUN echo "Downloading and installing FSL 6..." \
    && wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.4-sources.tar.gz \
    && tar zxf fsl-6.0.4-sources.tar.gz \
    && cd fsl \
    && ./build || true  \
    && cd .. \
    && mv -f fsl /opt/ && chmod a+x /opt/fsl/bin/* \
    && rm -rf fsl-6.0.4-sources.tar.gz

ADD ./fix4mcrv97/fix.tar.gz /opt/

ADD ./spmscripts/spm12_r7771_mcrv97.tar.gz  /opt/

COPY ./src /opt/

RUN echo "Installing miniconda and pipeline..." \
    && yum install -y -q tcsh graphviz  \
    && curl -sSL -o miniconda.sh http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash miniconda.sh -b -p /opt/miniconda \
    && rm -f miniconda.sh \
    && /opt/miniconda/bin/conda update -y -n base conda \
    && /opt/miniconda/bin/conda config --add channels conda-forge \
    && /opt/miniconda/bin/conda create -y -n default python=3.7 \
    && source /opt/miniconda/bin/activate default \
    && conda clean -y --all \
    && chmod a+x /opt/miniconda/bin/* \
    && pip install -U --no-cache-dir pip \
    && cd /opt \
    && pip install -r requirements.txt \
    && pip install -e .


RUN echo "Updating settings ..." \
    && chmod a+rx /opt/spm12_r7771_mcrv97/*spm12* \
    && chmod -R a+rx /opt/spm12_r7771_mcrv97/spm12_mcr/spm12/spm12_r7771_src4v97/* \
    && chmod -R a+rx /opt/fix/* \
    && sed -i 's/backend      : TkAgg/backend      : Agg/' /opt/miniconda/envs/default/lib/python3.7/site-packages/matplotlib/mpl-data/matplotlibrc \
    && STARTUP=/opt/startup.sh \
    && if [ ! -f "$STARTUP" ]; then \
       echo '#!/usr/bin/env bash' >> $STARTUP; \
       echo 'echo;echo "startup: setting environment"...' >> $STARTUP; \
       echo 'if [ -f /docker_environment ];then' >> $STARTUP; \
       echo '   source /docker_environment'  >> $STARTUP; \
       echo 'fi'  >> $STARTUP; \
       echo 'source /opt/miniconda/bin/activate default' >> $STARTUP; \
       echo 'IFS=";" read -ra COMMANDS <<< "$@"' >> $STARTUP; \
       echo 'for command in "${COMMANDS[@]}";do' >> $STARTUP; \
       echo 'echo "`date` : executing ${command}";' >> $STARTUP; \
       echo '  ${command}' >> $STARTUP; \
       echo 'done' >> $STARTUP; \
       fi \
    && chmod 755 $STARTUP


ENTRYPOINT ["/opt/startup.sh"]

