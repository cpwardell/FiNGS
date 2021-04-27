# Based on Miniconda3 image
FROM continuumio/miniconda3

# Update conda channels and install FiNGS
RUN conda config --add channels defaults && \
conda config --add channels bioconda && \
conda config --add channels conda-forge && \
conda install -y fings

# Default entrypoint is FiNGS, default command is "--help"
ENTRYPOINT ["fings"]
CMD ["--help"]
