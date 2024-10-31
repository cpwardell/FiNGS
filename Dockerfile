# Use the continuumio/miniconda3 base image
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Clone the repository directly
RUN apt-get update && apt-get install -y git && \
    git clone https://github.com/cpwardell/FiNGS.git . && \
    apt-get clean

# Install the package
RUN pip install .

# Set entrypoing and the default command
ENTRYPOINT ["fings"]
CMD ["--help"]
