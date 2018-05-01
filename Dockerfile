# Based on Ubuntu 16.04 LTS release for stability; supported until 2021/04
FROM ubuntu:16.04

# Update package repositories, install git,
# fetch FiNGS from GitHub and install all dependencies
RUN apt-get update && \
	apt-get -y install git && \
	git clone https://github.com/cpwardell/FiNGS.git && \
	cd FiNGS && \
	chmod +x install_FiNGS_dependencies.sh && \
	./install_FiNGS_dependencies.sh

