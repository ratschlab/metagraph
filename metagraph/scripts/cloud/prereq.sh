#!/bin/bash

# checks that a client machine has all the necessary software in order to successfully run metagraph

# displays a message to stderr in red
function echo_err() {
	RED='\033[0;31m'
	NC='\033[0m'
	echo -e "${RED}Error:${NC} $*" 1>&2;
}


# make sure the aspera client is available
# function check_ascp {
# 	if ! [ -x "$(command -v ascp)" ]; then
# 		echo_err "ascp executable not found. Please install by downloading it from" "https://download.asperasoft.com/download/sw/connect/3.9.7/ibm-aspera-connect-3.9.7.175481-linux-g2.12-64.tar.gz"
#     	exit 1
#   	fi
# }

function check_metagraph {
	if ! [ -x "$(command -v metagraph)" ]; then
		echo_err "metagraph executable not found in PATH. Would it be too much to ask to add it?"
    	exit 1
  	fi
}

function check_gsutil {
	if ! [ -x "$(command -v gsutil)" ]; then
		echo_err "gsutil executable not found in PATH. Please install it from https://cloud.google.com/storage/docs/gsutil_install"
    	exit 1
  	fi
}

############ Main Script ###############
#check_ascp
check_metagraph
check_gsutil