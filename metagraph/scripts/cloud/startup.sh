#!/bin/bash
# script that gets excecuted when a cloud instance is started
chmod a+rx $0

echo "Executing script as: $(whoami)"
# set up SSD disk
sudo mdadm --create /dev/md0 --level=0 --raid-devices=8 /dev/nvme0n1 /dev/nvme0n2 /dev/nvme0n3 /dev/nvme0n4 /dev/nvme0n5 /dev/nvme0n6 /dev/nvme0n7 /dev/nvme0n8
sudo mkfs.ext4 -F /dev/md0
sudo mkdir -p /mnt/disks/ssd
sudo mount /dev/md0 /mnt/disks/ssd
sudo chmod a+w /mnt/disks/ssd

# starts the cloud metagraph client
function start_client() {
  echo "Executing script as: $(whoami)"
  metagraph_path="$HOME/projects2014-metagenome/metagraph"
  PATH="$PATH:~/.aspera/connect/bin/:${metagraph_path}/build:/home/linuxbrew/.linuxbrew/bin/:/snap/bin:/usr/local/ncbi/sra-tools/bin/:${metagraph_path}/build/KMC/"

  cd "$metagraph_path"
  git fetch origin dd/cloud
  git checkout origin/dd/cloud  # TODO: remove this at some point

  mkdir -p build
  cd build
  cmake ..
  make -j metagraph

  cd "$metagraph_path/scripts/cloud"
  python3 ./client2.py --output_dir="/mnt/disks/ssd/"
}

export -f start_client
su ddanciu -c "bash -c start_client"