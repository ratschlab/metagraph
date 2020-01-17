#!/bin/bash
# script that gets excecuted when a cloud instance is started
chmod a+rx $0

echo "Executing script as: $(whoami)"
# set up SSD disk
sudo mkfs.ext4 -F /dev/nvme0n1
sudo mkdir -p /mnt/disks/ssd
sudo mount /dev/nvme0n1 /mnt/disks/ssd
sudo chmod a+w /mnt/disks/ssd
echo UUID=`sudo blkid -s UUID -o value /dev/disk/by-id/google-local-nvme-ssd-0` /mnt/disks/ssd ext4 discard,defaults,nofail 0 2 | sudo tee -a /etc/fstab

sudo -u ddanciu /home/ddanciu/projects2014-metagenome/metagraph/scripts/cloud/start_client.sh