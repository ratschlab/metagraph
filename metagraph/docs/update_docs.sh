rm -r build
make html
make html

# don't forget to set up your ~/.ssh/config accordingly
HOST=metagraph-docs

ssh $HOST <<'ENDSSH'
rm -r ~/projects2019-dna_web_search/static/docs
ENDSSH

scp -r build/html $HOST:~/projects2019-dna_web_search/static/docs

ssh $HOST <<'ENDSSH'
cd ~/projects2019-dna_web_search
docker build . -t dnaloc:prod
export $(cat ~/projects2021-metagraph-indices-deployment/lab-machine/.env | xargs); docker-compose -f ~/projects2021-metagraph-indices-deployment/lab-machine/docker-compose.yml up -d dnaloc
ENDSSH
