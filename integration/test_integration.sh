# import example calculations
nomad upload --unstage tests/data/proc/examles_vasp.zip

# create a new index with coe repoTool
docker run nomad/coe-repotool:latest

# start coe repoServer
docker run nomad/coe-repowebservice:latest

# try to search for new calculations
curl localhost:8111/...