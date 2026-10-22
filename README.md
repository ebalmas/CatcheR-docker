# CatcheR-docker

This repository describes all the Docker images available for [catcheR](https://github.com/ebalmas/catcheR)

hedgelab/ciriseqv4-hedgelab_1.1:base is our current image in the repository

Potentially, the Docker image can be rebuilt by changing the name of the image:

```bash
docker build --no-cache -t hedgelab/ciriseqv4-hedgelab_1.1:base .
docker push hedgelab/ciriseqv4-hedgelab_1.1:base 
```
Load the image with our script.
You can place it into your working/experiment folder.

Ideally, organize the experimental folder and data with:

Output folder
Scratch (input files)
scripts
start-ciri.sh

When the docker is lounched these folders will be created if not present.
/home/project 
will contain your shared experiment folder (shared within your computer and the Docker)

Load the Docker:

```bash
cd path/to/your/folder
sh start-ciri.sh .

Scegli la password di RStudio: #write your password when the prompt appears
RStudio: http://localhost:8081
Utente: rstudio
```
Open the browser and type: http://localhost:8081 you will be able to log in into the complete RStudio

Log in by entering username: rstudio and password: your password

If you use a different image, you should change the name inside the script start-ciri.sh




