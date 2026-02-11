## How can I access this app?
The easiest way is to check it out on Rshiny connect: [https://hopetownsend.shinyapps.io/hgsoc-gene-explorer/](https://hopetownsend.shinyapps.io/hgsoc-gene-explorer/).
WARNINGS:
* It runs best on Google chrome (sometimes graphs don't show up on Safari)
* Connect allows limited free monthly runtime so if you hit any problems, please follow the instructions below to run it locally (i.e. on your own computer).

## How can I run this app on my own computer?

1. Download [docker](https://www.docker.com/get-started/) if you don't have it downloaded already.
  * Troubleshooting Note: Docker version 28.4.0, build d8eb465 was used to build the current image.

### Prepare the image
#### Download from Dockerhub
1. Pull the most recent versioned image from Dockerhub (e.g. v2): `docker pull hope2925/gene-explorer`

#### OR Rebuild the image again
1. Clone this repo: `git clone https://github.com/Hope2925/NanoString-RNAseq-HGSOC.git` and go into the appropriate directory
`cd NanoString-RNAseq-HGSOC/hgsoc-gene-explorer/`
2. Build the image with docker: `docker build -t gene-explorer .`

### Run the image
1. Run the shiny app: 
  * If downloaded from Dockerhub: `docker run -p 3838:3838 hope2925/gene-explorer`
  * If built yourself: `docker run -p 3838:3838 gene-explorer`
2. Open the app in your browser: http://localhost:3838
3. Stop the app when you are finished with it: `docker stop $(docker ps -q)`
