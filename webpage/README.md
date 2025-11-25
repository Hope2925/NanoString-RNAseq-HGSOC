## How can I run this app on my own computer?

1. Download [docker](https://www.docker.com/get-started/) if you don't have it downloaded already.
  * Troubleshooting Note: Docker version 28.4.0, build d8eb465 was used to build the current image.

### Download the premade image
1. Pull the latest image from Github: `docker pull ghcr.io/<your-username>/gene-explorer:latest`
2. Run the shiny app: `docker run -d -p 3838:3838 ghcr.io/<your-username>/gene-explorer:latest`
3. Open the app in your browser: http://localhost:3838
4. Stop the app when you are finished with it: `docker stop $(docker ps -q)`


### Build the image again
1. Clone this repo: `git clone https://github.com/Hope2925/<repo-name>.git` and go into the appropriate directory
`cd <repo-name>/webpage/`
2. Build the image with docker: `docker build -t gene-explorer .`
3. Run it: `docker run -p 3838:3838 gene-explorer`
4. Open the app in your browser: http://localhost:3838
5. To end it:
  A. Run `docker ps -a` where you should see something like this
```
  CONTAINER ID   IMAGE               COMMAND                  CREATED         STATUS                       PORTS                                         NAMES
19f13e7b6c84   gene-explorer       "Rscript app.R"          6 minutes ago   Up 6 minutes                 0.0.0.0:3838->3838/tcp, [::]:3838->3838/tcp   eloquent_babbage
```

  Then take the CONTAINER ID and add it to `docker stop`. So in the above case I would run `docker stop 19f13e7b6c84`