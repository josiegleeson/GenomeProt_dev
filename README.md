# GenomeProt: an integrated proteogenomics analysis platform for long-read RNA-Seq datasets

### Run the shiny application with Docker
Make sure you have [Docker](https://docs.docker.com/engine/install/) installed and the application running in the background before you begin.

Open your terminal application and run:
```
docker run --rm -p 3838:3838 josieg/GenomeProt:latest
```
This will take approximately 10-20 minutes to download the Docker image the first time the app is run.
The --rm removes the container after it’s stopped and the -p 3838:3838 maps your local port 3838 to the same port inside the container.

To **access the local shiny application**, navigate to this link on your web browser http://0.0.0.0:3838.

You can now upload all files and run the steps in your web browser. Although the app is running through a web browser, no files are being uploaded to the internet and everything will be run locally.

To stop the container, close the web browser tab and head back to the terminal where Docker is running and press ctrl+c.


