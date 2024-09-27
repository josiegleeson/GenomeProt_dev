# GenomeProt: an integrated proteogenomics analysis platform for long-read RNA-Seq datasets

### Run the shiny application with Docker (recommended)
Make sure you have [Docker](https://docs.docker.com/engine/install/) installed and the application running in the background before you begin.

Open your terminal application and run:
```
docker run --rm -p 3838:3838 josieg/genomeprot:dev
```
This will take approximately 10-20 minutes to download the Docker image the first time the app is run.
The --rm removes the container after it’s stopped and the -p 3838:3838 maps your local port 3838 to the same port inside the container.

To **access the local shiny application**, navigate to this link on your web browser http://0.0.0.0:3838.

You can now upload all files and run the steps in your web browser. Although the app is running through a web browser, no files are being uploaded to the internet and everything will be run locally.

To stop the container, close the web browser tab and head back to the terminal where Docker is running and press ctrl+c.

### Locally install the shiny application

The application has substantial dependencies (R requirements found in global.R).

Python (tested with 3.9):
- biopython==1.77
- py-cdhit
- peptides

Command line:
- minimap2
- samtools
- cd-hit
- gffread
- gffcompare

Clone this repository:
```
git clone https://github.com/josiegleeson/GenomeProt_dev.git
```

Obtain the uniprot+openport reference file for your organism and save it into the GenomeProt/data directory.
```
GenomeProt/data/openprot_uniprotDb_hs.txt
```

Run the app from the command line:
```
Rscript -e "shiny::runApp('path/to/app/GenomeProt/', host='0.0.0.0', port=3838)"
```

Otherwise open the app (server.R, ui.R) in RStudio and click 'Run app'.

### Troubleshooting a local installation
The app executes command line tools. Sometimes the local app can't find these tools and it might be easier to install them all into a conda environemnt.
If you do this, you'll need to open the server.R file and replace all instances of 'genomeprot_env' with your conda env name, uncomment out these lines, and comment out the commands which don't include the conda command. You should be able to now run the app as usual, except it will execute command line arguments using your specified conda env. 




