### Install with anaconda with an environment

To install the code in your computer, you need first to install the anaconda
python (https://www.anaconda.com/download).

We advise to install the latest anaconda distribution.

conda create -n sospex python=3.10

conda activate sospex

conda install -c conda-forge -y lmfit  reproject  fitsio

conda install -c darioflute -y sospex

To exit the environment, type:

conda deactivate

### Note for windows users

If fitsio is not available for windows, sospex can be still used since it defaults to open fits is the astropy library.

### Updating with anaconda

If there is a new release, you will have only to update the package:

conda update -c darioflute sospex -y

Alternatively, add the channel to your anaconda:

conda config --add channels darioflute

and simply do:

conda update sospex

### Install with pip -e:

An alternative is to install the package with pip after cloning it.

#### A) clone the package

git clone https://github.com/darioflute/sospex.git

#### B) add the required libraries

conda install -c conda-forge lmfit reproject fitsio

#### C) pip install

From the directory containing the sospex package:

pip install -e sospex


#### D) to update with the latest release

cd sospex

git pull




