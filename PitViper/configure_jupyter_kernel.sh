conda env create -f pitviper_jupyter.yml -n pitviper_jupyter -y

conda activate pitviper_jupyter

pip install --user ipykernel

python -m ipykernel install --user --name=pitviper_jupyter
