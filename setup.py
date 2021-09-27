import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SpaDecon", 
    version="1.1.0",
    author="Kyle Coleman",
    author_email="kylecole@pennmedicine.upenn.edu",
    description="SpaDecon: Integrating gene expression and histology information to perform cell type deconvolution on spatial transcriptomics data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/kylepcoleman87/SpaDecon",
    packages=setuptools.find_packages(),
    install_requires=["keras==2.2.4","pandas==1.2.4","numpy==1.20.1","scipy==1.6.2","scanpy==1.7.0","anndata==0.7.6","sklearn","tensorflow==1.14.0"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

