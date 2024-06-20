import setuptools

setuptools.setup(
    name="melodia_py",
    version="0.1.2",
    url="https://github.com/rwmontalvao/Melodia_py",
    author="Rinaldo Wander Montalvão, PhD",
    author_email="rwmontalvao@gmail.com",
    description="Differential Geometry of Proteins Backbones",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    install_requires=['biopython>=1.83',
                      'pandas>=2.2.2',
                      'typing-extensions>=4.12',
                      'seaborn>=0.13.2',
                      'jupyterlab>=4.2.2',
                      'tqdm>=4.66.4',
                      'dill>=0.3.8',
                      'nglview>=3.1.2',
                      'pyarrow>=16.1.0',
                      'scikit-learn>=1.5.0',
                      'sty>=1.0.0'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.12',
    ],
    include_package_data=True,
)
