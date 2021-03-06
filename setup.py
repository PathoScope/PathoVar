from setuptools import setup, find_packages

setup(
    name='PathoVar',
    version='0.4.1',
    description='Pathogenic microorganism variant calling and annotation pipeline',
    author='J. Klein; D. Lancour; C. Wake',
    author_email='jaklein@bu.edu',
    url='https://bitbucket.org/pathovar/pathovar',
    packages=find_packages(),
    include_package_data=True,
    data_files=[('pathovar/visualize', ['pathovar/visualize/template.html'])],
    entry_points={
        "console_scripts": [
            "pathovar-config = pathovar.setup_external_data.__main__:main",
            "pathovar-call = pathovar.snp_caller.__main__:main",
            "pathovar-annotate = pathovar.snp_annotation.__main__:main",
            "pathovar = pathovar.__main__:main"
            ]
        },
    zip_safe=False,
    install_requires=[
      "PyVCF >= 0.6.4",
      "beautifulsoup4 >= 4.3.2",
      "requests >= 2.1.0"
    ]
)
