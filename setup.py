import setuptools
from glob import glob
from os.path import dirname, join

DIR = (dirname(__file__) or '.')

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pangolin",
    version="1.0.2",
    author="Tony Zeng",
    author_email="tkyzeng@gmail.com",
    description="Pangolin",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=['pangolin'],
    package_data={
        "pangolin": ["models/*"],
    },
    entry_points={
        "console_scripts": [
            "pangolin=pangolin.pangolin:main"
        ]
    },
    scripts=glob(join(DIR, 'scripts/*.py')),
)
