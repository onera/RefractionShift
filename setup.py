import setuptools

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    longDescription = f.read()

setuptools.setup(
    name="RefractionShift", 
    version="0.1.0",
    author="Hanae Labriji",
    author_email="hanae.labriji@gmail.com",
    description="A python module for computing the lateral shift due to atmospheric refraction ",
    long_description=longDescription,
    long_description_content_type="text/markdown",
    url="https://github.com/OneraHub/RefractionShift",
    project_urls={
        "Bug Tracker": "https://github.com/OneraHub/RefractionShift/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6")