import setuptools

with open("README.md", "r", encoding="utf-8") as fp:
    long_description = fp.read()

setuptools.setup(
    name="c14-witcher-PiscatorX",
    version="0.0.1a1",
    author="Andrew Ndhlovu",
    author_email="drewxdvst@outlook.com",
    description="A tool for in silico identification and classification of caspase-like homologs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/caspase-like-homolog-identifier/c14_witcher.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
