import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("src/streamseqs/__init__.py", "r") as fh:
    for line in fh:
        if line.startswith("__version__"):
            version = line.split("=")[1].strip()[1:-1].strip()
            break
    else:
        raise ValueError("No __version__ found.")


setuptools.setup(
    name = "streamseqs",
    version = version,
    description = "Read sequences in (compressed) FASTX files or from stdin",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/osthomas/streamseqs",
    project_urls = {
        "Bug Tracker": "https://github.com/osthomas/streamseqs/issues",
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_dir = {"": "src"},
    packages = setuptools.find_packages(where = "src"),
    python_requires = ">=3.0",
    install_requires = []
)
