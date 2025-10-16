import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="curvesketching",
    version="0.1",
    author="Isidro Gonzalez",
    author_email="isidro.gonzalez@uzh.ch",
    description="A simple package for Curve Sketching",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/isoUZH/curvesketching.git",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["curvesketching"],
    python_requires=">=3.0",
)
