from setuptools import setup, find_packages

setup(
    name="Rocketry_formulas",                  # Package name
    version="0.1.0",                    # Version
    author="Sebastian Król",
    # author_email="you@example.com",
    description="a set of commonly used formula during the project ABORT",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Clopeq/Rocketry_formulas.git",
    packages=find_packages(),           # Automatically finds packages
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[                  # Dependencies
        "numpy",
        "scipy",
        "csv"
    ],
)