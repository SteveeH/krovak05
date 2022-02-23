from setuptools import setup, find_packages


with open("README.md", "r") as f:
    readme_description = f.read()


setup(
    name="krovak05",
    version="0.0.9",
    author="SteveH",
    author_email="steveeh07@gmail.com",
    description="geodetic package for transformation ETRS coordinates to S-JTSK",
    package_dir={"krovak05": "src"},
    packages=["krovak05", "krovak05.difference_tables",
              "krovak05.kvazigeoids"],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    url="https://github.com/SteveeH/krovak05",
    long_description=readme_description,
    long_description_content_type="text/markdown",
    extras_require={
        "dev": [
            "pytest >=7.0",
        ]
    },
)
