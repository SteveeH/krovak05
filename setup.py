from setuptools import setup


with open("README.md","r") as f:
    readme_description = f.read()


setup(
    name="krovak05",
    version="0.0.2",
    author="SteveH",
    author_email="steveeh07@gmail.com",
    description="geodetic package for transformation ETRS coordinates to S-JTSK",
    py_modules=["transformation"],
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    url="https://github.com/SteveeH/krovak05",
    long_description=readme_description,
    long_description_content_type="text/markdown",
    extras_require = {
        "dev": [
            "pytest >=7.0",
        ]
    },
)