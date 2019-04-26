from skbuild import setup
import io
import os.path
import re

here = os.path.abspath(os.path.dirname(__file__))


def readfile(arg):
    with io.open(os.path.join(here, arg)) as f:
        return f.read()


setup(
        name="vipster",
        version=re.findall(r'project\(Vipster VERSION ([0-9.]*)',
                           readfile('CMakeLists.txt'))[0],
        author="Sebastian Gsänger",
        url="https://github.com/sgsaenger/vipster",
        description="A pre- and post-processing toolkit "
                    "for atomistic simulations.",
        long_description=readfile('README.md'),
        long_description_content_type="text/markdown",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Programming Language :: C++",
            "Programming Language :: Python :: Implementation :: CPython",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Scientific/Engineering :: Physics",
            ],
        keywords=['chemistry'],
        license="GPL",
        cmake_args=['-DCMAKE_BUILD_TYPE=Release',
                    '-DPYTHON=YES'],
        )
