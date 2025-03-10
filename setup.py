import os
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from pathlib import Path

class CMakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name, sources=[])

class CMakeBuild(build_ext):
    def build_extension(self, ext):
        ext_dir = Path(self.get_ext_fullpath(ext.name)).parent.resolve()
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={ext_dir}",
            "-DPYBIND11_FINDPYTHON=ON",
            "-DCMAKE_BUILD_TYPE=Release",
        ]
        build_args = ["--config", "Release"]

        os.makedirs(self.build_temp, exist_ok=True)
        subprocess.run(["cmake", "../../"] + cmake_args, cwd=self.build_temp, check=True)
        subprocess.run(["cmake", "--build", "."] + build_args, cwd=self.build_temp, check=True)

setup(
    name="SGEN_Py",
    version="0.0.2",
    author="Oleg Senkevich",
    author_email="oleg.r.senkevich@gmail.com",
    description="SGEN (Stochastic Gene Expression in Neurons) is a high-performance library for computing distributions of mRNAs and proteins in arbitrary neurons. SGEN_Py Python wrapper for SGEN.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    packages=["SGEN_Py"],
    ext_modules=[CMakeExtension("SGEN_Py")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    install_requires=[
        "numpy",
        "pybind11",
        "pyvista"
    ],
)
