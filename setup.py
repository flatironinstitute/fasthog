from skbuild import setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='fasthog',
      version='0.1.6',
      description='Reasonably fast implementation of HOG descriptor calculation',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Robert Blackwell',
      author_email='rblackwell@flatironinstitute.org',
      url='https://github.com/flatironinstitute/fasthog',
      packages=['fasthog'],
      package_dir={'fasthog': 'fasthog'},
      install_requires=['numpy'],
      cmake_args=['-DCMAKE_BUILD_TYPE=Release'],
      )
