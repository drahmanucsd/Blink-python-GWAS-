from setuptools import setup, find_packages

# Check if haptools is installed
try:
    import haptools
    py_modules = ['blink']
    console_script = 'blink = blink:main'
except ImportError:
    has_haptools = False
    py_modules = ['blink_wo_haptools']
    console_script = 'blink = blink_wo_haptools:main'

setup(
    name='blink',
    version=1.9,
    py_modules=py_modules,
    author='Abhishek Ganga, Daniyal Rahman, and Jung Tzen Liew',
    description='This tool, Blink, is used to perform genome-wide association studies (GWAS). GWAS is important as it helps scientists identify genes associated with diseases or traits. Our implementation will be compared to a GWAS tool called plink which is well-known in the bioinformatic field.',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            console_script
        ],
    },
)