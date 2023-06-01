from setuptools import setup, find_packages

setup(
    name='blink',
    author='Abhishek Ganga, Daniyal Rahman, and Jung Tzen Liew',
    description='This tool, Blink, is used to perform genome-wide association studies (GWAS). GWAS is important as it helps scientists identify genes associated with diseases or traits. Our implementation will compared to a GWAS tool called plink which is well-known in the bioinformatic field.',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'blink = blink.blink:main',
        ],
    },
)
